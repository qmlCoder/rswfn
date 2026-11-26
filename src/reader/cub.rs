/*
 * http://sobereva.com/125
 * 所有坐标、平移矢量单位都是bohr而不是埃
 */
use std::path;

use regex::Regex;

use crate::{
    base::{Atoms, Basis, Coefs},
    props::gridprop::CubeGrid,
    reader::{Reader, read_lines},
};

/// 格点文件(.cub)读取器
pub struct CubReader {
    pub path: String,
    /// 文件行
    pub lines: Vec<String>,
    /// 是否为分子轨道数据
    pub mobt: bool,
    /// 原子数量
    pub natm: usize,
}

impl CubReader {
    /// 创建CubReader实例
    pub fn new(path: &str) -> Self {
        let lines = read_lines(path);
        let path = path.to_string();
        let natm = lines[2][..5].trim().parse::<isize>().unwrap();
        let mobt = if natm < 0 { true } else { false };
        let natm = natm.abs() as usize;
        Self {
            path,
            lines,
            mobt,
            natm,
        }
    }

    /// 读取格点信息
    pub fn read_cube(&self) -> CubeGrid {
        let mut start = [0.0; 3];
        let px = self.lines[2][5..17].trim().parse::<f64>().unwrap();
        let py = self.lines[2][17..29].trim().parse::<f64>().unwrap();
        let pz = self.lines[2][29..41].trim().parse::<f64>().unwrap();
        let start = [px, py, pz];
        let mut shape = [0; 3];
        for i in 0..3 {
            let line = &self.lines[i + 3];
            shape[i] = line[..5].trim().parse::<usize>().unwrap();
        }
        let sx = self.lines[3][5..17].trim().parse::<f64>().unwrap();
        let sy = self.lines[4][17..29].trim().parse::<f64>().unwrap();
        let sz = self.lines[5][29..41].trim().parse::<f64>().unwrap();
        let steps = [sx, sy, sz];
        CubeGrid { shape, start, steps }
    }

    /// 读取原子信息
    pub fn read_atms(&self) -> (Vec<usize>, Vec<f64>, Vec<[f64; 3]>) {
        let natm = self.natm;
        let mut atms = Vec::with_capacity(natm);
        let mut crgs = Vec::with_capacity(natm);
        let mut xyzs = Vec::with_capacity(natm);
        for l in 6..6 + natm {
            let line = &self.lines[l];
            let atm = line[..5].trim().parse::<usize>().unwrap();
            let crg = line[5..17].trim().parse::<f64>().unwrap();
            let x = line[17..29].trim().parse::<f64>().unwrap();
            let y = line[29..41].trim().parse::<f64>().unwrap();
            let z = line[41..53].trim().parse::<f64>().unwrap();
            atms.push(atm);
            crgs.push(crg);
            xyzs.push([x, y, z]);
        }
        (atms, crgs, xyzs)
    }

    /// 读取格点数值
    pub fn read_vals(&self) -> (Vec<usize>, Vec<Vec<f64>>) {
        let (nobt, obts, skip) = match self.lines[6 + self.natm][..5].trim().parse::<usize>() {
            Ok(nobt) => {
                let obts: Vec<usize> = self.lines[6 + self.natm][5..]
                    .trim()
                    .split_whitespace()
                    .map(|e| e.parse::<usize>().unwrap())
                    .collect();
                (nobt, obts, 1)
            }
            Err(_) => (1, vec![0], 0),
        };
        let [nx, ny, nz] = self.read_cube().shape;
        let count = nx * ny * nz;
        let mut vals = vec![vec![0.0; count]; nobt];
        let mut idx = 0;
        let re = Regex::new(r"-?\d+.\d+E[+-]\d").unwrap();
        for line in self.lines.iter().skip(6 + self.natm + skip) {
            // 匹配改行文本中的所有数值
            let line_vals: Vec<f64> = re
                .find_iter(line)
                .map(|cap| cap.as_str().parse::<f64>().unwrap())
                .collect();
            for (i, val) in line_vals.iter().enumerate() {
                vals[idx % nobt][idx / nobt] = *val;
                idx += 1;
            }
        }
        (obts, vals)
    }
}

impl Reader for CubReader {
    fn get_atoms(&self) -> Atoms {
        let (atms, crgs, xyzs) = self.read_atms();
        Atoms::new(&atms, &xyzs, &crgs)
    }

    fn get_basis(&self) -> Option<Basis> {
        None
    }

    fn get_coefs(&self) -> Option<Coefs> {
        None
    }
}
