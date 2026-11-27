/*
 * http://sobereva.com/125
 * 所有坐标、平移矢量单位都是bohr而不是埃
 */

use crate::base::Mole;
use crate::datas::elems;
use crate::utils::fort_form;
use crate::writer::Writer;
use chrono::Local;
use std::fs::File;
use std::io::prelude::*;

/// Cub格式文件(.cub)写入器
pub struct CubWriter {
    /// 标题（写入文件第一行）
    pub title: String, // 标题
    /// 原子符号列表
    pub syms: Vec<String>,
    /// 原子坐标
    pub xyzs: Vec<[f64; 3]>,
    /// 轨道索引
    pub obts: Vec<isize>,
    /// 起始坐标
    pub pos0: [f64; 3],
    /// 三个方向的格点数量
    pub size: [usize; 3],
    /// 三个方向格点步长
    pub step: [f64; 3],
    /// 数值: `[轨道, 格点]`
    pub vals: Vec<Vec<f64>>,
}

impl CubWriter {
    /// 创建CubWriter实例
    pub fn new() -> Self {
        Self {
            title: "generate by rswfn".to_string(),
            syms: vec![],
            xyzs: vec![],
            obts: vec![],
            pos0: [0.0; 3],
            size: [0; 3],
            step: [0.0; 3],
            vals: vec![],
        }
    }

    fn write_grids(&mut self, file: &mut File) {
        let now = Local::now();
        let now = now.format("%Y-%m-%d %H:%M:%S");
        let nops = self.size[0] * self.size[1] * self.size[2];
        write!(file, "{}\n{} {}\n", self.title, now, nops * self.obts.len()).unwrap();
        let [nx, ny, nz] = self.size;
        let [x0, y0, z0] = self.pos0;
        let [sx, sy, sz] = self.step;
        let natm = self.syms.len() as isize;
        write!(file, "{:>5}{:>12.6}{:>12.6}{:>12.6}\n", -natm, x0, y0, z0).unwrap();
        write!(file, "{:>5}{:>12.6}{:>12.6}{:>12.6}\n", nx, sx, 0.0, 0.0).unwrap();
        write!(file, "{:>5}{:>12.6}{:>12.6}{:>12.6}\n", ny, 0.0, sy, 0.0).unwrap();
        write!(file, "{:>5}{:>12.6}{:>12.6}{:>12.6}\n", nz, 0.0, 0.0, sz).unwrap();
    }

    fn write_coord(&mut self, file: &mut File) {
        let natm = self.syms.len();
        for i in 0..natm {
            let sym = &self.syms[i];
            let atm = elems::sym2atm(sym);
            let [x, y, z] = self.xyzs[i];
            write!(file, "{:>5}{:>12.6}{:>12.6}{:>12.6}{:>12.6}\n", atm, atm as f64, x, y, z).unwrap();
        }
    }

    fn write_value(&mut self, file: &mut File) {
        let nobt = self.obts.len();
        let npos = self.size[0] * self.size[1] * self.size[2];
        let obts: Vec<isize> = self.obts.iter().map(|obt| *obt + 1).collect();
        let mut info: Vec<isize> = vec![nobt as isize];
        info.extend(obts.clone());
        for (i, obt) in info.iter().enumerate() {
            write!(file, "{:>5}", obt).unwrap();
            if (i + 1) % 10 == 0 {
                write!(file, "\n").unwrap();
            }
        }
        if obts.len() % 10 != 0 {
            write!(file, "\n").unwrap();
        }

        let nz = self.size[2];
        let mut index = 0;
        // 对每一个点进行循环
        for i in 0..npos {
            // 对每一个轨道进行循环
            for j in 0..nobt {
                let val = self.vals[j][i];
                write!(file, "{:>13}", fort_form(val, 5)).unwrap();
                index += 1;
                if index == nz * nobt {
                    write!(file, "\n").unwrap();
                    index = 0;
                    continue;
                }
                if index % 6 == 0 {
                    write!(file, "\n").unwrap();
                }
            }
        }
    }

    /// 将数据保存到文件
    pub fn save(&mut self, path: &str) {
        let mut file = File::create(path).unwrap();
        self.write_grids(&mut file);
        self.write_coord(&mut file);
        self.write_value(&mut file);
    }
}

impl Writer for CubWriter {
    fn read_mole(&mut self, mole: &Mole) {
        self.syms = mole.atoms.syms();
        self.xyzs = mole.atoms.xyzs().clone();
    }
}
