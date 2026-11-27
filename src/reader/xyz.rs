use std::sync::OnceLock;

use regex::Regex;

use crate::base::{Atoms, Basis, Coefs, Stm};
use crate::datas::consts::BOHR;
use crate::datas::elems::sym2atm;
use crate::reader::{Reader, read_lines};

/// Xyz格式文件(.xyz)读取器
pub struct XyzReader {
    /// 文件路径
    pub path: String,
    /// 缓存的多帧几何结构(惰性解析)
    pub geomes: OnceLock<Vec<Atoms>>,
    /// 文件行
    pub lines: Vec<String>,
}

impl XyzReader {
    /// 创建XyzReader实例
    pub fn new(path: &str) -> Self {
        Self {
            path: path.to_string(),
            lines: read_lines(path),
            geomes: OnceLock::new(),
        }
    }
    /// 读取所有几何结构(支持多个构象)
    pub fn read_geomes(&self) -> &Vec<Atoms> {
        let geomes = self.geomes.get_or_init(|| {
            let r2 = Regex::new(r"^ *([A-Za-z]+) +(-?\d+.\d+) +(-?\d+.\d+) +(-?\d+.\d+)$").expect("reg new error");
            let mut atms = Vec::new();
            let mut xyzs = Vec::new();
            let mut geomes = Vec::new();
            let mut loc = 0;
            for (l, line) in self.lines.iter().enumerate() {
                if r2.is_match(line) {
                    let caps = r2.captures(line).expect("caps error");

                    let sym = caps[1].to_string();
                    let x = caps[2].parse::<f64>().expect("f64 parse error") * BOHR;
                    let y = caps[3].parse::<f64>().expect("f64 parse error") * BOHR;
                    let z = caps[4].parse::<f64>().expect("f64 parse error") * BOHR;
                    atms.push(sym2atm(&sym));
                    xyzs.push([x, y, z]);
                    loc = l;
                } else {
                    if atms.len() != 0 && xyzs.len() != 0 {
                        let crgs: Vec<f64> = atms.iter().map(|atm| *atm as f64).collect();
                        let atoms = Atoms::new(&atms, &xyzs, &crgs);
                        geomes.push(atoms);
                    }
                    atms.clear();
                    xyzs.clear();
                }
            }
            if loc == self.lines.len() - 1 {
                let crgs: Vec<f64> = atms.iter().map(|atm| *atm as f64).collect();
                let atoms = Atoms::new(&atms, &xyzs, &crgs);
                geomes.push(atoms);
            }
            geomes
        });
        geomes
    }
}

impl Reader for XyzReader {
    fn get_atoms(&self) -> Atoms {
        let geomes = self.read_geomes();
        let count = geomes.len();
        assert!(count > 0, "没有找到几何结构");
        geomes[count - 1].clone()
    }
    fn get_basis(&self) -> Option<Basis> {
        None
    }
    fn get_coefs(&self) -> Option<Coefs> {
        None
    }
}
