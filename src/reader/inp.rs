// 读取orca的输入文件

use crate::base::{Atoms, Basis, Coefs, Stm};
use regex::Regex;

use crate::datas::elems::sym2atm;
use crate::reader::{Reader, read_lines};
use std::collections::HashMap;

/// Orca输入文件(.inp)读取器
pub struct InpReader {
    /// 文件中各数据块的标记位置(行号)
    pub marks: HashMap<String, usize>,
    /// 文件行
    pub lines: Vec<String>,
}

impl InpReader {
    /// 创建InpReader实例
    pub fn new(path: &str) -> InpReader {
        let mut marks = HashMap::new();
        let lines = read_lines(path);
        for (l, line) in lines.iter().enumerate() {
            if line.starts_with("%coords") {
                marks.entry("geome".to_string()).or_insert(l);
            }
        }
        InpReader { marks, lines }
    }

    /// 读取几何结构
    pub fn read_geome(&self) -> (Vec<String>, Vec<[f64; 3]>) {
        let mark = self.marks.get("geome").expect("获取结构行失败");
        let reg1 = Regex::new(r"([a-z]+) +(-?\d+.\d+) +(-?\d+.\d+) +(-?\d+.\d+)").expect("正则表达式构建失败");
        let mut syms = vec![];
        let mut xyzs = vec![];
        for line in self.lines.iter().skip(*mark) {
            if line.starts_with("end") {
                break;
            }
            if reg1.is_match(line) {
                let caps = reg1.captures(line).expect("正则表达式解析失败");
                let sym: String = caps[1].parse().expect("解析元素失败");
                let x: f64 = caps[2].parse().expect("解析坐标失败");
                let y: f64 = caps[3].parse().expect("解析坐标失败");
                let z: f64 = caps[4].parse().expect("解析坐标失败");
                syms.push(sym);
                xyzs.push([x, y, z]);
            }
        }
        (syms, xyzs)
    }
}

impl Reader for InpReader {
    fn get_atoms(&self) -> Atoms {
        let (syms, xyzs) = self.read_geome();
        let atms: Vec<usize> = syms.iter().map(|sym| sym2atm(&sym.to_uppercase())).collect();
        let crgs: Vec<f64> = atms.iter().map(|atm| *atm as f64).collect();
        Atoms::new(&atms, &xyzs, &crgs)
    }
    fn get_basis(&self) -> Option<Basis> {
        None
    }

    fn get_coefs(&self) -> Option<Coefs> {
        None
    }
}

/// 元素符号从小写到大写的转换(预留)
pub fn sym_s2b() {}
