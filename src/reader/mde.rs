use ndarray::Array2;
use regex::Regex;
use std::collections::HashMap;
use std::f64::consts::PI;
use std::sync::OnceLock;

use crate::base::{Atoms, Basis, Coefs};
use crate::datas::btrans;
use crate::reader::{BasData, Reader, read_lines, shl_basis, shl_forms};

/// Molden格式文件(.molden)读取器
pub struct MdeReader {
    path: String,
    lines: Vec<String>,
    marks: HashMap<String, usize>,
    btype: String,
    read_basis: OnceLock<(
        // 缓存函数执行结果
        Vec<usize>,
        Vec<usize>,
        Vec<String>,
        Vec<BasData>,
    )>,
}

impl MdeReader {
    /// 创建MdeReader实例
    pub fn new(path: &str) -> Self {
        let lines = read_lines(path);
        let mut marks = HashMap::new();
        let mut btype = String::new();
        for (l, line) in lines.iter().enumerate() {
            let line = line.trim();
            if line.starts_with("[Atoms]") {
                *marks.entry("geome".to_string()).or_insert(l as usize) = l as usize;
                continue;
            }
            if line == "[GTO]" {
                *marks.entry("GTO".to_string()).or_insert(l as usize) = l as usize;
                continue;
            }
            if line == "[MO]" {
                *marks.entry("MO".to_string()).or_insert(l as usize) = l as usize;
                continue;
            }
            if line.starts_with("[5D") {
                *marks.entry("sph".to_string()).or_insert(l as usize) = l as usize;
                btype = "sph".to_string();
                continue;
            }
        }
        log::debug!("marks: {marks:?}");
        Self {
            path: path.to_string(),
            lines: lines,
            marks: marks,
            btype: btype,
            read_basis: OnceLock::new(),
        }
    }

    /// 读取几何结构
    pub fn read_geome(&self) -> (Vec<usize>, Vec<[f64; 3]>) {
        let mark = self.marks.get("geome").expect("获取原子坐标标题行");
        let reg1 = Regex::new(r"[A-Za-z] +\d+ +(\d+) +(-?\d+.\d+) +(-?\d+.\d+) +(-?\d+.\d+)").unwrap();
        let mut atms = Vec::new();
        let mut xyzs = Vec::new();
        for line in self.lines.iter().skip(*mark + 1) {
            if reg1.is_match(line) {
                let caps = reg1.captures(line).expect("cap error");
                let atm: usize = caps[1].parse().unwrap();
                let x: f64 = caps[2].parse().unwrap();
                let y: f64 = caps[3].parse().unwrap();
                let z: f64 = caps[4].parse().unwrap();
                atms.push(atm);
                xyzs.push([x, y, z]);
            } else {
                break;
            }
        }
        (atms, xyzs)
    }

    /// 读取基组信息(带缓存)
    pub fn read_basis(&self) -> (&Vec<usize>, &Vec<usize>, &Vec<String>, &Vec<BasData>) {
        let (ato_atms, ato_shls, ato_syms, datas) = self.read_basis.get_or_init(|| {
            let reg1 = Regex::new(r"^ +(\d+)( +\d)?$").expect("正则表达式构建失败");
            let reg2 = Regex::new(r" *([spdf]) +(\d+)( 1.0)? *").expect("正则表达式构建失败");
            let reg3_1 = Regex::new(r"^ +(-?\d+.\d+E[+-]\d+) +(-?\d+.\d+E[+-]\d+)").expect("正则表达式构建失败");
            let reg3_2 = Regex::new(r" +(\d+.\d+) +(-?\d+.\d+)").expect("正则表达式构建失败");
            let mark = self.marks.get("GTO").expect("获取mark GTO 失败");
            let mut atm = 0;
            let mut shl = 0;
            let mut sym = "".to_string();
            let syms = ["s", "p", "d", "f"];
            let mut ato_syms = Vec::new();
            let mut ato_atms = Vec::new();
            let mut ato_shls = Vec::new();
            let mut datas: Vec<BasData> = Vec::new();
            for l in (*mark + 1)..self.lines.len() {
                let line = self.lines.get(l).unwrap();
                if line == "" {
                    continue;
                }
                if reg1.is_match(line) {
                    let caps = reg1.captures(line).expect("cap error");
                    atm = caps[1].parse::<usize>().expect("解析失败") - 1;
                    shl = 0;
                    continue;
                }
                if reg2.is_match(line) {
                    let caps = reg2.captures(line).expect("cap error");
                    sym = caps[1].parse::<String>().expect("解析失败");
                    shl += 1;
                    match (sym.as_str(), self.btype.as_str()) {
                        ("s", _) => {
                            ato_syms.extend(["S"]);
                            ato_atms.extend([atm; 1]);
                            ato_shls.extend([shl - 1; 1]);
                        }
                        ("p", _) => {
                            ato_syms.extend(["PX", "PY", "PZ"]);
                            ato_atms.extend([atm; 3]);
                            ato_shls.extend([shl - 1; 3]);
                        }
                        ("d", "sph") => {
                            ato_syms.extend(btrans::SPH_D_SYMS);
                            ato_atms.extend([atm; 5]);
                            ato_shls.extend([shl - 1; 5]);
                        }
                        ("d", "car") => {
                            ato_syms.extend(btrans::CAR_D_SYMS);
                            ato_atms.extend([atm; 6]);
                            ato_shls.extend([shl - 1; 6]);
                        }
                        ("f", "sph") => {
                            ato_syms.extend(btrans::SPH_F_SYMS);
                            ato_atms.extend([atm; 7]);
                            ato_shls.extend([shl - 1; 7]);
                        }
                        ("f", "car") => {
                            ato_syms.extend(btrans::CAR_F_SYMS);
                            ato_atms.extend([atm; 10]);
                            ato_shls.extend([shl - 1; 10]);
                        }
                        _ => panic!("未知类型"),
                    }

                    continue;
                }
                if reg3_1.is_match(line) {
                    let caps = reg3_1.captures(line).expect("cap error");
                    let alp = caps[1].parse::<f64>().expect("parse error");
                    let exp = caps[2].parse::<f64>().expect("parse error");
                    let ang = syms
                        .iter()
                        .position(|&x| x == sym)
                        .expect(&format!("get pos error: {}", sym));
                    datas.push((atm, shl - 1, ang, alp, exp));
                    continue;
                }
                if reg3_2.is_match(line) {
                    let caps = reg3_2.captures(line).expect("cap error");
                    let alp = caps[1].parse::<f64>().expect("parse error");
                    let exp = caps[2].parse::<f64>().expect("parse error");
                    let ang = syms.iter().position(|&x| x == sym).expect("get pos error");
                    datas.push((atm, shl - 1, ang, alp, exp));
                    continue;
                }

                if line.starts_with("[") {
                    log::debug!("basis 读取结束: {}|{}", l + 1, line);
                    break;
                }
            }
            if let Some(last) = ato_atms.last() {
                if *last == 0 {
                    log::warn!("最后一个原子轨道的索引为0，可能存在问题");
                }
            } else {
                panic!("未读取到任何原子轨道数据");
            };
            let ato_syms = ato_syms.iter().map(|sym| sym.to_string()).collect();
            (ato_atms, ato_shls, ato_syms, datas)
        });
        (ato_atms, ato_shls, ato_syms, datas)
    }

    /// 读取轨道系数
    pub fn read_coefs(&self) -> (Vec<f64>, Vec<f64>, Vec<Vec<f64>>) {
        let reg2 = Regex::new(r" Ene= +(-?\d+.\d+)").expect("reg new error");
        let reg4 = Regex::new(r" Occup= +(\d.\d+)").expect("reg new error");
        let reg5 = Regex::new(r"^ *(\d+) +(-?\d+.\d+)$").expect("reg new error");
        let mut obt_occs = Vec::new();
        let mut obt_engs = Vec::new();
        let mut cmat = Vec::new();
        let l_mo = self.marks.get("MO").expect("get MO line error");
        let lines = &self.lines[(l_mo + 1) as usize..];

        for (l, line) in lines.iter().enumerate() {
            if line == "" {
                continue;
            }
            if line == " " {
                continue;
            }
            if line.starts_with(" Sym=") {
                continue;
            }
            if line.starts_with(" Spin=") {
                continue;
            }
            if reg2.is_match(line) {
                let caps = reg2.captures(line).expect("cap error");
                let eng = caps[1].parse::<f64>().expect("parse error");
                obt_engs.push(eng);
                cmat.push(vec![]);
                continue;
            }
            if reg4.is_match(line) {
                let caps = reg4.captures(line).expect("cap error");
                let occ = caps[1].parse::<f64>().expect("parse error");
                obt_occs.push(occ);
                continue;
            }
            if reg5.is_match(line) {
                let caps = reg5.captures(line).expect("cap error");
                let coe = caps[2].parse::<f64>().expect("parse error");
                let row = cmat.len() - 1;
                cmat[row].push(coe);
                continue;
            }
            log::error!("未匹配类型错误");
            panic!("未匹配类型错误: {}:{}", l + l_mo + 1, line);
        }
        (obt_engs, obt_occs, cmat)
    }

    /// 读取电子数
    pub fn read_neles(&self) -> [usize; 2] {
        let occs = self.read_coefs().1;
        let nobt = occs.len();
        let mut neles = [0; 2];
        for occ in occs.iter() {
            if *occ == 0.0 {
                continue;
            }
            neles[0] += 1;
            neles[1] += 1;
        }
        neles
    }
}

impl Reader for MdeReader {
    fn get_atoms(&self) -> Atoms {
        let (atms, xyzs) = self.read_geome();
        let crgs: Vec<f64> = atms.iter().map(|atm| *atm as f64).collect();
        let mut atoms = Atoms::new(&atms, &xyzs, &crgs);
        atoms.neles = self.read_neles();
        atoms
    }

    fn get_basis(&self) -> Option<Basis> {
        let (ato_atms, ato_shls, ato_syms, datas) = self.read_basis();
        let mut datas = datas.clone();
        for data in datas.iter_mut() {
            let ang = data.2;
            let alp = data.3;
            let coe = data.4;
            let nm = (2.0 * alp / PI).powf(0.75) * (4.0 * alp).powf(ang as f64 / 2.0);
            data.4 = coe / nm;
        }

        datas.sort_by_key(|data| (data.0, data.1, data.2));

        let forms = shl_forms(ato_atms, ato_shls, ato_syms);
        Some(shl_basis(&datas, &forms))
    }

    fn get_coefs(&self) -> Option<Coefs> {
        let (ato_atms, ato_shls, ato_syms, _) = self.read_basis();
        let (obt_engs, obt_occs, cmat) = self.read_coefs();
        let nobt = cmat.len();
        let nato = cmat[0].len();
        let cmat = Array2::from_shape_vec((nobt, nato), cmat.into_iter().flatten().collect()).expect("cmat error");
        let cmat = cmat.t().to_owned();
        let coefs = Coefs::new(
            ato_atms.clone(),
            ato_shls.clone(),
            ato_syms.clone(),
            obt_engs,
            obt_occs,
            cmat,
        );
        Some(coefs)
    }
}
