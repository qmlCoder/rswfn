use ndarray::Array2;
use regex::Regex;

use crate::base::{Atoms, Basis, Coefs};
use crate::datas::elems::sym2atm;
use crate::reader::{Reader, read_lines, shl_basis, shl_forms};
use std::collections::HashMap;
use std::vec;

/// Orca输出文件(.out)读取器
pub struct OutReader {
    /// 文件中各数据块的位置标记
    pub marks: HashMap<String, usize>,
    /// 文件行
    pub lines: Vec<String>,
    /// 原子数量
    pub natm: usize,
    /// 原子轨道数量
    pub nato: usize,
}

impl OutReader {
    /// 创建OutReader实例,解析文件中的标记位置
    pub fn new(path: &str) -> OutReader {
        let mut marks = HashMap::new();
        // println!("读取out:{}", path);
        let lines = read_lines(path);
        let mut nato = 0;
        let mut natm = 0;
        for (l, line) in lines.iter().enumerate() {
            if line == "CARTESIAN COORDINATES (A.U.)" {
                marks.entry("geome".to_string()).or_insert(l);
            }
            if line == "BASIS SET IN INPUT FORMAT" {
                marks.entry("basis".to_string()).or_insert(l);
            }
            if line == "ORBITAL ENERGIES" {
                marks.entry("obtEngs".to_string()).or_insert(l);
            }
            if line == "MOLECULAR ORBITALS" {
                marks.entry("coefs".to_string()).or_insert(l);
            }
            if line.starts_with("Number of atoms                             ...") {
                natm = line[47..].trim().parse().unwrap();
            }
            if line.starts_with("Number of basis functions                   ...") {
                nato = line[47..].trim().parse().unwrap();
            }
        }
        OutReader { marks, lines, natm, nato }
    }

    /// 读取几何结构
    pub fn read_geome(&self) -> (Vec<String>, Vec<[f64; 3]>) {
        let mark = self.marks.get("geome").unwrap();
        let mut syms = vec![];
        let mut xyzs = vec![];
        for line in self.lines.iter().skip(*mark + 3) {
            if line == "" {
                break;
            }
            log::debug!("line: {line}");
            let vals: Vec<&str> = line.split_whitespace().into_iter().collect();
            let sym = vals[1].to_string();
            let x: f64 = vals[5].parse().unwrap();
            let y: f64 = vals[6].parse().unwrap();
            let z: f64 = vals[7].parse().unwrap();
            syms.push(sym);
            xyzs.push([x, y, z]);
        }
        (syms, xyzs)
    }

    /// 读取基组信息
    pub fn read_basis(&self) -> HashMap<usize, Vec<(i32, usize, f64, f64)>> {
        let mark = self.marks.get("basis").unwrap();
        let reg1 = Regex::new(r"^ ([SPD]) \d ").expect("正则表达式构建失败");
        let reg2 = Regex::new(r" +\d+ +(\d+.\d+) +(\d+.\d+)").expect("正则表达式构建失败");
        let mut data: HashMap<usize, Vec<(i32, usize, f64, f64)>> = HashMap::new();
        let mut atm = 0;
        let mut ang = 0;
        let mut shl = -1;
        for line in self.lines.iter().skip(*mark + 3) {
            if line.starts_with(" # Basis set for element :") {
                let elm = line[26..].trim().to_string();
                atm = sym2atm(&elm);
            }
            if reg1.is_match(line) {
                let caps = reg1.captures(line).unwrap();
                let sym = caps[1].trim().to_string();
                ang = sym2ang(&sym);
                shl += 1;
            }
            if reg2.is_match(line) {
                let caps = reg2.captures(line).unwrap();
                let alp: f64 = caps[1].parse().unwrap();
                let coe: f64 = caps[2].parse().unwrap();

                data.entry(atm).or_insert_with(Vec::new).push((shl, ang, alp, coe));
            }
            if line == "  end;" {
                shl = -1;
            }
            if line == "---------------------------------" {
                break;
            }
        }
        data
    }

    /// 读取轨道系数
    pub fn read_coefs(&self) -> (Vec<usize>, Vec<usize>, Vec<String>, Vec<f64>, Vec<f64>, Array2<f64>) {
        let mark = self.marks.get("coefs").unwrap();
        let nato = self.nato;
        let block_size = nato + 4;
        let block_numb = nato / 6 + if nato % 6 == 0 { 0 } else { 1 };

        let mut obt_engs = vec![];
        let mut obt_occs = vec![];
        let mut ato_atms = vec![];
        let mut ato_shls = vec![];
        let mut ato_syms = vec![];
        let mut cmat = Array2::zeros((nato, nato));
        let mut lidx = mark + 2;
        for ib in 0..block_numb {
            log::debug!("lidx: {}, {}", lidx, self.lines[lidx]);
            let idx_part = &self.lines[lidx + 0][16..];
            log::debug!("idx_part: {idx_part}");
            let idxs: Vec<usize> = idx_part.char_indices().step_by(10).map(|(i, _)| idx_part[i..i + 10].trim().parse().unwrap()).collect();
            let eng_part = &self.lines[lidx + 1][16..];
            let engs: Vec<f64> = eng_part.char_indices().step_by(10).map(|(i, _)| eng_part[i..i + 10].trim().parse().unwrap()).collect();
            obt_engs.extend(engs);
            let occ_part = &self.lines[lidx + 2][16..];
            let occs: Vec<f64> = occ_part.char_indices().step_by(10).map(|(i, _)| occ_part[i..i + 10].trim().parse().unwrap()).collect();
            obt_occs.extend(occs);
            for ato in 0..nato {
                let line = &self.lines[lidx + 4 + ato];
                let coef_part = &line[16..];
                let coefs: Vec<f64> = coef_part.char_indices().step_by(10).map(|(i, _)| coef_part[i..i + 10].trim().parse().unwrap()).collect();
                for (i, obt) in idxs.iter().enumerate() {
                    cmat[[ato, *obt]] = coefs[i];
                }
                if ib == 0 {
                    let atm: usize = line[0..3].trim().parse().unwrap();
                    ato_atms.push(atm);
                    let shl: usize = line[5..8].trim().parse().unwrap();
                    ato_shls.push(shl);
                    let sym = line[8..16].trim().to_string();
                    ato_syms.push(sym);
                }
            }
            lidx += block_size;
        }
        (ato_atms, ato_shls, ato_syms, obt_engs, obt_occs, cmat)
    }

    /// 读取轨道能量和占据数
    pub fn read_obtEngs(&self) -> (Vec<f64>, Vec<f64>) {
        let mark = self.marks.get("obtEngs").unwrap();
        let mut occs = vec![];
        let mut engs = vec![];
        for line in self.lines.iter().skip(*mark + 4) {
            let occ: f64 = line[7..13].parse().unwrap();
            let eng: f64 = line[13..28].parse().unwrap();
            occs.push(occ);
            engs.push(eng);
        }
        (occs, engs)
    }
}

impl Reader for OutReader {
    fn get_atoms(&self) -> Atoms {
        let (syms, xyzs) = self.read_geome();
        let atms: Vec<usize> = syms.iter().map(|sym| sym2atm(sym)).collect();
        let crgs: Vec<f64> = atms.iter().map(|atm| *atm as f64).collect();
        Atoms::new(&atms, &xyzs, &crgs)
    }

    fn get_basis(&self) -> Option<Basis> {
        let mut basis = Basis::new(vec![], None);
        let data = self.read_basis();
        let syms = self.read_geome().0;
        let atms: Vec<usize> = syms.iter().map(|sym| sym2atm(sym)).collect();
        let mut datas = vec![];
        for (idx, atm) in atms.iter().enumerate() {
            for (shl, ang, alp, coe) in data.get(atm).unwrap().iter() {
                let data = (idx, *shl as usize, *ang, *alp, *coe);
                datas.push(data);
            }
        }
        let (ato_atms, ato_shls, ato_syms, _, _, cmat) = self.read_coefs();
        let forms = shl_forms(&ato_atms, &ato_shls, &ato_syms);
        let basis = shl_basis(&datas, &forms);
        Some(basis)
    }

    fn get_coefs(&self) -> Option<Coefs> {
        let (ato_atms, ato_shls, ato_syms, obt_engs, obt_occs, cmat) = self.read_coefs();
        let obt_occs = obt_occs.iter().map(|occ| *occ).collect();
        let ato_syms = ato_syms.iter().map(|sym| sym2sym(sym)).collect();
        let coefs = Coefs::new(ato_atms, ato_shls, ato_syms, obt_engs, obt_occs, cmat);
        Some(coefs)
    }
}

fn sym2ang(sym: &str) -> usize {
    match sym {
        "S" => 0,
        "P" => 1,
        "D" => 2,
        "F" => 3,
        _ => panic!("error sym"),
    }
}

// 将orca的原子轨道符号转为gaussian的
fn sym2sym(sym: &str) -> String {
    let sym = match sym {
        "pz" => "PZ",
        "px" => "PX",
        "py" => "PY",
        "dz2" => "D 0",
        "dxz" => "D+1",
        "dyz" => "D-1",
        "dx2y2" => "D+2",
        "dxy" => "D-2",
        _ => panic!("error sym"),
    };
    sym.to_string()
}
