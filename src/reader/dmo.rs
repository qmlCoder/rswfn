use crate::base::{Atoms, Bas, Basis, Coefs, Shl, Stm, Vec3};
use crate::datas::consts::BOHR;
use crate::datas::elems::{self, sym2atm};
use crate::maths;
use crate::reader::Reader;
use crate::reader::read_lines;
use crate::utils;
use ndarray::{Array2, s};
use regex::Regex;
use std::collections::HashMap;
use std::sync::OnceLock;
/// dmi/dmo文件读取器
///
/// 读取自定义的dmi/dmo格式文件,包含原子轨道信息。
pub struct DmoReader {
    pub path: String,
    pub lines: Vec<String>,
    pub marks: HashMap<String, usize>,
    atms: OnceLock<(Vec<usize>, Vec<[f64; 3]>)>,
    dirs: OnceLock<Vec<(usize, usize, String, [f64; 3])>>,
}

impl DmoReader {
    /// 创建DmoReader实例,解析文件中的标记位置
    pub fn new(path: &str) -> Self {
        let mut marks = HashMap::new();
        let lines = read_lines(path);
        for (l, line) in lines.iter().enumerate() {
            if line.starts_with("#atms") {
                marks.entry("atms".to_string()).or_insert(l);
                continue;
            }
            if line.starts_with("#atos") {
                marks.entry("atos".to_string()).or_insert(l);
                continue;
            }
            if line.starts_with("#BM") {
                marks.entry("BM".to_string()).or_insert(l);
                continue;
            }
            if line.starts_with("#CM") {
                marks.entry("CM".to_string()).or_insert(l);
                continue;
            }
        }
        let path = path.to_string();
        let atms = OnceLock::new();
        let dirs = OnceLock::new();
        Self {
            path,
            lines,
            marks,
            atms,
            dirs,
        }
    }

    /// 读取原子信息(带缓存)
    pub fn read_atms(&self) -> &(Vec<usize>, Vec<[f64; 3]>) {
        self.atms.get_or_init(|| {
            let mut atms = vec![];
            let mut xyzs = vec![];
            let s1 = Regex::new(r" ([A-Za-z\d]+) +(-?\d+.\d+) +(-?\d+.\d+) +(-?\d+.\d+)").expect("正则表达式生成失败");
            for line in self.lines.iter().skip(self.marks["atms"] + 1) {
                if s1.is_match(line) {
                    let caps = s1.captures(line).expect("捕获失败");
                    let sym = caps[1].to_string();
                    let atm = sym2atm(&sym);
                    let x = caps[2].parse::<f64>().expect("解析浮点数失败") * BOHR;
                    let y = caps[3].parse::<f64>().expect("解析浮点数失败") * BOHR;
                    let z = caps[4].parse::<f64>().expect("解析浮点数失败") * BOHR;
                    atms.push(atm);
                    xyzs.push([x, y, z]);
                } else {
                    break;
                }
            }
            (atms, xyzs)
        })
    }

    /// 读取原子轨道信息(带缓存)
    pub fn read_atos(&self) -> &Vec<(usize, usize, String, [f64; 3])> {
        self.dirs.get_or_init(|| {
            let mut atos = vec![];
            if self.marks.contains_key("atos") {
                let s1 = Regex::new(r" *(\d+) *(\d+) *([A-Za-z\d]+) +(-?\d+.\d+) +(-?\d+.\d+) +(-?\d+.\d+)")
                    .expect("正则表达式生成失败");
                for line in self.lines.iter().skip(self.marks["atos"] + 1) {
                    if s1.is_match(line) {
                        let caps = s1.captures(line).expect("捕获失败");
                        let atm = caps[1].parse::<usize>().expect("转换失败");
                        let shl = caps[2].parse::<usize>().expect("转换失败");
                        let sym = caps[3].to_string();
                        let dx = caps[4].parse::<f64>().expect("解析浮点数失败") * BOHR;
                        let dy = caps[5].parse::<f64>().expect("解析浮点数失败") * BOHR;
                        let dz = caps[6].parse::<f64>().expect("解析浮点数失败") * BOHR;
                        let l = maths::vec_len(&[dx, dy, dz]);
                        atos.push((atm - 1, shl, sym, [dx / l, dy / l, dz / l]));
                    } else {
                        break;
                    }
                }
            }
            atos
        })
    }

    /// 读取键连矩阵
    pub fn read_BM(&self) -> Array2<f64> {
        let atos = &self.read_atos();
        let nato = atos.len();
        let mut BM = Array2::zeros((nato, nato));
        let s1 = Regex::new(r" *(\d+) +(\d+) +(\d+) +(\d+) +(-?\d+.\d+)").unwrap();
        for line in self.lines.iter().skip(self.marks["BM"] + 1) {
            if s1.is_match(line) {
                let caps = s1.captures(line).expect("捕获失败");
                let a1 = caps[1].parse::<usize>().expect("转换失败") - 1;
                let a2 = caps[2].parse::<usize>().expect("转换失败") - 1;
                let i1 = caps[3].parse::<usize>().expect("转换失败") - 1;
                let i2 = caps[4].parse::<usize>().expect("转换失败") - 1;
                let val = caps[5].parse::<f64>().expect("转换失败");
                BM[(i1, i2)] = val;
                BM[(i2, i1)] = val;
            } else {
                continue;
            }
        }
        BM
    }

    /// 读取轨道系数和能量
    pub fn read_coefs(&self) -> (Vec<f64>, Array2<f64>) {
        let atos = &self.read_atos();
        let nato = atos.len();
        let start = self.marks["CM"] + 1;
        let mut CM = Array2::zeros((nato, nato));
        let mut obt_engs = vec![];
        for (ib, lb) in (start..self.lines.len()).step_by(nato + 2).enumerate() {
            let engs: Vec<f64> = self.lines[lb + 1]
                .trim()
                .split_whitespace()
                .map(|eng| eng.parse::<f64>().unwrap())
                .collect();
            obt_engs.extend(engs);
            for ato in 0..nato {
                let line = &self.lines[lb + 2 + ato];
                // let atm = line[0..10].parse::<usize>().expect("转换失败");
                // let sym = line[10..20].parse::<usize>().expect("转换失败");
                let vals: Vec<f64> = line[10..]
                    .trim()
                    .split_whitespace()
                    .map(|val| val.parse::<f64>().unwrap())
                    .collect();
                for (j, val) in vals.iter().enumerate() {
                    CM[(ato, ib * 10 + j)] = *val;
                }
            }
        }
        (obt_engs, CM)
    }
}

impl Reader for DmoReader {
    fn get_atoms(&self) -> Atoms {
        let (atms, xyzs) = self.read_atms();
        let atos = self.read_atos();
        let crgs: Vec<f64> = atms.iter().map(|atm| *atm as f64).collect();
        let natm = atms.len();
        let mut stms = vec![Stm::blank(); natm];
        let nato = atos.len();
        for i in 0..nato {
            let [zx, xy, zz] = atos[i].3;
            let ez = Vec3::new(zx, xy, zz);
            let ex = Vec3::rand();
            let ey = ez.cross(&ex).normalize();
            let ex = ez.cross(&ey).normalize();
            let stm = Stm::new(ex.into_arr(), ey.into_arr(), ez.into_arr());
            stms[atos[i].0] = stm;
        }
        Atoms::new(&atms, &xyzs, &crgs)
    }
    fn get_basis(&self) -> Option<Basis> {
        if self.path.ends_with(".dmi") {
            return None;
        }
        let (atms, _) = self.read_atms();
        let atos = self.read_atos();
        let nato = atos.len();
        let mut shls = vec![];
        let mut keys = vec![];
        for i in 0..nato {
            let atm = atos[i].0;
            let shl = atos[i].1;
            let sym = &atos[i].2;
            let key = (atm, shl, sym);
            if !keys.contains(&key) {
                keys.push(key);
                let nlm = sym2nlm(atms[atm], shl, sym);
                let ang = 1;
                let shl = Shl::new(atm, shl, ang, "car".to_string(), vec![], vec![]);
                shls.push(shl);
            }
        }
        let coefs = self.get_coefs();
        Some(Basis::new(shls, None))
    }

    fn get_coefs(&self) -> Option<Coefs> {
        if self.path.ends_with(".dmi") {
            return None;
        }
        let atos = self.read_atos();
        let nato = atos.len();
        if nato == 0 {
            return None;
        }
        let (obt_engs, CM) = self.read_coefs();
        let mut obt_occs = vec![];
        obt_occs.extend(vec![2.0; nato / 2]);
        obt_occs.extend(vec![0.0; nato / 2]);

        // let ato_syms = ato_syms.clone();
        let mut atms = vec![];
        let mut shls = vec![];
        let mut syms = vec![];
        let mut cmat = Array2::zeros((nato * 3, nato));
        for i in 0..nato {
            let atm = atos[i].0;
            let shl = atos[i].1;
            // let sym = ato_syms[i].clone();
            atms.extend(vec![atm; 3]);
            shls.extend(vec![shl; 3]);
            syms.extend(vec!["PX".to_owned(), "PY".to_owned(), "PZ".to_owned()]);
            let info = s![i * 3 + 2, ..];
            cmat.slice_mut(info).assign(&CM.slice(s![i, ..]));
        }
        let coefs = Coefs::new(atms, shls, syms, obt_engs, obt_occs, cmat);
        Some(coefs)
    }
}

// 原子轨道符号转为量子数
fn sym2nlm(nuc: usize, shl: usize, sym: &str) -> (usize, usize, isize) {
    let n = elems::atm2row(nuc);
    match (shl, sym) {
        (1, "S") => (n, 0, 0),
        (2, "PX") => (n, 1, 1),
        (2, "PY") => (n, 1, -1),
        (2, "PZ") => (n, 1, 0),
        _ => panic!("无法识别的符号"),
    }
}
