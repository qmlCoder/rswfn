#![allow(non_snake_case)]

use ndarray::{Array2, Axis, concatenate};
use std::collections::HashMap;
use std::fmt::Debug;
use std::sync::LazyLock;

use crate::base::{Atoms, Bas, Basis, Coefs, Shl, Stm};
use crate::datas::btrans;
use crate::reader::{Reader, read_lines, shl_basis, shl_forms};

static SHL2SYM: LazyLock<HashMap<isize, Vec<&str>>> = LazyLock::new(|| {
    let mut map = HashMap::new();
    map.insert(-4, btrans::SPH_G_SYMS.to_vec());
    map.insert(-3, btrans::SPH_F_SYMS.to_vec());
    map.insert(-2, btrans::SPH_D_SYMS.to_vec());
    map.insert(-1, vec!["S", "PX", "PY", "PZ"]);
    map.insert(0, vec!["S"]);
    map.insert(1, vec!["PX", "PY", "PZ"]);
    map.insert(2, btrans::CAR_D_SYMS.to_vec());
    map.insert(3, btrans::CAR_F_SYMS.to_vec());
    map.insert(4, btrans::CAR_G_SYMS.to_vec());
    map
});

/// 检查点文件(.fch/.fchk)中存储的值
pub struct Vals {
    I: Vec<isize>,
    R: Vec<f64>,
}

#[derive(Debug, Clone)]
/// fch文件的原子轨道计数信息
pub struct Counts {
    /// 原子轨道数
    pub nato: usize,
    /// 分子轨道数
    pub nobt: usize,
    /// α电子数
    pub nela: usize,
    /// β电子数
    pub nelb: usize,
}

impl Counts {
    /// 创建空的Counts实例
    pub fn blank() -> Self {
        Counts { nato: 0, nobt: 0, nela: 0, nelb: 0 }
    }
}

/// 高斯检查点文件(.fch/.fchk)读取器
pub struct FchReader {
    /// 文件路径
    pub path: String,
    /// 各数据块的标题行行号
    pub marks: HashMap<String, usize>,
    /// 原子轨道与电子计数信息
    pub count: Counts,
    /// 文件行
    pub lines: Vec<String>,
}

impl FchReader {
    /// 创建FchReader实例,解析文件中的标记行
    pub fn new(path: &str) -> Self {
        let mut count = Counts::blank();
        let mut marks = HashMap::new();
        let lines = read_lines(path);
        let keys = vec![
            "Atomic numbers",
            "Current cartesian coordinates",
            "Primitive exponents",
            "Contraction coefficients",
            "P(S=P) Contraction coefficients",
            "Alpha MO coefficients",
            "Beta MO coefficients",
            "Alpha Orbital Energies",
            "Beta Orbital Energies",
            "Shell types",
            "Shell to atom map",
            "Number of primitives per shell",
            "Nuclear charges",
            "Local Coord Systems",
        ];
        for (l, line) in lines.iter().enumerate() {
            if line.len() < 40 {
                continue;
            }
            let mark = line[..40].trim().to_string();
            if mark == "Number of alpha electrons" {
                count.nela = line[49..61].trim().to_string().parse().expect("转换失败");
                continue;
            }
            if mark == "Number of beta electrons" {
                count.nelb = line[49..61].trim().to_string().parse().expect("转换失败");
                continue;
            }
            if mark == "Number of basis functions" {
                count.nato = line[49..61].trim().to_string().parse().expect("转换失败");
                continue;
            }
            if mark == "Number of independent functions" {
                count.nobt = line[49..61].trim().to_string().parse().expect("转换失败");
                continue;
            }
            for key in keys.iter() {
                let key = (*key).to_string();
                if mark == key {
                    marks.entry(key).or_insert(l as usize);
                }
            }
        }

        FchReader {
            path: path.to_string(),
            marks: marks,
            count: count,
            lines: lines,
        }
    }

    /// 解析fch文件中的数据块
    ///
    /// 根据标记名称读取对应的数据块,返回整数或浮点数向量。
    pub fn parse_block(&self, mark: &str) -> Vals {
        let mark = mark.to_string();
        if !self.marks.contains_key(&mark) {
            return Vals { I: vec![], R: vec![] };
        }
        let mark_num = self.marks.get(&mark).unwrap_or_else(|| panic!("获取行失败{}", mark));
        let mark_txt = self.lines.get(*mark_num as usize).expect("获取标题行失败");
        let vtype = &mark_txt[43..44];
        let count = mark_txt[49..61].trim().parse::<usize>().expect("转换失败");
        let mut r_vals = Vec::new();
        let mut i_vals = Vec::new();
        match vtype {
            "I" => {
                let num_line = count / 6 + if count % 6 == 0 { 0 } else { 1 };
                for l in mark_num + 1..mark_num + num_line + 1 {
                    let line = self.lines.get(l as usize).expect("获取行失败");
                    let vals = line.split_whitespace().map(|s| s.parse().expect("转换失败")).collect::<Vec<isize>>();
                    i_vals.extend(vals);
                }
            }
            "R" => {
                let num_line = count / 5 + if count % 5 == 0 { 0 } else { 1 };
                for l in mark_num + 1..mark_num + num_line + 1 {
                    let line = self.lines.get(l as usize).expect("获取行失败");
                    let vals = line.split_whitespace().map(|s| s.parse().expect("转换失败")).collect::<Vec<f64>>();
                    r_vals.extend(vals);
                }
            }
            _ => panic!("未知类型:{}", vtype),
        }
        Vals { I: i_vals, R: r_vals }
    }

    /// 读取原子序数
    pub fn read_atms(&self) -> Vec<isize> {
        self.parse_block("Atomic numbers").I
    }

    /// 读取原子核电荷(无ECP时为原子序数)
    pub fn read_crgs(&self) -> Vec<f64> {
        self.parse_block("Nuclear charges").R
    }

    /// 读取笛卡尔坐标
    pub fn read_xyzs(&self) -> Vec<[f64; 3]> {
        let vals = self.parse_block("Current cartesian coordinates").R;
        let mut xyzs = Vec::new();
        for i in (0..vals.len()).step_by(3) {
            let x = vals[i + 0];
            let y = vals[i + 1];
            let z = vals[i + 2];
            xyzs.push([x, y, z]);
        }
        xyzs
    }

    /// 读取局部坐标系
    pub fn read_stms(&self) -> Option<Vec<Stm>> {
        let vals = self.parse_block("Local Coord Systems").R;
        if vals.len() == 0 {
            return None;
        }
        let mut stms = vec![];

        for i in (0..vals.len()).step_by(9) {
            let xx = vals[i + 0];
            let xy = vals[i + 1];
            let xz = vals[i + 2];
            let yx = vals[i + 3];
            let yy = vals[i + 4];
            let yz = vals[i + 5];
            let zx = vals[i + 6];
            let zy = vals[i + 7];
            let zz = vals[i + 8];
            let ex = [xx, xy, xz];
            let ey = [yx, yy, yz];
            let ez = [zx, zy, zz];
            let stm = Stm::new(ex, ey, ez);
            stms.push(stm);
        }
        Some(stms)
    }

    /// 读取轨道系数及其他轨道相关信息
    pub fn read_coefs(&self) -> (Vec<usize>, Vec<usize>, Vec<String>, Vec<f64>, Vec<f64>, Array2<f64>) {
        let nato = self.count.nato;
        let nobt = self.count.nobt;
        let coefs_a = self.parse_block("Alpha MO coefficients").R;
        assert_eq!(coefs_a.len(), nato * nobt);
        let cmat_a = Array2::from_shape_vec((nobt as usize, nato as usize), coefs_a).expect("数组转矩阵失败").t().to_owned();
        let open = self.marks.contains_key("Beta MO coefficients");
        let cmat = if open {
            let coefs_b = self.parse_block("Beta MO coefficients").R;

            let cmat_b = Array2::from_shape_vec((nobt as usize, nato as usize), coefs_b).expect("数组转矩阵失败").t().to_owned();
            let cmat = concatenate(Axis(1), &[cmat_a.view(), cmat_b.view()]).expect("矩阵拼接失败");
            cmat
        } else {
            cmat_a
        };
        let shl_type = self.parse_block("Shell types").I;
        let shl_atms = self.parse_block("Shell to atom map").I;
        let shl_atms = shl_atms.iter().map(|x| (*x - 1) as usize).collect::<Vec<usize>>();
        let mut ato_atms = Vec::new();
        let mut ato_shls = Vec::new();
        let mut ato_syms = Vec::new();

        let mut atm = -1;
        let mut shl = 0;

        assert_ne!(shl_type.len(), 0);
        for (i, st) in shl_type.iter().enumerate() {
            let syms = SHL2SYM.get(st).expect(&format!("获取符号失败:{}", st));
            if atm != shl_atms[i] as isize {
                shl = 0;
                atm = shl_atms[i] as isize;
            } else {
                shl += 1;
            }
            for s in 0..syms.len() {
                ato_atms.push(shl_atms[i]);
                ato_syms.push(syms[s].to_string());
                ato_shls.push(shl);
            }
        }
        let mut engs = self.parse_block("Alpha Orbital Energies").R;
        if open {
            let engs_b = self.parse_block("Beta Orbital Energies").R;
            engs.extend(engs_b);
        }
        let mut occs = Vec::new();
        let nela = self.count.nela;
        let nelb = self.count.nelb;
        if open {
            for i in 0..self.count.nobt {
                occs.push(if i < nela { 1.0 } else { 0.0 });
            }
            for i in 0..self.count.nobt {
                occs.push(if i < nelb { 1.0 } else { 0.0 });
            }
        } else {
            for i in 0..self.count.nobt {
                occs.push(if i < nela { 2.0 } else { 0.0 });
            }
        }
        assert_ne!(engs.len(), 0);
        assert_ne!(occs.len(), 0);
        assert_ne!(ato_atms.len(), 0);
        assert_ne!(ato_shls.len(), 0);
        assert_ne!(ato_syms.len(), 0);
        (ato_atms, ato_shls, ato_syms, engs, occs, cmat)
    }
}

impl Reader for FchReader {
    fn get_atoms(&self) -> Atoms {
        let atms = self.read_atms().iter().map(|x| *x as usize).collect::<Vec<usize>>();
        let xyzs = self.read_xyzs();
        let crgs = self.read_crgs();
        let mut atoms = Atoms::new(&atms, &xyzs, &crgs);
        atoms.neles = [self.count.nela, self.count.nelb];
        atoms
    }

    fn get_basis(&self) -> Option<Basis> {
        let shl_type = self.parse_block("Shell types").I;
        let shl_nums = self.parse_block("Number of primitives per shell").I;
        let shl_atms = self.parse_block("Shell to atom map").I;
        let shl_atms = shl_atms.iter().map(|x| *x as usize - 1).collect::<Vec<usize>>();
        let pm_alps = self.parse_block("Primitive exponents").R;
        let pm_coes = self.parse_block("Contraction coefficients").R;
        let mut sp_coes = Vec::new();
        if self.marks.contains_key("P(S=P) Contraction coefficients") {
            sp_coes = self.parse_block("P(S=P) Contraction coefficients").R;
        }
        let mut atm = 0usize;
        let mut shl = 0usize;
        let mut idx = 0;
        let mut datas = Vec::new();
        for i in 0..shl_type.len() {
            let shl_atm = shl_atms[i] as usize;
            let shl_num = shl_nums[i] as usize;
            if atm != shl_atm {
                atm = shl_atm;
                shl = 0;
            }
            for j in 0..shl_num as usize {
                if shl_type[i] == -1 {
                    datas.push((shl_atm, shl, 0, pm_alps[idx], pm_coes[idx]));
                    datas.push((shl_atm, shl, 1, pm_alps[idx], sp_coes[idx]));
                } else {
                    let ang = shl_type[i].abs() as usize;
                    datas.push((shl_atm, shl, ang, pm_alps[idx], pm_coes[idx]));
                }
                idx += 1;
            }
            shl += 1;
        }
        datas.sort_by_key(|data| (data.0, data.1, data.2));
        let (ato_atms, ato_shls, ato_syms, _, _, _) = self.read_coefs();
        let forms = shl_forms(&ato_atms, &ato_shls, &ato_syms);
        Some(shl_basis(&datas, &forms))
    }

    fn get_coefs(&self) -> Option<Coefs> {
        let (ato_atms, ato_shls, ato_syms, engs, occs, cmat) = self.read_coefs();
        let coefs = Coefs::new(ato_atms, ato_shls, ato_syms, engs, occs, cmat);
        Some(coefs)
    }
}
