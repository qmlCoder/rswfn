use ndarray::Array2;
use std::collections::HashMap;

use crate::base::{Atoms, Basis, Coefs};
use crate::datas::btrans;
use crate::datas::consts::BOHR;
use crate::reader::{BasData, Reader, read_lines, shl_basis};

/// Multiwfn 波函数文件(.mwfn)读取器
///
/// 按照 "mwfn: A Strict, Concise and Extensible Format for Electronic
/// Wavefunction Storage and Exchange" 中定义的格式, 依次读取系统信息、
/// 原子信息、基函数信息、轨道信息, 以及可选的矩阵信息之外的其余字段。
pub struct MfnReader {
    path: String,
    lines: Vec<String>,
    marks: HashMap<String, usize>,
}

impl MfnReader {
    /// 创建MfnReader实例,解析文件中各字段的位置
    pub fn new(path: &str) -> Self {
        let lines = read_lines(path);
        let mut marks = HashMap::new();
        let keys = [
            "Wfntype",
            "Charge",
            "Naelec",
            "Nbelec",
            "E_tot",
            "VT_ratio",
            "Ndim",
            "Ncenter",
            "Nbasis",
            "Nindbasis",
            "Nprims",
            "Nshell",
            "Nprimshell",
        ];
        for (l, line) in lines.iter().enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let block = if line.starts_with("$Centers") {
                Some("Centers")
            } else if line.starts_with("$Shell types") {
                Some("Shell types")
            } else if line.starts_with("$Shell centers") {
                Some("Shell centers")
            } else if line.starts_with("$Shell contraction degrees") {
                Some("Shell degrees")
            } else if line.starts_with("$Primitive exponents") {
                Some("Primitive exponents")
            } else if line.starts_with("$Contraction coefficients") {
                Some("Contraction coefficients")
            } else {
                None
            };
            if let Some(key) = block {
                marks.entry(key.to_string()).or_insert(l);
                continue;
            }
            if let Some((key, _)) = line.split_once('=') {
                let key = key.trim();
                if keys.contains(&key) {
                    marks.entry(key.to_string()).or_insert(l);
                }
            }
        }
        Self {
            path: path.to_string(),
            lines,
            marks,
        }
    }

    /// 读取标量字段
    fn scalar(&self, key: &str) -> Option<f64> {
        let l = self.marks.get(key)?;
        let line = self.lines.get(*l)?;
        let (_, val) = line.split_once('=')?;
        Some(val.trim().parse::<f64>().expect(&format!("解析{}失败", key)))
    }

    /// 读取整数标量字段
    fn scalar_usize(&self, key: &str) -> usize {
        self.scalar(key).unwrap_or_else(|| panic!("缺少{}字段", key)) as usize
    }

    /// 读取向量/矩阵/列表数据块
    ///
    /// 数据块从标记行之后开始, 直到遇到空行、注释、下一个标记或标量字段为止。
    fn vals(&self, mark: &str) -> Vec<f64> {
        let l = self.marks.get(mark).unwrap_or_else(|| panic!("缺少${}数据块", mark));
        let mut vals = Vec::new();
        for line in self.lines.iter().skip(*l + 1) {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') || line.starts_with('$') || line.contains('=') {
                break;
            }
            vals.extend(line.split_whitespace().map(|s| s.parse::<f64>().expect("解析数值失败")));
        }
        vals
    }

    /// 读取几何结构
    ///
    /// 返回原子序数、笛卡尔坐标(内部使用Bohr)和有效核电荷数。
    pub fn read_geome(&self) -> (Vec<usize>, Vec<[f64; 3]>, Vec<f64>) {
        let l = self.marks.get("Centers").expect("获取$Centers失败");
        let mut atms = Vec::new();
        let mut xyzs = Vec::new();
        let mut crgs = Vec::new();
        for line in self.lines.iter().skip(*l + 1) {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') || line.starts_with('$') || line.contains('=') {
                break;
            }
            let parts: Vec<&str> = line.split_whitespace().collect();
            assert!(parts.len() >= 7, "原子信息行格式错误:{}", line);
            // 依次为: 序号 名称 元素序数 核电荷 x y z
            let atm = parts[2].parse::<usize>().expect("解析原子序数失败");
            let crg = parts[3].parse::<f64>().expect("解析核电荷失败");
            let x = parts[4].parse::<f64>().expect("解析坐标失败") * BOHR;
            let y = parts[5].parse::<f64>().expect("解析坐标失败") * BOHR;
            let z = parts[6].parse::<f64>().expect("解析坐标失败") * BOHR;
            atms.push(atm);
            crgs.push(crg);
            xyzs.push([x, y, z]);
        }
        (atms, xyzs, crgs)
    }

    /// 读取基函数信息
    ///
    /// 返回每个原子轨道对应的原子索引、壳层索引、符号, 原始基组数据以及
    /// 每个壳层的基函数形式("car"或"sph")。
    pub fn read_basis(&self) -> (Vec<usize>, Vec<usize>, Vec<String>, Vec<BasData>, Vec<String>) {
        let types: Vec<isize> = self.vals("Shell types").iter().map(|x| *x as isize).collect();
        let centers: Vec<usize> = self.vals("Shell centers").iter().map(|x| *x as usize - 1).collect();
        let degs: Vec<usize> = self.vals("Shell degrees").iter().map(|x| *x as usize).collect();
        let alps = self.vals("Primitive exponents");
        let coes = self.vals("Contraction coefficients");

        let nshl = types.len();
        assert_eq!(centers.len(), nshl, "$Shell centers数量与$Shell types不一致");
        assert_eq!(degs.len(), nshl, "$Shell contraction degrees数量与$Shell types不一致");
        assert_eq!(alps.len(), coes.len(), "原始指数与收缩系数数量不一致");
        assert_eq!(alps.len(), degs.iter().sum::<usize>(), "原始函数数量与收缩度不一致");

        let mut ato_atms = Vec::new();
        let mut ato_shls = Vec::new();
        let mut ato_syms = Vec::new();
        let mut datas: Vec<BasData> = Vec::new();
        let mut forms = Vec::with_capacity(nshl);
        let mut shl_counts: HashMap<usize, usize> = HashMap::new();
        let mut idx = 0;

        for i in 0..nshl {
            let atm = centers[i];
            let typ = types[i];
            // 正数表示笛卡尔壳层, 负数表示纯壳层(球谐)
            let ang = typ.unsigned_abs();
            let form = if typ < 0 { "sph" } else { "car" };
            let shl = *shl_counts.entry(atm).or_insert(0);
            *shl_counts.get_mut(&atm).unwrap() += 1;

            for sym in btrans::ang2syms(ang, form) {
                ato_atms.push(atm);
                ato_shls.push(shl);
                ato_syms.push(sym.to_string());
            }
            forms.push(form.to_string());

            for _ in 0..degs[i] {
                datas.push((atm, shl, ang, alps[idx], coes[idx]));
                idx += 1;
            }
        }

        (ato_atms, ato_shls, ato_syms, datas, forms)
    }

    /// 读取轨道系数、能量和占据数
    pub fn read_coefs(&self) -> (Vec<usize>, Vec<usize>, Vec<String>, Vec<f64>, Vec<f64>, Array2<f64>) {
        let (ato_atms, ato_shls, ato_syms, _, _) = self.read_basis();
        let nbasis = self.scalar_usize("Nbasis");
        assert_eq!(ato_atms.len(), nbasis, "基函数数量与Nbasis不一致");

        let start = self
            .lines
            .iter()
            .position(|line| line.trim_start().starts_with("Index="))
            .expect("未找到轨道信息");
        let mut obt_engs = Vec::new();
        let mut obt_occs = Vec::new();
        let mut cols: Vec<Vec<f64>> = Vec::new();
        let mut in_coeff = false;

        for line in self.lines.iter().skip(start) {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            if line.starts_with("Index=") {
                cols.push(Vec::with_capacity(nbasis));
                in_coeff = false;
                continue;
            }
            if line.starts_with("Type=") || line.starts_with("Sym=") {
                continue;
            }
            if line.starts_with("Energy=") {
                let eng = line
                    .split_once('=')
                    .unwrap()
                    .1
                    .trim()
                    .parse::<f64>()
                    .expect("解析轨道能量失败");
                obt_engs.push(eng);
                continue;
            }
            if line.starts_with("Occ=") {
                let occ = line
                    .split_once('=')
                    .unwrap()
                    .1
                    .trim()
                    .parse::<f64>()
                    .expect("解析轨道占据数失败");
                obt_occs.push(occ);
                continue;
            }
            if line.starts_with("$Coeff") {
                in_coeff = true;
                continue;
            }
            if line.starts_with('$') {
                // 轨道信息之后的其它矩阵字段, 与轨道系数无关
                in_coeff = false;
                continue;
            }
            if in_coeff && !cols.is_empty() {
                let vals = line
                    .split_whitespace()
                    .map(|s| s.parse::<f64>().expect("解析轨道系数失败"));
                cols.last_mut().unwrap().extend(vals);
            }
        }

        let nobt = cols.len();
        assert_eq!(obt_engs.len(), nobt, "轨道能量数量与轨道数量不一致");
        assert_eq!(obt_occs.len(), nobt, "轨道占据数数量与轨道数量不一致");

        let mut cmat = Array2::zeros((nbasis, nobt));
        for (j, col) in cols.iter().enumerate() {
            assert_eq!(col.len(), nbasis, "第{}个轨道的系数数量与Nbasis不一致", j + 1);
            for (i, val) in col.iter().enumerate() {
                cmat[(i, j)] = *val;
            }
        }
        (ato_atms, ato_shls, ato_syms, obt_engs, obt_occs, cmat)
    }

    /// 读取α和β电子数
    pub fn read_neles(&self) -> [usize; 2] {
        let nela = self.scalar("Naelec").expect("缺少Naelec字段").round() as usize;
        let nelb = self.scalar("Nbelec").expect("缺少Nbelec字段").round() as usize;
        [nela, nelb]
    }
}

impl Reader for MfnReader {
    fn get_atoms(&self) -> Atoms {
        let (atms, xyzs, crgs) = self.read_geome();
        let mut atoms = Atoms::new(&atms, &xyzs, &crgs);
        atoms.neles = self.read_neles();
        atoms
    }

    fn get_basis(&self) -> Option<Basis> {
        if !self.marks.contains_key("Shell types") {
            return None;
        }
        let (_, _, _, datas, forms) = self.read_basis();
        Some(shl_basis(&datas, &forms))
    }

    fn get_coefs(&self) -> Option<Coefs> {
        if !self.marks.contains_key("Shell types") {
            return None;
        }
        let (ato_atms, ato_shls, ato_syms, obt_engs, obt_occs, cmat) = self.read_coefs();
        let coefs = Coefs::new(ato_atms, ato_shls, ato_syms, obt_engs, obt_occs, cmat);
        Some(coefs)
    }
}
