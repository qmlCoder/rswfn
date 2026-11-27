use ndarray::{Array2, s};

use crate::base::Mole;
use crate::datas::FACS;
use crate::datas::btrans;
use crate::datas::btrans::ang2syms;
use crate::datas::btrans::sym2lmn;
use crate::rsint::RsInt;
use core::fmt;
use std::collections::HashMap;
use std::f64::consts::PI;
use std::sync::Arc;
use std::sync::Mutex;
use std::sync::OnceLock;

/// 单个原子轨道（基组中的一个具体函数，如 s、px、dxy 等）
#[derive(Debug, Clone)]
pub struct Ato {
    /// 所属原子的索引
    pub atm: usize,
    /// 所属壳层的索引
    pub shl: usize,
    /// 角动量量子数 l
    pub ang: usize,
    /// 轨道符号（如 s、px、dxy）
    pub sym: String,
    /// 定域性（0 表示高指数定域，1 表示弥散）
    pub loc: usize,
    /// 轨道形式（"raw", "sph", "car"）
    pub form: String,
    /// 高斯指数列表
    pub alps: Vec<f64>,
    /// 收缩系数列表（与 `alps` 一一对应）
    pub coes: Vec<f64>,
    coes_norm: OnceLock<Vec<f64>>,
    lmn: OnceLock<[usize; 3]>,
}

impl Ato {
    /// 创建一个新的原子轨道对象
    ///
    /// # Parameters
    /// - `atm`: 原子索引
    /// - `shl`: 壳层索引
    /// - `ang`: 角动量
    /// - `sym`: 轨道符号（如 s, px, py, pz 等）
    /// - `loc`: 定域性（0 表示高指数定域，1 表示弥散）
    /// - `form`: 轨道形式（"raw", "sph", "car"）
    /// - `alps`: 指数列表
    /// - `coes`: 收缩系数列表
    pub fn new(atm: usize, shl: usize, ang: usize, sym: String, loc: usize, form: String, alps: &Vec<f64>, coes: &Vec<f64>) -> Self {
        let alps = alps.clone();
        let coes = coes.clone();
        let coes_norm = OnceLock::new();
        let lmn = OnceLock::new();
        Self {
            atm,
            shl,
            ang,
            sym,
            loc,
            form,
            alps,
            coes,
            coes_norm,
            lmn,
        }
    }
    /// 返回归一化后的收缩系数（惰性初始化）
    ///
    /// 对原始收缩系数进行归一化处理，考虑 (l, m, n) 阶乘因子
    pub fn coes_norm(&self) -> &Vec<f64> {
        self.coes_norm.get_or_init(|| {
            let [l, m, n] = self.lmn();
            let mut coes = self.coes.clone();
            for i in 0..coes.len() {
                let fac = FACS[*l] * FACS[*m] * FACS[*n];
                let alp = self.alps[i];
                let coe = self.coes[i];
                let nm = (2.0 * alp / PI).powf(0.75) * f64::sqrt((4.0 * alp).powi(self.ang as i32) / fac);
                coes[i] = nm * coe;
            }
            coes
        })
    }

    /// 返回轨道对称性对应的 (l, m, n) 指数（惰性初始化）
    pub fn lmn(&self) -> &[usize; 3] {
        self.lmn.get_or_init(|| sym2lmn(&self.sym))
    }
}

/// 原子轨道集合，包含 raw、sph、car 三种形式
#[derive(Debug, Clone)]
pub struct Atos {
    /// 原始形式的原子轨道
    pub raws: Vec<Ato>,
    /// 球谐形式的原子轨道
    pub sphs: Vec<Ato>,
    /// 笛卡尔形式的原子轨道
    pub cars: Vec<Ato>,
}

impl Atos {
    /// 返回指定形式的原子轨道数量
    ///
    /// # Parameters
    /// - `form`: 轨道形式（"raw", "sph", "car"）
    pub fn len(&self, form: &str) -> usize {
        match form {
            "raw" => self.raws.len(),
            "sph" => self.sphs.len(),
            "car" => self.cars.len(),
            _ => panic!("form error"),
        }
    }

    /// 返回指定形式的原子索引列表
    ///
    /// # Parameters
    /// - `form`: 轨道形式（"raw", "sph", "car"）
    pub fn atms(&self, form: &str) -> Vec<usize> {
        match form {
            "raw" => self.raws.iter().map(|ato| ato.atm).collect(),
            "sph" => self.sphs.iter().map(|ato| ato.atm).collect(),
            "car" => self.cars.iter().map(|ato| ato.atm).collect(),
            _ => panic!("form error"),
        }
    }

    /// 返回指定形式的壳层索引列表
    ///
    /// # Parameters
    /// - `form`: 轨道形式（"raw", "sph", "car"）
    pub fn shls(&self, form: &str) -> Vec<usize> {
        match form {
            "raw" => self.raws.iter().map(|ato| ato.shl).collect(),
            "sph" => self.sphs.iter().map(|ato| ato.shl).collect(),
            "car" => self.cars.iter().map(|ato| ato.shl).collect(),
            _ => panic!("form error"),
        }
    }

    /// 返回指定形式的轨道符号列表
    ///
    /// # Parameters
    /// - `form`: 轨道形式（"raw", "sph", "car"）
    pub fn syms(&self, form: &str) -> Vec<String> {
        match form {
            "raw" => self.raws.iter().map(|ato| ato.sym.clone()).collect(),
            "sph" => self.sphs.iter().map(|ato| ato.sym.clone()).collect(),
            "car" => self.cars.iter().map(|ato| ato.sym.clone()).collect(),
            _ => panic!("form error"),
        }
    }

    /// 返回指定形式的定域性索引列表
    ///
    /// # Parameters
    /// - `form`: 轨道形式（"raw", "sph", "car"）
    pub fn locs(&self, form: &str) -> Vec<usize> {
        match form {
            "raw" => self.raws.iter().map(|ato| ato.loc).collect(),
            "sph" => self.sphs.iter().map(|ato| ato.loc).collect(),
            "car" => self.cars.iter().map(|ato| ato.loc).collect(),
            _ => panic!("form error"),
        }
    }
}

/// 基组中的一个壳层（一组共享收缩系数的原始高斯函数）
#[derive(Debug, Clone)]
pub struct Shl {
    /// 所属原子的索引
    pub atm: usize,
    /// 壳层索引
    pub shl: usize,
    /// 角动量量子数 l
    pub ang: usize,
    /// 函数形式（"sph" 球谐 / "car" 笛卡尔）
    pub form: String,
    /// 高斯指数列表
    pub alps: Vec<f64>,
    /// 收缩系数列表（与 `alps` 一一对应）
    pub coes: Vec<f64>,
}

impl Shl {
    /// 创建一个新的壳层
    ///
    /// # Parameters
    /// - `atm`: 原子索引
    /// - `shl`: 壳层索引
    /// - `ang`: 角动量
    /// - `form`: 函数形式（sph/car)
    /// - `alps`: 指数列表
    /// - `coes`: 收缩系数列表
    pub fn new(atm: usize, shl: usize, ang: usize, form: String, alps: Vec<f64>, coes: Vec<f64>) -> Self {
        Self { atm, shl, ang, form, alps, coes }
    }
    /// 返回归一化后的收缩系数
    pub fn coes_norm(&self) -> Vec<f64> {
        let mut coes = self.coes.clone();
        for i in 0..coes.len() {
            let alp = self.alps[i];
            let coe = self.coes[i];
            let nm = (2.0 * alp / PI).powf(0.75) * f64::sqrt((4.0 * alp).powi(self.ang as i32));
            coes[i] = nm * coe;
        }
        coes
    }

    /// 判断壳层的定域性
    ///
    /// # Returns
    /// - `0`: 高指数定域（最大指数 > 10.0）
    /// - `1`: 弥散（最大指数 <= 10.0）
    pub fn loc(&self) -> usize {
        let alp_max = self.alps.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
        if *alp_max > 10.0 { 0 } else { 1 }
    }
}

/// 壳层集合，按原始顺序存储基组的所有壳层
#[derive(Debug, Clone)]
pub struct Shls {
    /// 壳层列表
    pub shls: Vec<Shl>,
    /// 各壳层在 raw/sph/car 三种形式下的起止索引缓存，依次为 [raw, sph, car]
    pub uplo: OnceLock<[Vec<[usize; 2]>; 3]>,
}

const ANG_NUM_SPH: [usize; 6] = [1, 3, 5, 7, 9, 11];
const ANG_NUM_CAR: [usize; 6] = [1, 3, 6, 10, 15, 21];

impl Shls {
    /// 创建一个新的壳层集合
    ///
    /// # Parameters
    /// - `shls`: 壳层列表
    pub fn new(shls: Vec<Shl>) -> Shls {
        Shls { shls, uplo: OnceLock::new() }
    }

    /// 获取指定索引的壳层
    ///
    /// # Parameters
    /// - `i`: 壳层索引
    pub fn get(&self, i: usize) -> &Shl {
        &self.shls[i]
    }

    /// 返回壳层总数
    pub fn len(&self) -> usize {
        self.shls.len()
    }
    /// 获取每个壳层的起止索引范围（惰性初始化）
    ///
    /// # Parameters
    /// - `form`: 轨道形式（"raw", "sph", "car"）
    ///
    /// # Returns
    /// 每个壳层对应的 [起始索引, 结束索引) 列表
    pub fn uplo(&self, form: &str) -> &Vec<[usize; 2]> {
        let uplo = self.uplo.get_or_init(|| {
            let mut uplo_raw = vec![];
            let mut uplo_sph = vec![];
            let mut uplo_car = vec![];
            let mut idx_raw = 0;
            let mut idx_sph = 0;
            let mut idx_car = 0;
            for shl in self.shls.iter() {
                let ang = shl.ang;
                let num_raw = match shl.form.as_str() {
                    "car" => ANG_NUM_CAR[ang],
                    "sph" => ANG_NUM_SPH[ang],
                    _ => panic!("form error"),
                };
                let num_sph = ANG_NUM_SPH[ang];
                let num_car = ANG_NUM_CAR[ang];
                uplo_raw.push([idx_raw, idx_raw + num_raw]);
                uplo_sph.push([idx_sph, idx_sph + num_sph]);
                uplo_car.push([idx_car, idx_car + num_car]);
                idx_raw += num_raw;
                idx_sph += num_sph;
                idx_car += num_car;
            }
            [uplo_raw, uplo_sph, uplo_car]
        });
        match form {
            "raw" => &uplo[0],
            "sph" => &uplo[1],
            "car" => &uplo[2],
            _ => panic!("form error"),
        }
    }

    /// 返回指定形式的原子轨道总数
    ///
    /// # Parameters
    /// - `form`: 轨道形式（"raw", "sph", "car"）
    pub fn nato(&self, form: &str) -> usize {
        let uplo = self.uplo(form);
        let mut nato = 0;
        for [up, lo] in uplo.iter() {
            nato += lo - up;
        }
        nato
    }

    /// 按原子分组返回每个原子的起止索引范围
    ///
    /// # Parameters
    /// - `form`: 轨道形式（"raw", "sph", "car"）
    pub fn atm_uplo(&self, form: &str) -> Vec<[usize; 2]> {
        let uplo = self.uplo(form);
        let mut atm_uplo = Vec::new();
        let mut prev_atm = None;
        for (i, shl) in self.shls.iter().enumerate() {
            let [start, end] = uplo[i];
            if Some(shl.atm) != prev_atm {
                atm_uplo.push([start, end]);
                prev_atm = Some(shl.atm);
            } else {
                atm_uplo.last_mut().unwrap()[1] = end;
            }
        }
        atm_uplo
    }
}

/// 单个原始基函数（一个未收缩的高斯函数）
#[derive(Debug, Clone)]
pub struct Bas {
    /// 所属原子的索引
    pub atm: usize,
    /// 所属壳层的索引
    pub shl: usize,
    /// 角动量量子数 l
    pub ang: usize,
    /// 高斯指数
    pub alp: f64,
    /// 收缩系数
    pub coe: f64,
    // pub nlm: (usize, usize, isize), //n,l,m 兼容真实原子轨道，先这么写，后面应该要搞成枚举类型
}

impl Bas {
    /// 创建一个原始基函数
    ///
    /// # Parameters
    /// - `atm`: 原子索引
    /// - `shl`: 壳层索引
    /// - `ang`: 角动量
    /// - `alp`: 高斯指数
    /// - `coe`: 收缩系数
    pub fn new(atm: usize, shl: usize, ang: usize, alp: f64, coe: f64) -> Self {
        Self { atm, shl, ang, alp, coe }
    }
}
/// 原始基函数集合
#[derive(Debug, Clone)]
pub struct Bass {
    bass: Vec<Bas>,
}

/// 分子基组，包含壳层定义以及惰性构建的原子轨道、原始基函数和重叠积分
#[derive(Debug)]
pub struct Basis {
    /// 壳层集合，每个壳层记录 [atm, shl, ang, form, alps, coes]
    pub shls: Shls,
    /// 基组名称（若已知）
    pub name: Option<String>,
    /// 惰性构建的原子轨道集合（raw/sph/car 三种形式）
    pub atos: OnceLock<Atos>,
    /// 惰性构建的原始基函数集合
    pub bass: OnceLock<Bass>,
    /// 惰性计算的重叠积分矩阵，依次为 [raw, sph, car]
    pub smat: OnceLock<[Array2<f64>; 3]>,
}

impl Basis {
    /// 从壳层列表创建新的基组
    ///
    /// # Parameters
    /// - `shls`: 壳层列表
    pub fn new(shls: Vec<Shl>, name: Option<String>) -> Basis {
        let shls = Shls::new(shls);
        let atos = OnceLock::new();
        let bass = OnceLock::new();
        let smat = OnceLock::new();
        Basis { shls, name, atos, bass, smat }
    }

    /// 用新的壳层列表重新构建基组
    ///
    /// # Parameters
    /// - `shls`: 新的壳层列表
    pub fn build(&mut self, shls: Vec<Shl>) {
        self.shls = Shls::new(shls);
    }

    /// 判断基组是否为空（没有壳层）
    pub fn is_none(&self) -> bool {
        self.shls.len() == 0
    }

    /// 判断基组在积分计算中应当采用的轨道形式
    ///
    /// s、p 壳层的球谐与笛卡尔形式完全一致，因此只需看 l >= 2 的壳层：
    /// 只要存在一个笛卡尔壳层就返回 `"car"`，否则返回 `"sph"`。
    /// 对于 l >= 2 形式混合的基组（例如 6D 7F），积分只能统一按一种形式计算，
    /// 这也是 `unify_form` 存在的原因。
    pub fn int_form(&self) -> &'static str {
        let car = self.shls.shls.iter().any(|shl| shl.ang >= 2 && shl.form == "car");
        if car { "car" } else { "sph" }
    }

    /// 返回一个把 l >= 2 壳层形式统一为 `form` 的基组副本
    ///
    /// s、p 壳层保持原样（两种形式在数值和顺序上完全等价），
    /// 这样在 `form` 下算出的系数矩阵与该基组声明的壳层形式一一对应，
    /// 后续 `Coefs` 的形式转换（raw/sph/car）才不会错位。
    pub fn unify_form(&self, form: &str) -> Basis {
        let shls = self
            .shls
            .shls
            .iter()
            .map(|shl| {
                let mut shl = shl.clone();
                if shl.ang >= 2 {
                    shl.form = form.to_string();
                }
                shl
            })
            .collect();
        Basis::new(shls, self.name.clone())
    }

    /// 获取指定原子的所有壳层
    ///
    /// # Parameters
    /// - `atm`: 原子索引
    pub fn get_by_atm(&self, atm: usize) -> Vec<Shl> {
        let mut shls = vec![];
        let nshl = self.shls.len();
        for i in 0..nshl {
            let shl = self.shls.get(i);
            if shl.atm == atm {
                shls.push(shl.clone());
            }
        }
        shls
    }

    /// 返回所有原子轨道的唯一键值列表
    ///
    /// 格式为 "atm,shl,sym"
    pub fn ato_keys(&self) -> Vec<String> {
        let mut keys = vec![];
        let nshl = self.shls.len();
        for i in 0..nshl {
            let shl = self.shls.get(i);
            let syms = ang2syms(shl.ang, &shl.form);
            for sym in syms.iter() {
                let key = format!("{},{},{}", shl.atm, shl.shl, sym);
                if !keys.contains(&key) {
                    keys.push(key);
                }
            }
        }
        keys
    }

    /// 生成并返回原子轨道列表（惰性初始化）
    ///
    /// 包含 raw、sph、car 三种形式的原子轨道
    pub fn atos(&self) -> &Atos {
        let shls = &self.shls;
        let nshl = shls.len();
        self.atos.get_or_init(|| {
            let mut raws = vec![];
            let mut sphs = vec![];
            let mut cars = vec![];
            for i in 0..nshl {
                let atm = shls.get(i).atm;
                let shl = shls.get(i).shl;
                let ang = shls.get(i).ang;
                let loc = shls.get(i).loc();
                let form = &shls.get(i).form;
                let alps = &shls.get(i).alps;
                let coes = &shls.get(i).coes;
                let syms = ang2syms(ang, form);
                for sym in syms {
                    let ato = Ato::new(atm, shl, ang, sym.to_owned(), loc, form.to_owned(), alps, coes);
                    raws.push(ato);
                }
                let syms = ang2syms(ang, "sph");
                for sym in syms {
                    let ato = Ato::new(atm, shl, ang, sym.to_owned(), loc, form.to_owned(), alps, coes);
                    sphs.push(ato);
                }
                let syms = ang2syms(ang, "car");
                for sym in syms {
                    let ato = Ato::new(atm, shl, ang, sym.to_owned(), loc, form.to_owned(), alps, coes);
                    cars.push(ato);
                }
            }
            Atos { raws, sphs, cars }
        })
    }

    /// 返回原始基函数的数量
    pub fn nbas(&self) -> usize {
        let mut nbas = 0;
        let nshl = self.shls.len();
        for i in 0..nshl {
            let shl = self.shls.get(i);
            nbas += match shl.form.as_str() {
                "sph" => ANG_NUM_SPH[shl.ang],
                "car" => ANG_NUM_CAR[shl.ang],
                _ => 0,
            };
            nbas += self.shls.get(i).alps.len();
        }
        nbas
    }

    /// 返回壳层数量
    pub fn nshl(&self) -> usize {
        self.shls.len()
    }

    /// 判断两个基组信息是否相等
    ///
    /// 比较原子索引、壳层索引和角动量
    pub fn equal(&self, other: &Basis) -> bool {
        let datas_0 = &self.shls;
        let datas_1 = &other.shls;
        if datas_0.len() != datas_1.len() {
            log::error!("两个基组信息长度不相等");
            return false;
        }
        for i in 0..datas_0.len() {
            let each_1 = datas_0.get(i);
            let each_2 = datas_1.get(i);
            if each_1.atm != each_2.atm {
                log::error!("{},atm 不相等:{}!={}", i, each_1.atm, each_2.atm);
                return false;
            }
            if each_1.shl != each_2.shl {
                log::error!("{},shl 不相等:{}!={}", i, each_1.shl, each_2.shl);
                return false;
            }
            if each_1.ang != each_2.ang {
                log::error!("{},ang 不相等:{}!={}", i, each_1.ang, each_2.ang);
                return false;
            }
        }
        return true;
    }

    /// 获取指定形式的重叠积分矩阵（惰性初始化）
    ///
    /// # Parameters
    /// - `form`: 轨道形式（"raw", "sph", "car"）
    /// - `mole`: 分子信息，提供原子结构用于积分计算
    pub fn get_smat(&self, form: &str, mole: &Mole) -> &Array2<f64> {
        let smats = self.smat.get_or_init(|| {
            let int_caler = RsInt::new(mole);
            let smat_car = int_caler.smat("car");
            let smat_sph = int_caler.mat1e_sph(&smat_car);
            let smat_raw = int_caler.mat1e_raw(&smat_car);
            [smat_raw, smat_sph, smat_car]
        });
        match form {
            "raw" => &smats[0],
            "sph" => &smats[1],
            "car" => &smats[2],
            _ => panic!("Invalid form"),
        }
    }

    /// 返回按原子分组的起止索引范围
    ///
    /// # Parameters
    /// - `form`: 轨道形式（"raw", "sph", "car"）
    pub fn atm_uplo(&self, form: &str) -> Vec<[usize; 2]> {
        self.shls.atm_uplo(form)
    }
}

impl Clone for Basis {
    fn clone(&self) -> Self {
        Self {
            shls: self.shls.clone(),
            name: self.name.clone(),
            atos: self.atos.clone(),
            bass: self.bass.clone(),
            smat: self.smat.clone(),
        }
    }
}

impl fmt::Display for Basis {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "{:<5}{:<5}{:<5}{:<5}{:>10}", "idx", "atm", "shl", "ang", "form")?;
        let nshl = self.shls.len();
        for i in 0..nshl {
            let shl = self.shls.get(i);
            let alps = &shl.alps;
            let coes = &shl.coes;
            writeln!(f, "{:<5}{:<5}{:<5}{:<5}{:>10}({})", i, shl.atm, shl.shl, shl.ang, shl.form, alps.len())?;
            for i in 0..alps.len() {
                writeln!(f, "     {:>12.6E}  {:>12.6E}", alps[i], coes[i])?;
            }
        }
        Ok(())
    }
}

unsafe impl Send for Basis {}
unsafe impl Sync for Basis {}
