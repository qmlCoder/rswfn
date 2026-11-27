//! 键级（Mayer、Wiberg 等）计算。

use faer::stats::DistributionExt;
use ndarray::{Array2, Axis, s};
use std::collections::HashMap;

use crate::base::{Mole, Stm, Vec3};
use crate::maths::calc_dmat;
use crate::props::orbtprop::mocv;
use crate::rsint::RsInt;
// use crate::props::orbtprop::obtmat::Mocv;
use crate::maths;
use crate::props::orbtprop::mocv::{MocvDest, mocv_atos};
use crate::props::{atomprop, orbtprop};

/// 键级计算器，支持 Mayer 和 Wiberg 键级分析
pub struct Calculator<'m> {
    /// 分子引用
    pub mole: &'m Mole,
    /// 全部分子的原子索引列表
    pub atms: &'m Vec<usize>,
}

impl<'m> Calculator<'m> {
    /// 创建一个新的 Calculator 实例
    pub fn new(mole: &'m Mole) -> Self {
        let atms = mole.atoms().idxs();
        Calculator { mole, atms }
    }

    /// 计算 Mayer 键级矩阵
    pub fn mayer(&self) -> Array2<f64> {
        self.order("mayer")
    }

    /// 计算 Wiberg 键级矩阵（Löwdin 正交化基组下）
    pub fn wiberg(&self) -> Array2<f64> {
        self.order("wiberg")
    }

    /// 计算指定类型的键级矩阵
    ///
    /// # 参数
    /// * `otype` - 键级类型，"mayer" 或 "wiberg"
    pub fn order(&self, otype: &str) -> Array2<f64> {
        let smat = self.mole.basis().get_smat("raw", self.mole);
        let dmat = self.mole.coefs().get_dmat("raw", self.mole);
        assert_eq!(smat.shape(), dmat.shape());
        let uplo = self.mole.basis().atm_uplo("raw");
        let omat_ato = match otype {
            "mayer" => mayer(dmat, smat),
            "wiberg" => wiberg(dmat, smat),
            _ => panic!("未知键级类型:{}", otype),
        };
        ato2atm(&omat_ato, &uplo)
    }

    /// 计算某个分子轨道对键级的贡献
    ///
    /// 对于占据轨道，删去该轨道后键级减少表示正贡献；
    /// 对于空轨道，添加该轨道后键级增加表示正贡献。
    pub fn order_obt(&self, otype: &str, iobt: usize) -> Array2<f64> {
        let smat = self.mole.basis().get_smat("raw", self.mole);
        let old_cmat = self.mole.coefs().get_cmat("raw", self.mole);
        let old_dmat = self.mole.coefs().get_dmat("raw", self.mole);
        let nato = old_cmat.nrows();
        let nobt = old_cmat.ncols();
        let obt_occs = self.mole.coefs().obt_occs.clone();
        let new_dmat = if obt_occs[iobt] != 0.0 {
            let mut new_cmat = old_cmat.clone();
            new_cmat.slice_mut(s![.., iobt]).fill(0.0);
            maths::calc_dmat(&new_cmat, &obt_occs)
        } else {
            let mut obt_occs = obt_occs.clone();
            obt_occs[iobt] = if self.mole.is_open() { 1.0 } else { 2.0 };
            maths::calc_dmat(&old_cmat, &obt_occs)
        };
        let old_omat = order(&old_dmat, smat, otype);
        let new_omat = order(&new_dmat, smat, otype);
        let omat_ato = if obt_occs[iobt] != 0.0 {
            old_omat - new_omat
        } else {
            new_omat - old_omat
        };
        let uplo = &self.mole.basis().atm_uplo("raw");
        ato2atm(&omat_ato, uplo)
    }

    /// 使用 POCV（投影轨道系数向量）方法计算键级
    pub fn pocv(
        &self,
        dirs: &HashMap<usize, Vec3>,
        keep_other_atm: bool,
        keep_other_sym: bool,
        otype: &str,
    ) -> Result<Array2<f64>, String> {
        let obt_caler = orbtprop::obtmat::Calculator::new(self.mole);
        let cmat_pocv = obt_caler.pocv(dirs, keep_other_atm, keep_other_sym)?;

        let mut mole_pocv = self.mole.clone();
        mole_pocv.coefs.as_mut().unwrap().set_cmat("raw", &cmat_pocv, self.mole);
        let caler = Calculator::new(&mole_pocv);
        Ok(caler.order(otype))
    }

    /// 使用 MOCV（分子轨道系数向量）方法计算键级
    pub fn mocv(&self, stms: &HashMap<usize, Stm>, atos: &Vec<usize>, otype: &str) -> Array2<f64> {
        let obt_caler = orbtprop::obtmat::Calculator::new(self.mole);
        let cmat_mocv = obt_caler.mocv(stms, atos);
        let mut mole_mocv = self.mole.clone();
        mole_mocv.coefs.as_mut().unwrap().set_cmat("sph", &cmat_mocv, self.mole);
        let caler = Calculator::new(&mole_mocv);
        caler.order(otype)
    }

    /// 使用 POCV 方法计算 π 键级
    pub fn pi_pocv(&self, otype: &str) -> Result<(HashMap<usize, Vec3>, Array2<f64>), String> {
        let mut obt_caler = orbtprop::obtmat::Calculator::new(self.mole);
        obt_caler.atms = self.atms;
        let (dirs, cmat_pocv) = obt_caler.pi_pocv()?;

        let mut mole_pocv = self.mole.clone();
        mole_pocv.coefs.as_mut().unwrap().set_cmat("raw", &cmat_pocv, self.mole);

        let caler = Calculator::new(&mole_pocv);
        let mut omat = caler.order(otype);
        omat = omat.abs().sqrt();
        Ok((dirs, omat))
    }

    /// 使用 MOCV 方法计算 π 键级
    pub fn pi_mocv(&self, otype: &str) -> (HashMap<usize, Stm>, Vec<usize>, Array2<f64>) {
        let obt_caler = orbtprop::obtmat::Calculator::new(self.mole);
        let dest = MocvDest::PiProperty;
        let (stms, atos, cmat_mocv) = obt_caler.pi_mocv(&dest);
        let mut mole_mocv = self.mole.clone();
        mole_mocv.coefs.as_mut().unwrap().set_cmat("sph", &cmat_mocv, self.mole);
        let caler = Calculator::new(&mole_mocv);
        let mut omat = caler.order(otype);
        omat = omat.abs().sqrt(); // 取绝对值再开根号
        (stms, atos, omat)
    }

    /// 键级分解 获得σ,pi和δ键级
    pub fn decom(&self, ai: usize, aj: usize, otype: &str) -> Vec<f64> {
        let ez = self.mole.atoms.vec(ai, aj).normalize();
        let er = Vec3::rand();
        let ey = ez.cross(&er);
        let ex = ey.cross(&ez);
        let stm = Stm::new(ex.into_arr(), ey.into_arr(), ez.into_arr());
        let mut stms = HashMap::new();
        stms.insert(ai, stm.clone());
        stms.insert(aj, stm.clone());
        let dests = ["OrderDecom-SI", "OrderDecom-PI", "OrderDecom-DE", "OrderDecom-PH"];
        let sqrts = [false, false, false, false];
        let dests: Vec<MocvDest> = dests.iter().map(|dest| MocvDest::from_str(*dest).unwrap()).collect();
        let vals: Vec<f64> = dests
            .iter()
            .enumerate()
            .map(|(i, dest)| {
                let atos = mocv_atos(self.mole, dest);
                let omat = self.mocv(&stms, &atos, otype);
                let val = omat[(ai, aj)];
                if sqrts[i] { val.sqrt() } else { val }
            })
            .collect();
        vals
    }

    /// 计算束缚键级 使用pocv或mocv
    pub fn bound(
        &self,
        atm: usize,
        dir: &Vec3,
        nebs: &Vec<usize>,
        otype: &str,
        mtype: &str,
    ) -> Result<Vec<f64>, String> {
        let mut info = format!("计算束缚键级:{} {} {}\n", atm + 1, dir, otype);

        let omat = match mtype {
            "pocv" => {
                let mut dirs = HashMap::new();
                dirs.insert(atm, dir.clone());
                self.pocv(&dirs, true, true, otype)?
            }
            "mocv" => {
                let mut stms = HashMap::new();
                let stm = Stm::rand_ez(dir.into_arr());
                stms.insert(atm, stm);
                let natm = self.mole.atoms.len();
                for i in 0..natm {
                    if i != atm {
                        stms.insert(i, Stm::blank());
                    }
                }
                let ato_atms = self.mole.basis().atos().atms("sph");
                let ato_syms = self.mole.basis().atos().syms("sph");
                let dest = mocv::MocvDest::from_str("OrderDecom-SI").unwrap();
                let keep_syms = mocv::mocv_syms(&dest, 0);
                let mut atos = vec![]; // 局部坐标系下保留的原子轨道
                for (i, idx) in ato_atms.iter().enumerate() {
                    if *idx == atm {
                        if keep_syms.contains(&ato_syms[i]) {
                            atos.push(i);
                        }
                    } else {
                        atos.push(i);
                    }
                }
                self.mocv(&stms, &atos, otype)
            }
            _ => panic!("方法只可以使用pocv或mocv"),
        };
        let mut orders = Vec::new();
        for idx in nebs.iter() {
            let order = omat[(atm as usize, *idx as usize)];
            info.push_str(&format!("order.bound:{}->{}:{:>10.4}\n", atm + 1, idx + 1, order));
            orders.push(order);
        }
        log::info!("{}", info);
        Ok(orders)
    }
}

/// 根据指定类型计算键级矩阵
pub fn order(dmat: &Array2<f64>, smat: &Array2<f64>, otype: &str) -> Array2<f64> {
    match otype {
        "mayer" => mayer(dmat, smat),
        "wiberg" => wiberg(dmat, smat),
        _ => mayer(dmat, smat),
    }
}

/// Mayer 键级分析
pub fn mayer(dmat: &Array2<f64>, smat: &Array2<f64>) -> Array2<f64> {
    let psmat = dmat.dot(smat);
    &psmat * &psmat.t()
}

/// Wiberg bond order analysis in Löwdin orthogonalized basis
pub fn wiberg(dmat: &Array2<f64>, smat: &Array2<f64>) -> Array2<f64> {
    let smat_sqrt = maths::matrix_phalf(smat);
    let pmat_lowd = smat_sqrt.dot(dmat).dot(&smat_sqrt);
    &pmat_lowd * &pmat_lowd
}

/// 将原子轨道（AO）键级矩阵约化为原子（ATM）键级矩阵
pub fn ato2atm(omat_ato: &Array2<f64>, uplo: &Vec<[usize; 2]>) -> Array2<f64> {
    let natm = uplo.len();
    let mut omat_atm = Array2::zeros((natm, natm));
    for i in 0..natm {
        let [ui, li] = uplo[i];
        for j in 0..natm {
            let [uj, lj] = uplo[j];
            let order = omat_ato.slice(s![ui..li, uj..lj]).sum();
            omat_atm[(i, j)] = order;
        }
    }
    omat_atm
}
