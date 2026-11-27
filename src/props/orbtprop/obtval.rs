use std::collections::HashMap;

use ndarray::Axis;

use crate::base::Mole;
use crate::maths;
use crate::props::atomprop;
use crate::props::orbtprop;
use crate::props::orbtprop::mocv::MocvDest;
use crate::rsint;
use crate::utils;
/// 分子轨道属性计算器，计算轨道电子数和相关性质
pub struct Calculator<'m> {
    /// 分子信息
    pub mole: &'m Mole,
    /// 参与计算的原子索引
    pub atms: &'m Vec<usize>,
}

impl<'m> Calculator<'m> {
    /// 创建一个新的 Calculator 实例
    pub fn new(mole: &'m Mole) -> Self {
        let atms = mole.atoms().idxs();
        Calculator { mole, atms }
    }

    // 计算分子轨道的pi电子数
    /// 使用 POCV 方法计算每个分子轨道的 π 电子数
    pub fn pi_eles_pocv(&self) -> Vec<f64> {
        let mut caler = orbtprop::obtmat::Calculator::new(self.mole);
        caler.vobt = true;
        caler.atms = self.atms;
        let cmat = self.mole.coefs().get_cmat("raw", self.mole);
        let (_, cmat_p) = caler.pi_pocv().expect("POCV 投影失败");
        let smat = self.mole.basis().get_smat("raw", self.mole);
        let k = if self.mole.is_open() { 1.0 } else { 2.0 };
        let cmat_s = cmat - &cmat_p;
        let nobt = cmat.ncols();
        let mut eles = Vec::new();
        for i in 0..nobt {
            let c_pi = cmat_p.column(i);
            let c_si = cmat_s.column(i);
            let pp = c_pi.dot(&smat.dot(&c_pi)); // c_πᵀ S c_π
            let ss = c_si.dot(&smat.dot(&c_si)); // c_σᵀ S c_σ
            let sum = pp + ss;
            if sum.abs() < f64::EPSILON {
                eles.push(0.0);
            } else {
                eles.push(pp / sum * k);
            }
        }
        eles
    }

    /// 使用 MOCV 方法计算每个分子轨道的 π 电子数
    pub fn pi_eles_mocv(&self) -> Vec<f64> {
        let mut caler = orbtprop::obtmat::Calculator::new(self.mole);
        caler.vobt = true;
        caler.atms = self.atms;
        let dest = MocvDest::PiProperty;
        let (_, _, cmat_p) = caler.pi_mocv(&dest);
        let cmat = self.mole.coefs().get_cmat("sph", self.mole);
        let smat = self.mole.basis().get_smat("sph", self.mole);
        let k = if self.mole.is_open() { 1.0 } else { 2.0 };
        let cmat_s = cmat - &cmat_p;
        let nobt = cmat.ncols();
        let mut eles = Vec::new();
        for i in 0..nobt {
            let c_pi = cmat_p.column(i);
            let c_si = cmat_s.column(i);
            let pp = c_pi.dot(&smat.dot(&c_pi)); // c_πᵀ S c_π
            let ss = c_si.dot(&smat.dot(&c_si)); // c_σᵀ S c_σ
            let sum = pp + ss;
            if sum.abs() < f64::EPSILON {
                eles.push(0.0);
            } else {
                eles.push(pp / sum * k);
            }
        }
        eles
    }

    /// 计算每个分子轨道的电子数
    pub fn eles(&self) -> Vec<f64> {
        // let mut caler = orbtprop::obtmat::Calculator::new(self.mole);
        // caler.vobt = true;
        // let atms = self.mole.atoms().idxs();
        // let (dirs, cmat) = caler.pi_pocv(atms);
        let cmat = self.mole.coefs().get_cmat("sph", self.mole);
        let smat = self.mole.basis().get_smat("sph", self.mole);
        let mut eles = Vec::new();
        let nobt = cmat.ncols();
        let k = if self.mole.is_open() { 1.0 } else { 2.0 };
        for i in 0..nobt {
            let cvec = cmat.column(i).insert_axis(Axis(1));
            let dmat = cvec.dot(&cvec.t());
            let qmat = dmat * smat;
            eles.push(qmat.sum() * k);
        }
        eles
    }
}
