/*分子变体 */

use ndarray::{Array1, Array2, s};

use crate::base::Mole;
use crate::props::moleprop;
use crate::props::orbtprop;

/// 新分子构造器，支持自旋分离和 π 电子投影
pub struct Calculator<'m> {
    /// 分子信息
    pub mole: &'m Mole,
    /// 参与计算的原子索引
    pub atms: &'m Vec<usize>,
}

impl<'m> Calculator<'m> {
    /// 创建一个新的 Calculator 实例
    pub fn new(mole: &'m Mole) -> Calculator {
        let atms = mole.atoms().idxs();
        Calculator { mole, atms }
    }

    /// 将分子分离为 α 和 β 自旋两个独立的分子对象
    pub fn mole_spin(&self) -> (Mole, Mole) {
        let cmat = self.mole.coefs().get_cmat("raw", self.mole);
        let [nato, nobt] = cmat.shape().try_into().expect("获取矩阵形状失败");
        let occs = &self.mole.coefs().obt_occs;
        if self.mole.is_open() {
            let mut mole_a = self.mole.clone();
            let mut mole_b = self.mole.clone();

            let mut cmat_a = cmat.clone();
            let mut cmat_b = cmat.clone();

            cmat_a.slice_mut(s![.., nobt / 2..]).fill(0.0);
            cmat_b.slice_mut(s![.., ..nobt / 2]).fill(0.0);

            mole_a.coefs.as_mut().unwrap().set_cmat("raw", &cmat_a, self.mole);
            mole_b.coefs.as_mut().unwrap().set_cmat("raw", &cmat_b, self.mole);
            (mole_a, mole_b)
        } else {
            let occs_a: Vec<f64> = occs.iter().map(|occ| if *occ > 0.0 { 1.0 } else { 0.0 }).collect();
            let occs_b: Vec<f64> = occs.iter().map(|occ| if *occ > 0.0 { 1.0 } else { 0.0 }).collect();
            let mut mole_a = self.mole.clone();
            let mut mole_b = self.mole.clone();
            mole_a.coefs.as_mut().unwrap().obt_occs = occs_a;
            mole_b.coefs.as_mut().unwrap().obt_occs = occs_b;
            (mole_a, mole_b)
        }
    }

    /// 生成仅包含 π 电子成分的分子对象
    pub fn mole_pi_pocv(&self) -> Result<Mole, String> {
        let mut mole_caler = orbtprop::obtmat::Calculator::new(self.mole);
        mole_caler.atms = self.atms;
        let (_dirs, cmat_pocv) = mole_caler.pi_pocv()?;
        let mut mole_pocv = self.mole.clone();
        mole_pocv.coefs.as_mut().unwrap().set_cmat("raw", &cmat_pocv, self.mole);
        Ok(mole_pocv)
    }

    /// 生成仅包含 π 电子成分（基于 MOCV 方法）的分子对象
    pub fn mole_pi_mocv(&self) -> Result<Mole, String> {
        let mut mole_caler = orbtprop::obtmat::Calculator::new(self.mole);
        mole_caler.atms = self.atms;
        let dest = orbtprop::mocv::MocvDest::PiProperty;
        let (_, _, cmat_mocv) = mole_caler.pi_mocv(&dest);
        let mut mole_pocv = self.mole.clone();
        mole_pocv.coefs.as_mut().unwrap().set_cmat("sph", &cmat_mocv, self.mole);
        Ok(mole_pocv)
    }
}
