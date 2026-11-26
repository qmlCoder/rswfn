#![allow(non_snake_case)]
use std::collections::HashMap;

use ndarray::Array2;
use ndarray::{Array1, s};
use rayon::prelude::*;

use crate::props::atomprop;
use crate::datas::rdens;
use crate::props::gridprop;
use crate::rsint;

use crate::base::{Mole, Stm, Vec3};
use crate::maths;
use crate::props::moleprop;
use crate::props::orbtprop;
use crate::props::orbtprop::mocv::MocvDest;

/// 原子电荷计算器。
///
/// 提供多种原子电荷/电子布居分析方法，包括 Mulliken、Lowdin、Hirshfeld
/// 以及 POÇV、MOÇV 等轨道改造方法。
pub struct Calculator<'m> {
    pub mole: &'m Mole,
    pub atms: &'m Vec<usize>,
}

impl<'m> Calculator<'m> {
    /// 创建一个新的电荷计算器。
    ///
    /// # Parameters
    /// - `mole` — 分子引用
    pub fn new(mole: &'m Mole) -> Self {
        let atms = mole.atoms().idxs();
        Calculator { mole, atms }
    }

    /// 计算自旋密度布居。
    ///
    /// 仅适用于开壳层分子，分别计算 α 和 β 电子的布居并求差。
    ///
    /// # Parameters
    /// - `ctype` — 电荷类型: `"mulliken"` | `"lowdin"` | `"hirshfeld"`
    ///
    /// # Returns
    /// 每个原子 3 个值: `[α电子, β电子, α - β]`
    pub fn spin(&self, ctype: &str) -> Vec<[f64; 3]> {
        let natm = self.mole.atoms().len();
        let mole_caler = moleprop::variants::Calculator::new(self.mole);
        let (mole_a, mole_b) = mole_caler.mole_spin();

        let caler_a = Calculator::new(&mole_a);
        let caler_b = Calculator::new(&mole_b);

        let (eles_a, eles_b) = (caler_a.eles(ctype), caler_b.eles(ctype));

        let mut spin = vec![[0.0; 3]; natm];
        for i in 0..natm {
            spin[i][0] = eles_a[i];
            spin[i][1] = eles_b[i];
            spin[i][2] = eles_a[i] - eles_b[i];
        }
        spin
    }

    /// 计算 Mulliken 电子布居（电子数）。
    ///
    /// # Returns
    /// 每个原子的电子数
    pub fn mulliken(&self) -> Vec<f64> {
        let basis = self.mole.basis();
        let coefs = self.mole.coefs();
        let dmat = coefs.get_dmat("raw", self.mole);
        let smat = basis.get_smat("raw", self.mole); // 重叠矩阵
        let uplo = basis.atm_uplo("raw");
        mulliken(dmat, smat, &uplo)
    }

    /// 计算 Löwdin 电子布居（电子数）。
    ///
    /// # Returns
    /// 每个原子的电子数
    pub fn lowdin(&self) -> Vec<f64> {
        let dmat = self.mole.coefs().get_dmat("raw", self.mole);
        let smat = self.mole.basis().get_smat("raw", self.mole);
        let uplo = self.mole.basis().atm_uplo("raw");
        lowdin(dmat, smat, &uplo)
    }

    /// 计算 Hirshfeld 电子布居（电子数）。
    ///
    /// 基于格点积分和原子密度权重分配电子密度。
    ///
    /// # Returns
    /// 每个原子的电子数
    pub fn hirshfeld(&self) -> Vec<f64> {
        let gridClaer = gridprop::dftgrid::Calculator::new(self.mole);
        let densCaler = gridprop::density::Calculator::new(self.mole);
        let atoms = self.mole.atoms();
        let natm = atoms.len();
        let mut eles = vec![0.0; natm]; // 原子的电子数
        let (mole_grids, mole_weits) = gridClaer.mol_grids(75, 434);
        let ngrid = mole_grids.len();
        let mole_dens: Vec<f64> = mole_grids
            .par_iter()
            .map(|grid| densCaler.mol_rho_dm(grid, 0).val0)
            .collect();
        let pro_mole_dens: Vec<f64> = mole_grids.par_iter().map(|grid| densCaler.pro_rho(grid)).collect();
        let mut pro_atom_dens = vec![0.0; ngrid];
        let mut atom_weits = vec![0.0; ngrid];

        // 计算原子与格点之间的距离
        for i in 0..natm {
            let nuc = atoms.get(i as usize).nuc;
            let [ax, ay, az] = atoms.get(i as usize).xyz;
            pro_atom_dens = mole_grids
                .par_iter()
                .map(|grid| {
                    let [gx, gy, gz] = grid;
                    let [dx, dy, dz] = [gx - ax, gy - ay, gz - az];
                    let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                    rdens::get_rad_dens_v1(nuc, dist)
                })
                .collect();
            atom_weits = (0..ngrid)
                .map(|g| {
                    if pro_mole_dens[g] > 1e-12 {
                        pro_atom_dens[g] / pro_mole_dens[g]
                    } else {
                        0.0
                    }
                })
                .collect();
            for g in 0..ngrid {
                eles[i] += mole_dens[g] * mole_weits[g] * atom_weits[g];
            }

            log::info!("原子{:>3}  电子数:{:>10.4}", i, eles[i]);
        }
        eles
    }

    /// 按指定方法计算电子布居。
    ///
    /// # Parameters
    /// - `ctype` — 电荷类型: `"mulliken"` | `"lowdin"` | `"hirshfeld"`
    ///
    /// # Returns
    /// 每个原子的电子数
    pub fn eles(&self, ctype: &str) -> Vec<f64> {
        match ctype {
            "mulliken" => self.mulliken(),
            "lowdin" => self.lowdin(),
            "hirshfeld" => self.hirshfeld(),
            _ => panic!("未知的计算方法:{}", ctype),
        }
    }

    /// 计算指定分子轨道内的原子电子布居。
    ///
    /// # Parameters
    /// - `ctype` — 电荷类型
    /// - `obt` — 分子轨道索引
    ///
    /// # Returns
    /// 每个原子在该轨道中的电子数
    pub fn eles_obt(&self, ctype: &str, obt: usize) -> Vec<f64> {
        let smat = self.mole.basis().get_smat("raw", self.mole);
        let old_cmat = self.mole.coefs().get_cmat("raw", self.mole);
        let nato = old_cmat.nrows();
        let nobt = old_cmat.ncols();
        let mut new_cmat = Array2::zeros((nato, nobt));
        let nocc = if self.mole.is_open() { 1.0 } else { 2.0 };
        new_cmat.slice_mut(s![.., obt]).assign(&old_cmat.slice(s![.., obt]));
        let dmat = old_cmat.dot(&new_cmat.t()) * nocc;
        let uplo = self.mole.basis().atm_uplo("raw");
        mulliken(&dmat, smat, &uplo)
    }

    /// 计算经 POÇV（逐轨道正交化修正）方法改造后的电子布居。
    ///
    /// # Parameters
    /// - `dirs` — 原子方向映射
    /// - `keep_other_atm` — 是否保留其他原子的轨道
    /// - `keep_other_sym` — 是否保留其他对称性的轨道
    /// - `ctype` — 电荷类型
    ///
    /// # Returns
    /// 每个原子的电子数
    pub fn pocv(
        &self,
        dirs: &HashMap<usize, Vec3>,
        keep_other_atm: bool,
        keep_other_sym: bool,
        ctype: &str,
    ) -> Result<Vec<f64>, String> {
        let obt_caler = orbtprop::obtmat::Calculator::new(self.mole);
        let cmat_pocv = obt_caler.pocv(dirs, keep_other_atm, keep_other_sym)?;

        let mut mole_pocv = self.mole.clone();
        mole_pocv.coefs.as_mut().unwrap().set_cmat("raw", &cmat_pocv, self.mole);
        let ele_caler = Calculator::new(&mole_pocv);
        let eles = match ctype {
            "mulliken" => ele_caler.mulliken(),
            "lowdin" => ele_caler.lowdin(),
            "hirshfeld" => ele_caler.hirshfeld(),
            _ => panic!("未知的ctype"),
        };
        Ok(eles)
    }

    /// 计算经 MOÇV（分子轨道正交化修正）方法改造后的电子布居。
    ///
    /// # Parameters
    /// - `stms` — 局部坐标系映射
    /// - `atos` — 需要保留的原子轨道索引
    /// - `ctype` — 电荷类型
    ///
    /// # Returns
    /// 每个原子的电子数
    pub fn mocv(&self, stms: &HashMap<usize, Stm>, atos: &Vec<usize>, ctype: &str) -> Vec<f64> {
        let obt_caler = orbtprop::obtmat::Calculator::new(self.mole);
        let cmat_mocv = obt_caler.mocv(stms, atos);

        let mut mole_mocv = self.mole.clone();
        mole_mocv.coefs.as_mut().unwrap().set_cmat("sph", &cmat_mocv, self.mole);
        let ele_caler = Calculator::new(&mole_mocv);
        match ctype {
            "mulliken" => ele_caler.mulliken(),
            "lowdin" => ele_caler.lowdin(),
            "hirshfeld" => ele_caler.hirshfeld(),
            _ => panic!("未知的ctype"),
        }
    }

    /// 计算 π 轨道的 POÇV 电子布居。
    ///
    /// # Parameters
    /// - `ctype` — 电荷类型
    ///
    /// # Returns
    /// `(dirs, eles)`，其中 `dirs` 是 π 方向映射，`eles` 是每个原子的 π 电子数
    pub fn pi_pocv(&self, ctype: &str) -> Result<(HashMap<usize, Vec3>, Vec<f64>), String> {
        let mut obt_caler = orbtprop::obtmat::Calculator::new(self.mole);
        obt_caler.atms = self.atms;
        let (dirs, cmat_pocv) = obt_caler.pi_pocv()?;

        let mut mole_pocv = self.mole.clone();
        mole_pocv.coefs.as_mut().unwrap().set_cmat("raw", &cmat_pocv, self.mole);
        let ele_caler = Calculator::new(&mole_pocv);
        let eles = match ctype {
            "mulliken" => ele_caler.mulliken(),
            "lowdin" => ele_caler.lowdin(),
            "hirshfeld" => ele_caler.hirshfeld(),
            _ => panic!("未知的ctype"),
        };
        Ok((dirs, eles))
    }

    /// 计算 π 轨道的 MOÇV 电子布居。
    ///
    /// # Parameters
    /// - `ctype` — 电荷类型
    ///
    /// # Returns
    /// `(stms, atos, eles)`，分别是局部坐标系映射、原子轨道索引和每个原子的 π 电子数
    pub fn pi_mocv(&self, ctype: &str) -> Result<(HashMap<usize, Stm>, Vec<usize>, Vec<f64>), String> {
        let obt_caler = orbtprop::obtmat::Calculator::new(self.mole);
        let dest = MocvDest::PiProperty;
        let (stms, atos, cmat_mocv) = obt_caler.pi_mocv(&dest);
        let mut mole_mocv = self.mole.clone();
        mole_mocv.coefs.as_mut().unwrap().set_cmat("sph", &cmat_mocv, self.mole);
        let caler = Calculator::new(&mole_mocv);
        let eles = match ctype {
            "mulliken" => caler.mulliken(),
            "lowdin" => caler.lowdin(),
            "hirshfeld" => caler.hirshfeld(),
            _ => panic!("未知的ctype"),
        };
        Ok((stms, atos, eles))
    }

    pub fn pi_mocv_v2(&self, ctype: &str) -> Result<(HashMap<usize, Stm>, Vec<usize>, Vec<f64>), String> {
        let obt_caler = orbtprop::obtmat::Calculator::new(self.mole);
        let cmat = self.mole.coefs().get_cmat("sph", self.mole);
        let eles = self.eles(ctype);
        let dest = MocvDest::PiProperty;
        let (stms, atos, cmat_pi) = obt_caler.pi_mocv(&dest);
        let mut mole_pi = self.mole.clone();
        mole_pi.coefs.as_mut().unwrap().set_cmat("sph", &cmat_pi, self.mole);
        let caler_pi = Calculator::new(&mole_pi);
        let eles_pp = caler_pi.eles(ctype);

        let cmat_si = cmat - &cmat_pi;
        let mut mole_si = self.mole.clone();
        mole_si.coefs.as_mut().unwrap().set_cmat("sph", &cmat_si, self.mole);
        let caler_si = Calculator::new(&mole_si);
        let eles_ss = caler_si.eles(ctype);

        let natm = eles.len();
        let mut eles_pi = vec![0.0; natm];
        for i in 0..natm {
            let ele_pp = eles_pp[i];
            let ele_ss = eles_ss[i];
            let esum = ele_pp + ele_ss;
            eles_pi[i] = if esum.abs() < f64::EPSILON {
                0.0
            } else {
                ele_pp / esum * eles[i]
            };
        }

        Ok((stms, atos, eles_pi))
    }
}

/// 执行 Mulliken 布居分析的数值计算。
///
/// # Parameters
/// - `dmat` — 密度矩阵
/// - `smat` — 重叠矩阵
/// - `uplo` — 每个原子的轨道起止索引
///
/// # Returns
/// 每个原子的电子数
pub fn mulliken(dmat: &Array2<f64>, smat: &Array2<f64>, uplo: &Vec<[usize; 2]>) -> Vec<f64> {
    assert_eq!(dmat.shape(), smat.shape());
    let qmat = dmat * smat;
    qmat2eles(&qmat, uplo)
}

/// 计算lowdin 电荷，返回电子数量
/// cmat：分子轨道系数矩阵
/// smat：重叠矩阵
/// nums：每个原子的轨道数量
/// idxs：占据轨道的索引
/// atos：原子轨道是否保留
pub fn lowdin(dmat: &Array2<f64>, smat: &Array2<f64>, uplo: &Vec<[usize; 2]>) -> Vec<f64> {
    let nato = dmat.nrows();
    let smat_half = maths::matrix_phalf(&smat); // 重叠矩阵的平方根
    let qmat = smat_half.dot(dmat).dot(&smat_half);
    qmat2eles_lowdin(&qmat, uplo)
}

/// 将布居矩阵（Q 矩阵）转换为每个原子的电子数。
///
/// # Parameters
/// - `qmat` — 布居矩阵
/// - `uplo` — 每个原子的轨道起止索引
///
/// # Returns
/// 每个原子的电子数
pub fn qmat2eles(qmat: &Array2<f64>, uplo: &Vec<[usize; 2]>) -> Vec<f64> {
    let natm = uplo.len();
    let nato = qmat.nrows();
    let mut eles = vec![0.0; natm];
    for a in 0..natm {
        let [ua, la] = uplo[a];
        for i in ua..la {
            for j in 0..nato {
                eles[a] += qmat[(i, j)];
            }
        }
    }
    eles
}

/// Lowdin 布居只取对角元（在 S^{1/2} 正交化基下，非对角项不贡献）
fn qmat2eles_lowdin(qmat: &Array2<f64>, uplo: &Vec<[usize; 2]>) -> Vec<f64> {
    let natm = uplo.len();
    let mut eles = vec![0.0; natm];
    for a in 0..natm {
        let [ua, la] = uplo[a];
        for i in ua..la {
            eles[a] += qmat[(i, i)];
        }
    }
    eles
}

/// 执行 Hirshfeld 布居分析的数值计算（独立函数版本）。
///
/// # Parameters
/// - `atom_nucs` — 各原子的核电荷数
/// - `atom_xyzs` — 各原子的坐标
/// - `mole_grids` — 分子格点坐标
/// - `mole_weits` — 格点权重
/// - `mole_dens` — 格点处的分子电子密度
/// - `pro_mole_dens` — 格点处的游离原子密度之和
///
/// # Returns
/// 每个原子的电子数
pub fn hirshfeld(
    atom_nucs: Vec<usize>,
    atom_xyzs: Vec<[f64; 3]>,
    mole_grids: &Vec<[f64; 3]>,
    mole_weits: &Vec<f64>,
    mole_dens: &Vec<f64>,
    pro_mole_dens: &Vec<f64>,
) -> Vec<f64> {
    let natm = atom_nucs.len();
    let mut eles = vec![0.0; natm]; // 原子的电子数
    let ngrid = mole_grids.len();
    let mut pro_atom_dens = vec![0.0; ngrid];
    let mut atom_weits = vec![0.0; ngrid];

    // 计算原子与格点之间的距离
    for i in 0..natm {
        let nuc = atom_nucs[i];
        let [ax, ay, az] = atom_xyzs[i];
        pro_atom_dens = mole_grids
            .par_iter()
            .map(|grid| {
                let [gx, gy, gz] = grid;
                let [dx, dy, dz] = [gx - ax, gy - ay, gz - az];
                let dist = (dx * dx + dy * dy + dz * dz).sqrt();
                rdens::get_rad_dens_v1(nuc, dist)
            })
            .collect();
        atom_weits = (0..ngrid)
            .map(|g| {
                if pro_mole_dens[g] > 1e-12 {
                    pro_atom_dens[g] / pro_mole_dens[g]
                } else {
                    0.0
                }
            })
            .collect();
        for g in 0..ngrid {
            eles[i] += mole_dens[g] * mole_weits[g] * atom_weits[g];
        }

        log::info!("原子{:>3}  电子数:{:>10.4}", i, eles[i]);
    }
    eles
}
