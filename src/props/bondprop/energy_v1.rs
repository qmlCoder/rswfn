//! 键能量分解的实现（v1）。
//!
//! 只保留两个原子的原子轨道系数得到键密度，按 `D = D_bond + D_rest` 做可加分解。
//! 对外接口见 [`super::energy`]。

use ndarray::{Array2, Axis};

use crate::base::Mole;
use crate::props::moleprop::energy as mole_energy;
use crate::qmcal::ks;
use crate::rsint::{RsInt, jk_mat};
use crate::rsxcs::XcFunc;

use super::energy::BondEnergy;

/// 校验原子对、准备泛函/基组，并计算 bond / rest / cross 分解
pub(super) fn bond_parts(mole: &Mole, atm1: usize, atm2: usize) -> Result<BondEnergy, String> {
    let natm = mole.atoms().len();
    if atm1 >= natm || atm2 >= natm {
        return Err(format!("原子序号越界: {atm1}, {atm2}（分子共有 {natm} 个原子）"));
    }
    if atm1 == atm2 {
        return Err("键能量需要两个不同的原子".to_string());
    }
    let (owned_mole, func) = mole_energy::Calculator::new(mole).prepare_log()?;
    let mole = owned_mole.as_ref().map_or(mole, |mole| mole);
    Ok(decompose(mole, func, atm1, atm2))
}

/// 由原子对构造键密度并计算能量分解
fn decompose(mole: &Mole, func: XcFunc, atm1: usize, atm2: usize) -> BondEnergy {
    let form = mole.basis().int_form();
    let occ = mole.coefs().obt_idxs(true).clone();

    // 完整密度矩阵 D = 2·C_occ·C_occᵀ
    let cmat_full = mole.coefs().get_cmat(form, mole);
    let coef_occ = cmat_full.select(Axis(1), &occ);
    let dmat = coef_occ.dot(&coef_occ.t()) * 2.0;

    // 只保留两个原子的原子轨道系数，得到键密度 D_bond = 2·C_bond·C_bondᵀ
    let ato_atms = mole.basis().atos().atms(form);
    let mut coef_bond = Array2::zeros(coef_occ.raw_dim());
    for (ato, atm) in ato_atms.iter().enumerate() {
        if *atm == atm1 || *atm == atm2 {
            coef_bond.row_mut(ato).assign(&coef_occ.row(ato));
        }
    }
    let dmat_bond = coef_bond.dot(&coef_bond.t()) * 2.0;
    // 补密度（保证 D = D_bond + D_rest 精确成立）
    let dmat_rest = &dmat - &dmat_bond;

    // 单/双电子积分（form 形式）
    let int_caler = RsInt::new(mole);
    let [_smat, tmat, vmat] = int_caler.mat1e([1, 1, 1], form);
    let hmat = &tmat + &vmat;
    let emat = int_caler.mat2e(form);

    // J、K 对 D 线性：分别对 D_bond、D_rest 求，再相加得总
    let (j_bond, k_bond) = jk_mat(&dmat_bond, &emat);
    let (j_rest, k_rest) = jk_mat(&dmat_rest, &emat);
    let j_tot = &j_bond + &j_rest;
    let k_tot = &k_bond + &k_rest;

    // 交换相关能（数值格点积分）：键、其余、总
    let rks = ks::RKS::new(mole, func);
    let xc_bond = rks.xc_eng(&dmat_bond);
    let xc_rest = rks.xc_eng(&dmat_rest);
    let xc_tot = rks.xc_eng(&dmat);
    let e_nuc = rks.nuc_eng();

    let ax = func.exact_x();
    let tr = |a: &Array2<f64>, b: &Array2<f64>| (a * b).sum();

    let e_bond = tr(&dmat_bond, &hmat) + 0.5 * tr(&dmat_bond, &j_bond) - ax * 0.25 * tr(&dmat_bond, &k_bond) + xc_bond;
    let e_rest = tr(&dmat_rest, &hmat) + 0.5 * tr(&dmat_rest, &j_rest) - ax * 0.25 * tr(&dmat_rest, &k_rest) + xc_rest;
    let e_cross_coul = tr(&dmat_bond, &j_rest);
    let e_cross_exch = -ax * 0.5 * tr(&dmat_bond, &k_rest);
    let e_cross_xc = xc_tot - xc_bond - xc_rest;
    let e_cross = e_cross_coul + e_cross_exch + e_cross_xc;
    let e_elec = tr(&dmat, &hmat) + 0.5 * tr(&dmat, &j_tot) - ax * 0.25 * tr(&dmat, &k_tot) + xc_tot;

    BondEnergy {
        e_bond,
        e_rest,
        e_cross,
        e_cross_coul,
        e_cross_exch,
        e_cross_xc,
        e_elec,
        e_nuc,
        e_total: e_elec + e_nuc,
    }
}
