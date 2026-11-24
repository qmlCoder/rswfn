use std::collections::HashMap;
use std::time::{Duration, Instant};

use ndarray::{Array1, Array2, s};

use crate::base::{Basis, Coefs, Mole};
// 实现HF方法，根据结构和基组信息，计算得到波函数
use crate::rsint::{Mat2e, RsInt};
use crate::{atomprop, maths, rsint};

use super::diis::Diis;

/// 闭壳层Hartree-Fock计算结果
pub struct RHFResult {
    /// 电子能量（电子动能、电子-核吸引、电子-电子排斥之和）
    pub ele_eng: f64,
    /// 核-核排斥能
    pub nuc_eng: f64,
    /// 是否自洽收敛
    pub converged: bool,
    /// 实际迭代次数
    pub niter: usize,
    /// 最后一次迭代的能量变化（Hartree）
    pub delta_eng: f64,
    /// 单电子积分耗时
    pub time_int1e: Duration,
    /// 双电子积分耗时
    pub time_int2e: Duration,
    /// SCF 迭代耗时
    pub time_scf: Duration,
    /// 分子轨道系数
    pub coefs: Coefs,
    /// 与系数形式一致、可直接用于后续分析的基组
    pub basis: Basis,
    /// 实际使用的积分形式（"sph" 或 "car"）
    pub form: String,
}

impl RHFResult {
    /// 分子总能量（电子能量 + 核排斥能）
    pub fn total_eng(&self) -> f64 {
        self.ele_eng + self.nuc_eng
    }
}

/// 闭壳层Hartree-Fock计算器
pub struct RHF<'m> {
    /// 分子结构引用
    pub mole: &'m Mole,
    /// 是否使用局部坐标系
    pub use_lstm: bool,
    /// 积分形式("sph"或"car")，由基组自动决定
    pub int_form: String,
    /// 最大SCF迭代次数
    pub max_iter: usize,
    /// 收敛阈值（相邻两次迭代的能量变化，Hartree）
    pub tol: f64,
}

/// 根据闭壳层占据轨道系数与电子积分矩阵计算电子能量
///
/// - `dmat`: 由占据轨道系数构造的密度矩阵，`D = 2 C_occ C_occᵀ`
/// - `hmat`: 核心哈密顿量（动能 + 电子-核吸引）
/// - `fmat`: 由单、双电子积分构造的 Fock 矩阵
///
/// 返回电子能量 `E = ½ Tr[D(H + F)]`。
pub fn rhf_ele_eng(dmat: &Array2<f64>, hmat: &Array2<f64>, fmat: &Array2<f64>) -> f64 {
    (dmat * (hmat + fmat)).sum() * 0.5
}

impl<'m> RHF<'m> {
    pub fn new(mole: &'m Mole) -> Self {
        Self {
            mole,
            use_lstm: false,
            int_form: mole.basis().int_form().to_string(),
            max_iter: 100,
            tol: 1e-8,
        }
    }

    /// 计算核-核排斥能
    pub fn nuc_eng(&self) -> f64 {
        let natm = self.mole.atoms.len();
        let nucs = self.mole.atoms.nucs();
        let mut eng = 0.0;
        for i in 0..natm.saturating_sub(1) {
            for j in i + 1..natm {
                eng += (nucs[i] * nucs[j]) as f64 / self.mole.atoms.dist(i, j);
            }
        }
        eng
    }

    // 计算闭壳层密度矩阵 D = 2 * C_occ C_occ^T
    fn calc_rdmat(&self, nele: usize, cmat: &Array2<f64>) -> Array2<f64> {
        let coef = cmat.slice(s![.., ..nele / 2]).to_owned();
        coef.dot(&coef.t()) * 2.0
    }

    // 由双电子积分和密度矩阵构造库伦矩阵 J 与交换矩阵 K
    fn calc_jk(&self, dmat: &Array2<f64>, emat: &Mat2e) -> (Array2<f64>, Array2<f64>) {
        rsint::jk_mat(dmat, emat)
    }

    /// 执行闭壳层HF自洽场计算
    ///
    /// 通过迭代求解Roothaan方程,收敛后返回能量、分子轨道系数和对应的基组。
    ///
    /// # 参数
    /// - `emat`: 可选的双电子积分,不提供时内部自动计算
    pub fn run(&self, emat: Option<Mat2e>) -> RHFResult {
        let mole = self.mole;
        let [nela, nelb] = mole.atoms.neles;
        assert_eq!(nela, nelb, "RHF 仅支持闭壳层体系：α电子数({nela}) ≠ β电子数({nelb})");
        let nele = nela + nelb;
        let nocc = nele / 2;
        assert!(nocc > 0, "分子没有电子，无法进行 HF 计算");

        let ints = RsInt::new(mole);
        let (mat1e, emat, time_int1e, time_int2e) = if self.use_lstm {
            println!("使用局部坐标系");
            let mut dir_caler = atomprop::direction::Calculator::new(mole);
            let mut stms = HashMap::new();
            for atm in 0..mole.atoms.len() {
                let stm = dir_caler.LCS(atm, None).expect("局部坐标系生成失败");
                stms.insert(atm, stm);
            }
            let t_int1e = Instant::now();
            let mat1e = ints.mat1e_lcs(&stms, [1, 1, 1]);
            let time_int1e = t_int1e.elapsed();
            let (mat2e, time_int2e) = match emat {
                Some(emat) => (emat, Duration::ZERO),
                None => {
                    let t_int2e = Instant::now();
                    let mat2e = ints.mat2e_lcs(&stms);
                    (mat2e, t_int2e.elapsed())
                }
            };
            (mat1e, mat2e, time_int1e, time_int2e)
        } else {
            let t_int1e = Instant::now();
            let mat1e = ints.mat1e([1, 1, 1], &self.int_form);
            let time_int1e = t_int1e.elapsed();
            let (mat2e, time_int2e) = match emat {
                Some(emat) => (emat, Duration::ZERO),
                None => {
                    let t_int2e = Instant::now();
                    let mat2e = ints.mat2e(&self.int_form);
                    (mat2e, t_int2e.elapsed())
                }
            };
            (mat1e, mat2e, time_int1e, time_int2e)
        };
        let [smat, tmat, vmat] = mat1e;
        let nato = smat.shape()[0];
        assert!(nocc <= nato, "基函数数量({nato})不足以容纳 {nocc} 个占据轨道");
        let xmat = maths::matrix_nhalf(&smat);
        let hmat = &tmat + &vmat;

        // 核心哈密顿量作为初猜
        let (mut obt_engs, cmat_guess) = maths::eigh(&xmat.dot(&hmat).dot(&xmat));
        let mut cmat = xmat.dot(&cmat_guess);
        let mut dmat = self.calc_rdmat(nele, &cmat);
        let nuc_eng = self.nuc_eng();

        let mut ele_eng = f64::INFINITY;
        let mut delta_eng = f64::INFINITY;
        let mut converged = false;
        let mut niter = 0;
        let mut diis = Diis::new(8);
        let t_scf = Instant::now();
        for i in 1..=self.max_iter {
            niter = i;
            let (jmat, kmat) = self.calc_jk(&dmat, &emat);
            let fmat = &hmat + &jmat - &kmat * 0.5;
            // 能量始终用与当前密度一致的 Fock 矩阵计算
            let eng = rhf_ele_eng(&dmat, &hmat, &fmat);
            delta_eng = (eng - ele_eng).abs();
            println!("SCF 迭代 {i:>3}: 电子能量 = {eng:>18.10}  ΔE = {delta_eng:>9.2e}");
            ele_eng = eng;
            if delta_eng < self.tol {
                converged = true;
                break;
            }
            // DIIS 外推；历史不足时退回 0.5 阻尼
            let err = Diis::commutator(&fmat, &dmat, &smat);
            let extrap = diis.extrapolate(&fmat, &err);
            let fock = extrap.as_ref().unwrap_or(&fmat);
            let (engs, cmat_new) = maths::eigh(&xmat.dot(fock).dot(&xmat));
            let cmat_step = xmat.dot(&cmat_new);
            let dmat_step = self.calc_rdmat(nele, &cmat_step);
            obt_engs = engs;
            cmat = cmat_step;
            if extrap.is_some() {
                dmat = dmat_step;
            } else {
                dmat = &dmat * 0.5 + &dmat_step * 0.5;
            }
        }
        let time_scf = t_scf.elapsed();

        // 检查密度矩阵对应的电子数
        let dmat_fin = self.calc_rdmat(nele, &cmat);
        let nele_fin = (&dmat_fin * &smat).sum();
        if (nele_fin - nele as f64).abs() > 1.0e-4 {
            eprintln!("警告: 密度矩阵电子数 {nele_fin} 与预期 {nele} 不符");
        }

        // 构造与积分形式一致的基组和系数
        let basis = mole.basis().unify_form(&self.int_form);
        let atos = basis.atos();
        let ato_list = if self.int_form == "car" { &atos.cars } else { &atos.sphs };
        let ato_atms: Vec<usize> = ato_list.iter().map(|ato| ato.atm).collect();
        let ato_shls: Vec<usize> = ato_list.iter().map(|ato| ato.shl).collect();
        let ato_syms: Vec<String> = ato_list.iter().map(|ato| ato.sym.clone()).collect();
        let mut obt_occs = vec![0.0; nato];
        for occ in obt_occs.iter_mut().take(nocc) {
            *occ = 2.0;
        }
        let coefs = Coefs::new(ato_atms, ato_shls, ato_syms, obt_engs.to_vec(), obt_occs, cmat);

        RHFResult {
            ele_eng,
            nuc_eng,
            converged,
            niter,
            delta_eng,
            time_int1e,
            time_int2e,
            time_scf,
            coefs,
            basis,
            form: self.int_form.clone(),
        }
    }
}

/// 开壳层(自旋非限制)Hartree-Fock计算器
pub struct UHF<'m> {
    /// 分子结构引用
    pub mole: &'m Mole,
}

/// 根据开壳层 α/β 密度矩阵与电子积分矩阵计算电子能量
///
/// - `dmat`: α/β 密度矩阵，由各自的占据轨道系数构造
/// - `tmat` / `vmat`: 动能、电子-核吸引矩阵
/// - `fmat`: α/β 的 Fock 矩阵
///
/// 返回 `[E_α, E_β]`，其中 `E_σ = ½ Tr[D_σ(H + F_σ)]`、`H = T + V`。
pub fn uhf_ele_eng(
    dmat: &[Array2<f64>; 2],
    tmat: &Array2<f64>,
    vmat: &Array2<f64>,
    fmat: &[Array2<f64>; 2],
) -> [f64; 2] {
    let hmat = tmat + vmat;
    let eng_a = ((&hmat + &fmat[0]) * &dmat[0]).sum() * 0.5;
    let eng_b = ((&hmat + &fmat[1]) * &dmat[1]).sum() * 0.5;
    [eng_a, eng_b]
}

impl<'m> UHF<'m> {
    /// 创建新的UHF计算器
    pub fn new(mole: &'m Mole) -> Self {
        Self { mole }
    }
    /// 开壳层密度矩阵计算,对α和β电子分别计算
    pub fn calc_udmat(
        &self,
        nele: [usize; 2],        // alpha 和 beta 电子数
        cmat: &[Array2<f64>; 2], // 两个系数矩阵
    ) -> [Array2<f64>; 2] {
        let coef_a = cmat[0].slice(s![.., ..nele[0]]).to_owned();
        let coef_b = cmat[1].slice(s![.., ..nele[1]]).to_owned();
        let dmat_a = coef_a.dot(&coef_a.t());
        let dmat_b = coef_b.dot(&coef_b.t());
        [dmat_a, dmat_b]
    }

    fn calc_ufork(
        // 计算fork矩阵
        &self,
        nmat: usize,
        tmat: &Array2<f64>,      // 动能矩阵
        vmat: &Array2<f64>,      // 势能矩阵
        dmat: &[Array2<f64>; 2], // 密度矩阵
        emat: &Mat2e,            // 双电子积分矩阵
    ) -> [Array2<f64>; 2] {
        // 计算Fork矩阵

        let fmat_a = Array2::zeros((nmat, nmat));
        let fmat_b = Array2::zeros((nmat, nmat));
        let mut fmat = [fmat_a, fmat_b];
        let hmat = tmat + vmat;
        let dmat_sum = &dmat[0] + &dmat[1];
        let mut jsum;
        let mut csum;
        for s in 0..2usize {
            for i in 0..nmat {
                for j in 0..nmat {
                    jsum = 0.0; // 库伦矩阵元
                    csum = 0.0; // 交换矩阵元
                    for k in 0..nmat {
                        for l in 0..nmat {
                            jsum += emat.get(i, j, k, l) * dmat_sum[(k, l)];
                            csum += emat.get(i, k, j, l) * dmat[s][(k, l)];
                        }
                    }
                    fmat[s][(i, j)] = hmat[(i, j)] + jsum - csum;
                }
            }
        }

        return fmat;
    }

    /// 计算开壳层HF的自洽场
    ///
    /// 迭代求解Pople-Nesbet方程,返回α和β电子的能量、轨道能量和系数矩阵。
    pub fn calc_uhf(&self) -> ([f64; 2], [Array1<f64>; 2], [Array2<f64>; 2]) {
        let atos = &self.mole.basis().atos().sphs;
        let nato = atos.len();
        let cmat_a = Array2::zeros((nato, nato));
        let cmat_b = Array2::zeros((nato, nato));
        let mut cmat = [cmat_a, cmat_b];

        let dmat_a = Array2::zeros((nato, nato));
        let dmat_b = Array2::zeros((nato, nato));
        let mut dmat = [dmat_a, dmat_b];

        let rs_int = RsInt::new(self.mole);
        let [smat, tmat, vmat] = rs_int.mat1e([1, 1, 1], "car");
        let emat = rs_int.mat2e("car");

        let xmat = maths::matrix_nhalf(&smat);
        let mut eng = [0.0; 2];
        let mut obt_engs = [Array1::zeros(nato), Array1::zeros(nato)];
        let mut beta = 0.5; // 新旧密度矩阵的混合系数
        let mut old_eng = 100.0;
        let [nela, nelb] = self.mole.atoms.neles;
        log::info!("UHF计算");
        for i in 0..50 {
            let new_dmat = self.calc_udmat([nela, nelb], &cmat); // 计算密度矩阵
            for c in 0..2 {
                dmat[c] = &dmat[c] * (1.0 - beta) + &new_dmat[c] * beta;
            }
            let fmat = self.calc_ufork(nato, &tmat, &vmat, &dmat, &emat);
            for c in 0..2usize {
                let slice = &fmat[c];
                let fmat_ = xmat.dot(slice).dot(&xmat);
                let (_engs, _cmat) = maths::eigh(&fmat_);
                obt_engs[c] = _engs;
                cmat[c] = xmat.dot(&_cmat);
            }

            eng = uhf_ele_eng(&dmat, &tmat, &vmat, &fmat);

            if (eng[0] + eng[1] - old_eng).abs() < 1e-2 && beta != 1.0 {
                beta = 1.0;
            }
            if i > 10 {
                beta = 1.0;
            }

            log::info!(
                "epoch:{i:>2},电子能量:{:>10.4} + {:>10.4} = {:>14.8}",
                eng[0],
                eng[1],
                eng[0] + eng[1]
            );
            if (eng[0] + eng[1] - old_eng).abs() < 1e-6 {
                break;
            }
            old_eng = eng[0] + eng[1];
        }
        (eng, obt_engs, cmat)
    }
}
