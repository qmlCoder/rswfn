use crate::base::{Basis, Coefs, Mole};
use crate::maths;
use crate::props::gridprop;
use crate::props::gridprop::wfnfunc;
use crate::rsint::{Mat2e, RsInt};
use crate::rsxcs::{RsXcs, XcFunc};
use ndarray::{Array2, s};
use rayon::prelude::*;
use std::time::{Duration, Instant};

use super::diis::Diis;

/// B3LYP 中精确（Hartree-Fock）交换的混合系数
///
/// 等于 [`XcFunc::B3LYP.exact_x`]；保留为常量以兼容旧代码。
/// 其他泛函的系数应使用 [`XcFunc::exact_x`]。
pub const B3LYP_EXACT_X: f64 = 0.20;

/// 闭壳层 Kohn-Sham（DFT）计算结果
pub struct RKSResult {
    /// 电子能量（包含交换相关能）
    pub ele_eng: f64,
    /// 交换相关能
    pub xc_eng: f64,
    /// 核-核排斥能
    pub nuc_eng: f64,
    /// 交换相关泛函
    pub func: XcFunc,
    /// 是否自洽收敛
    pub converged: bool,
    /// 实际迭代次数
    pub niter: usize,
    /// 最后一次迭代的能量变化（Hartree）
    pub delta_eng: f64,
    /// 分子轨道系数
    pub coefs: Coefs,
    /// 与系数形式一致、可直接用于后续分析的基组
    pub basis: Basis,
    /// 实际使用的积分形式（"sph" 或 "car"）
    pub form: String,
    /// 数值积分格点数量
    pub ngrid: usize,
    /// 单电子积分耗时
    pub time_int1e: Duration,
    /// 双电子积分耗时
    pub time_int2e: Duration,
    /// SCF 迭代耗时
    pub time_scf: Duration,
}

impl RKSResult {
    /// 分子总能量（电子能量 + 核排斥能）
    pub fn total_eng(&self) -> f64 {
        self.ele_eng + self.nuc_eng
    }
}

/// 闭壳层 Kohn-Sham（DFT）计算器
///
/// 双电子部分用解析高斯积分构造库伦/交换矩阵，
/// 交换相关部分在原子中心格点上做数值积分。
/// 泛函由构造时的 [`XcFunc`] 指定；B3LYP 的局部相关采用 VWN3，
/// 与 Gaussian 一致。
pub struct RKS<'m> {
    /// 分子结构引用
    pub mole: &'m Mole,
    /// 交换相关泛函
    pub func: XcFunc,
    /// 积分形式（"sph" 或 "car"），由基组自动决定
    pub int_form: String,
    /// 径向格点数
    pub radial: usize,
    /// 球面格点数
    pub angular: usize,
    /// 最大SCF迭代次数
    pub max_iter: usize,
    /// 收敛阈值（相邻两次迭代的能量变化，Hartree）
    pub tol: f64,
    /// 所有格点坐标
    pub grids: Vec<[f64; 3]>,
    /// 所有格点权重
    pub weits: Vec<f64>,
    /// 各格点处所有原子轨道的取值
    pub wfn0s: Vec<Vec<f64>>,
    /// 各格点处所有原子轨道的一阶导数
    pub wfn1s: Vec<Vec<[f64; 3]>>,
}

/// 计算某个格点处所有原子轨道（指定形式）的取值与一阶导数
///
/// `wfnfunc::ato_wfn` 固定使用笛卡尔形式，而 KS 计算需要与积分形式
/// （`sph` 或 `car`）一致，否则系数矩阵与格点上的轨道值会错位。
fn ato_wfn_form(mole: &Mole, grid: &[f64; 3], form: &str) -> (Vec<f64>, Vec<[f64; 3]>) {
    let atos = mole.basis().atos();
    let atos = if form == "car" { &atos.cars } else { &atos.sphs };
    let atoms = mole.atoms();
    let nato = atos.len();
    let mut wfn0 = vec![0.0; nato];
    let mut wfn1 = vec![[0.0; 3]; nato];
    for i in 0..nato {
        let atm = atos[i].atm;
        let [ax, ay, az] = atoms.get(atm).xyz;
        let xyz = [grid[0] - ax, grid[1] - ay, grid[2] - az];
        let coes = atos[i].coes_norm();
        let alps = &atos[i].alps;
        let lmn = atos[i].lmn();
        let (val0, val1, _val2) = wfnfunc::cgf(&xyz, lmn, alps, coes, 1);
        wfn0[i] = val0;
        wfn1[i] = val1;
    }
    (wfn0, wfn1)
}

/// 并行遍历格点时的最小分块长度
///
/// `fold` 会为每个分块分配一份 `nato²` 累加器；分块太小会制造大量
/// 临时分配（每轮数百 MB），所以按线程数限制分块数量。
fn grid_min_len(ngrid: usize) -> usize {
    let chunks = rayon::current_num_threads() * 8;
    (ngrid / chunks).max(64)
}

/// 根据闭壳层占据轨道系数与电子积分矩阵计算 KS 电子能量
///
/// - `dmat`: 由占据轨道系数构造的密度矩阵，`D = 2 C_occ C_occᵀ`
/// - `hmat`: 核心哈密顿量（动能 + 电子-核吸引）
/// - `jmat` / `kmat`: 由双电子积分构造的库伦/交换矩阵
/// - `exact_x`: 泛函中精确交换的混合系数 `a_x`
/// - `xc_eng`: 格点积分得到的交换相关能
///
/// 返回 `E = Tr(DH) + ½Tr(DJ) − ¼ a_x Tr(DK) + E_xc`。
pub fn rks_ele_eng(dmat: &Array2<f64>, hmat: &Array2<f64>, jmat: &Array2<f64>, kmat: &Array2<f64>, exact_x: f64, xc_eng: f64) -> f64 {
    (dmat * hmat).sum() + (dmat * jmat).sum() * 0.5 - (dmat * kmat).sum() * exact_x * 0.25 + xc_eng
}

impl<'m> RKS<'m> {
    /// 创建新的 RKS 计算器（默认格点：100 径向 × 74 球面）
    ///
    /// `func` 指定交换相关泛函（如 [`XcFunc::B3LYP`]）。
    /// 初始化分子格点、权重，以及所有原子轨道在格点上的取值和一阶导数。
    pub fn new(mole: &'m Mole, func: XcFunc) -> Self {
        Self::with_grid(mole, 100, 74, func)
    }

    /// 按泛函名创建（如 `"b3lyp"`、`"pbe0"`）；未知名称返回 `None`
    pub fn from_name(mole: &'m Mole, name: &str) -> Option<Self> {
        XcFunc::from_name(name).map(|func| Self::new(mole, func))
    }

    /// 用指定的径向/球面格点数与泛函创建 RKS 计算器
    pub fn with_grid(mole: &'m Mole, radial: usize, angular: usize, func: XcFunc) -> Self {
        let int_form = mole.basis().int_form().to_string();
        let grid_caler = gridprop::dftgrid::Calculator::new(mole);
        let (grids, weits) = grid_caler.mol_grids(radial, angular);
        let wfns: Vec<(Vec<f64>, Vec<[f64; 3]>)> = grids.par_iter().map(|grid| ato_wfn_form(mole, grid, &int_form)).collect();
        let (wfn0s, wfn1s): (Vec<_>, Vec<_>) = wfns.into_iter().unzip();
        Self {
            mole,
            func,
            int_form,
            radial,
            angular,
            max_iter: 60,
            tol: 1e-8,
            grids,
            weits,
            wfn0s,
            wfn1s,
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

    /// 计算所有格点处的电子密度、密度梯度与 |∇ρ|²
    pub fn calc_rho(&self, dmat: &Array2<f64>) -> (Vec<f64>, Vec<[f64; 3]>, Vec<f64>) {
        let nato = dmat.nrows();
        // 摊平密度矩阵，避免热循环里反复做 ndarray 的边界检查与步长计算
        let dvals: Vec<f64> = dmat.iter().copied().collect();
        let vals: Vec<(f64, [f64; 3])> = self
            .grids
            .par_iter()
            .enumerate()
            .map(|(p, _grid)| {
                let wfn0 = &self.wfn0s[p];
                let wfn1 = &self.wfn1s[p];
                let mut rho0 = 0.0;
                let mut rho1 = [0.0; 3];
                for i in 0..nato {
                    let fi = wfn0[i];
                    let gi = wfn1[i];
                    if fi == 0.0 && gi == [0.0; 3] {
                        continue;
                    }
                    let irow = i * nato;
                    for j in 0..nato {
                        let d = dvals[irow + j];
                        if d == 0.0 {
                            continue;
                        }
                        let fj = wfn0[j];
                        let gj = wfn1[j];
                        rho0 += d * fi * fj;
                        rho1[0] += d * (fi * gj[0] + fj * gi[0]);
                        rho1[1] += d * (fi * gj[1] + fj * gi[1]);
                        rho1[2] += d * (fi * gj[2] + fj * gi[2]);
                    }
                }
                (rho0, rho1)
            })
            .collect();

        let mut rho0s = Vec::with_capacity(vals.len());
        let mut rho1s = Vec::with_capacity(vals.len());
        let mut sigmas = Vec::with_capacity(vals.len());
        for (rho0, rho1) in vals {
            let [dx, dy, dz] = rho1;
            sigmas.push(dx * dx + dy * dy + dz * dz);
            rho0s.push(rho0);
            rho1s.push(rho1);
        }
        (rho0s, rho1s, sigmas)
    }

    // 计算交换相关矩阵
    //
    // V_ij = Σ_p w_p [ φ_i φ_j v_ρ + 2 v_σ ∇ρ·(φ_i ∇φ_j + φ_j ∇φ_i) ]
    //
    // 按格点并行，每个任务累加一个局部矩阵，最后求和。
    fn calc_xcmat(&self, vrhos: &[f64], vsigmas: &[f64], rho1s: &[[f64; 3]]) -> Array2<f64> {
        let nato = self.wfn0s[0].len();
        // 用扁平缓冲区累加，避免在热循环里做 ndarray 下标运算
        let zero = || vec![0.0f64; nato * nato];
        let mut acc = (0..self.grids.len())
            .into_par_iter()
            .with_min_len(grid_min_len(self.grids.len()))
            .fold(zero, |mut acc, p| {
                let weit = self.weits[p];
                if weit == 0.0 {
                    return acc;
                }
                let vrho = vrhos[p];
                let vsigma = vsigmas[p];
                let [gx, gy, gz] = rho1s[p];
                let wfn0 = &self.wfn0s[p];
                let wfn1 = &self.wfn1s[p];
                for i in 0..nato {
                    let fi = wfn0[i];
                    let gi = wfn1[i];
                    if fi == 0.0 && gi == [0.0; 3] {
                        continue;
                    }
                    let irow = i * nato;
                    for j in 0..=i {
                        let fj = wfn0[j];
                        let gj = wfn1[j];
                        let val_lda = fi * fj * vrho;
                        let term = (fi * gj[0] + fj * gi[0]) * gx + (fi * gj[1] + fj * gi[1]) * gy + (fi * gj[2] + fj * gi[2]) * gz;
                        let val_gga = 2.0 * vsigma * term;
                        acc[irow + j] += (val_lda + val_gga) * weit;
                    }
                }
                acc
            })
            .reduce(zero, |mut a, b| {
                for (x, y) in a.iter_mut().zip(b.iter()) {
                    *x += y;
                }
                a
            });
        // 只累加了下三角，这里补全对称部分
        for i in 0..nato {
            for j in 0..i {
                acc[j * nato + i] = acc[i * nato + j];
            }
        }
        Array2::from_shape_vec((nato, nato), acc).expect("xcmat 维度不匹配")
    }

    /// 给定密度矩阵，返回格点上积分得到的交换相关能
    ///
    /// 只用于诊断/验证（自洽场内部会一并计算密度、势与能量）。
    pub fn xc_eng(&self, dmat: &Array2<f64>) -> f64 {
        let (rho0s, _rho1s, sigmas) = self.calc_rho(dmat);
        let (zk, _vrho, _vsigma) = RsXcs::with_func(self.func).run(&rho0s, &sigmas);
        (0..rho0s.len()).map(|p| zk[p] * rho0s[p] * self.weits[p]).sum()
    }

    /// 给定密度矩阵，返回格点上积分得到的电子数
    pub fn grid_ele(&self, dmat: &Array2<f64>) -> f64 {
        let (rho0s, _rho1s, _sigmas) = self.calc_rho(dmat);
        (0..rho0s.len()).map(|p| rho0s[p] * self.weits[p]).sum()
    }

    /// 执行闭壳层 KS-DFT 自洽场计算（泛函由构造时指定）
    pub fn run(&self) -> RKSResult {
        let mole = self.mole;
        let [nela, nelb] = mole.atoms.neles;
        assert_eq!(nela, nelb, "RKS 仅支持闭壳层体系：α电子数({nela}) ≠ β电子数({nelb})");
        let nele = nela + nelb;
        let nocc = nele / 2;
        assert!(nocc > 0, "分子没有电子，无法进行 DFT 计算");

        // 电子积分
        let int_caler = RsInt::new(mole);
        let t_int1e = Instant::now();
        let [smat, tmat, vmat] = int_caler.mat1e([1, 1, 1], &self.int_form);
        let time_int1e = t_int1e.elapsed();
        let t_int2e = Instant::now();
        let emat = int_caler.mat2e(&self.int_form);
        let time_int2e = t_int2e.elapsed();
        let hmat = &tmat + &vmat;
        let nato = smat.shape()[0];
        assert!(nocc <= nato, "基函数数量({nato})不足以容纳 {nocc} 个占据轨道");
        let xmat = maths::matrix_nhalf(&smat);

        let rsfun = RsXcs::with_func(self.func);
        let exact_x = self.func.exact_x();
        let nuc_eng = self.nuc_eng();

        // 核心哈密顿量作为初猜
        let (_engs0, cmat_guess) = maths::eigh(&xmat.dot(&hmat).dot(&xmat));
        let mut cmat = xmat.dot(&cmat_guess);
        let mut dmat = self.calc_rdmat(nele, &cmat);
        let mut obt_engs = Vec::new();

        let mut ele_eng = f64::INFINITY;
        let mut xc_eng = 0.0;
        let mut delta_eng = f64::INFINITY;
        let mut converged = false;
        let mut niter = 0;
        let mut diis = Diis::new(8);
        let t_scf = Instant::now();
        for i in 1..=self.max_iter {
            niter = i;
            // 格点上的密度与梯度
            let (rho0s, rho1s, sigmas) = self.calc_rho(&dmat);
            let (zk, vrho, vsigma) = rsfun.run(&rho0s, &sigmas);
            let xcmat = self.calc_xcmat(&vrho, &vsigma, &rho1s);
            xc_eng = (0..rho0s.len()).map(|p| zk[p] * rho0s[p] * self.weits[p]).sum();

            // 库伦/交换矩阵与 Fock 矩阵
            // 精确交换按泛函系数 a_x 加入：F += -a_x * K/2
            let (jmat, kmat) = crate::rsint::jk_mat(&dmat, &emat);
            let fmat = &hmat + &jmat - &kmat * (exact_x * 0.5) + &xcmat;

            // 能量：Tr(DH) + 1/2 Tr(DJ) - a_x/4 Tr(DK) + E_xc（与当前密度一致）
            let eng = rks_ele_eng(&dmat, &hmat, &jmat, &kmat, exact_x, xc_eng);
            delta_eng = (eng - ele_eng).abs();
            log::info!("SCF 迭代 {i:>3}: 电子能量 = {eng:>18.10}  ΔE = {delta_eng:>9.2e}  E_xc = {xc_eng:>16.8}");
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
            obt_engs = engs.to_vec();
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
        let coefs = Coefs::new(ato_atms, ato_shls, ato_syms, obt_engs, obt_occs, cmat);

        RKSResult {
            ele_eng,
            xc_eng,
            nuc_eng,
            func: self.func,
            converged,
            niter,
            delta_eng,
            coefs,
            basis,
            form: self.int_form.clone(),
            ngrid: self.grids.len(),
            time_int1e,
            time_int2e,
            time_scf,
        }
    }
}

/// 自旋分辨的格点密度信息（α/β 密度、梯度及梯度内积）
pub struct SpinRho {
    /// α 电子密度在各格点的值
    pub rho_a: Vec<f64>,
    /// β 电子密度在各格点的值
    pub rho_b: Vec<f64>,
    /// α 电子密度梯度
    pub grad_a: Vec<[f64; 3]>,
    /// β 电子密度梯度
    pub grad_b: Vec<[f64; 3]>,
    /// σ_αα = |∇ρ_α|²
    pub saa: Vec<f64>,
    /// σ_ββ = |∇ρ_β|²
    pub sbb: Vec<f64>,
    /// σ_αβ = ∇ρ_α·∇ρ_β
    pub sab: Vec<f64>,
}

/// 由自旋轨道系数构造分自旋密度矩阵 `D_σ = C_occ · C_occ^T`（每个自旋轨道占据数 1）
fn calc_spin_dmat(nocc: usize, cmat: &Array2<f64>) -> Array2<f64> {
    let coef = cmat.slice(s![.., ..nocc]).to_owned();
    coef.dot(&coef.t())
}

/// 自旋极化（开壳层）Kohn-Sham（DFT）计算结果
pub struct UKSResult {
    /// 电子能量（包含交换相关能）
    pub ele_eng: f64,
    /// 交换相关能
    pub xc_eng: f64,
    /// 核-核排斥能
    pub nuc_eng: f64,
    /// 交换相关泛函
    pub func: XcFunc,
    /// 是否自洽收敛
    pub converged: bool,
    /// 实际迭代次数
    pub niter: usize,
    /// 最后一次迭代的能量变化（Hartree）
    pub delta_eng: f64,
    /// 分子轨道系数（α 块在前、β 块在后）
    pub coefs: Coefs,
    /// 与系数形式一致、可直接用于后续分析的基组
    pub basis: Basis,
    /// 实际使用的积分形式（"sph" 或 "car"）
    pub form: String,
    /// 数值积分格点数量
    pub ngrid: usize,
    /// α 电子数
    pub nela: usize,
    /// β 电子数
    pub nelb: usize,
    /// ⟨S²⟩ 期望值
    pub spin_sq: f64,
    /// 单电子积分耗时
    pub time_int1e: Duration,
    /// 双电子积分耗时
    pub time_int2e: Duration,
    /// SCF 迭代耗时
    pub time_scf: Duration,
}

impl UKSResult {
    /// 分子总能量（电子能量 + 核排斥能）
    pub fn total_eng(&self) -> f64 {
        self.ele_eng + self.nuc_eng
    }

    /// 自旋多重度 `2S + 1 = n_α − n_β + 1`
    pub fn multip(&self) -> usize {
        self.nela - self.nelb + 1
    }
}

/// 自旋极化 Kohn-Sham（DFT）计算器（开壳层）
///
/// 与 [`RKS`] 共用同一套原子中心数值积分格点，α、β 两套密度各自自洽。
/// 库伦项用总密度，交换项按自旋分别处理，交换相关泛函走
/// [`RsXcs::run_spin`]（自旋极化版本，B3LYP 的局部相关仍为 VWN3）。
pub struct UKS<'m> {
    /// 复用闭壳层计算器的格点与积分设置
    pub base: RKS<'m>,
    /// α 电子数
    pub nela: usize,
    /// β 电子数
    pub nelb: usize,
}

/// 根据开壳层 α/β 轨道密度与电子积分矩阵计算 KS 电子能量
///
/// - `dmat`: 总密度矩阵，`D = D_α + D_β`
/// - `dmat_a` / `dmat_b`: α/β 分自旋密度矩阵，由各自的占据轨道系数构造
/// - `hmat`: 核心哈密顿量（动能 + 电子-核吸引）
/// - `jmat`: 由总密度构造的库伦矩阵
/// - `kmat_a` / `kmat_b`: α/β 分自旋交换矩阵
/// - `exact_x`: 精确交换混合系数 `a_x`
/// - `xc_eng`: 格点积分得到的交换相关能
///
/// 返回 `E = Tr(DH) + ½Tr(DJ) − ½ a_x [Tr(D_α K_α) + Tr(D_β K_β)] + E_xc`。
pub fn uks_ele_eng(dmat: &Array2<f64>, hmat: &Array2<f64>, jmat: &Array2<f64>, dmat_a: &Array2<f64>, dmat_b: &Array2<f64>, kmat_a: &Array2<f64>, kmat_b: &Array2<f64>, exact_x: f64, xc_eng: f64) -> f64 {
    (dmat * hmat).sum() + (dmat * jmat).sum() * 0.5 - ((dmat_a * kmat_a).sum() + (dmat_b * kmat_b).sum()) * exact_x * 0.5 + xc_eng
}

impl<'m> UKS<'m> {
    /// 创建新的 UKS 计算器（默认格点：100 径向 × 74 球面）
    ///
    /// `func` 指定交换相关泛函（如 [`XcFunc::B3LYP`]）。
    pub fn new(mole: &'m Mole, func: XcFunc) -> Self {
        Self::with_grid(mole, 100, 74, func)
    }

    /// 按泛函名创建（如 `"b3lyp"`、`"pbe0"`）；未知名称返回 `None`
    pub fn from_name(mole: &'m Mole, name: &str) -> Option<Self> {
        XcFunc::from_name(name).map(|func| Self::new(mole, func))
    }

    /// 用指定的径向/球面格点数与泛函创建 UKS 计算器
    pub fn with_grid(mole: &'m Mole, radial: usize, angular: usize, func: XcFunc) -> Self {
        let [nela, nelb] = mole.atoms.neles;
        assert!(nela >= nelb, "UKS 要求 α 电子数 ≥ β 电子数（{nela} < {nelb}）");
        Self {
            base: RKS::with_grid(mole, radial, angular, func),
            nela,
            nelb,
        }
    }

    /// 计算核-核排斥能
    pub fn nuc_eng(&self) -> f64 {
        self.base.nuc_eng()
    }

    /// 计算所有格点处的 α/β 密度、梯度与梯度内积
    ///
    /// α、β 在同一次格点遍历中完成，避免把 `φ`、`∇φ` 读两遍。
    pub fn calc_rho_spin(&self, dmat_a: &Array2<f64>, dmat_b: &Array2<f64>) -> SpinRho {
        let base = &self.base;
        let nato = dmat_a.nrows();
        // 摊平密度矩阵，避免热循环里反复做 ndarray 的边界检查与步长计算
        let da: Vec<f64> = dmat_a.iter().copied().collect();
        let db: Vec<f64> = dmat_b.iter().copied().collect();
        let vals: Vec<(f64, [f64; 3], f64, [f64; 3])> = base
            .grids
            .par_iter()
            .enumerate()
            .map(|(p, _grid)| {
                let wfn0 = &base.wfn0s[p];
                let wfn1 = &base.wfn1s[p];
                let mut rho_a = 0.0;
                let mut rho_b = 0.0;
                let mut grad_a = [0.0; 3];
                let mut grad_b = [0.0; 3];
                for i in 0..nato {
                    let fi = wfn0[i];
                    let gi = wfn1[i];
                    if fi == 0.0 && gi == [0.0; 3] {
                        continue;
                    }
                    let irow = i * nato;
                    for j in 0..nato {
                        let d_a = da[irow + j];
                        let d_b = db[irow + j];
                        if d_a == 0.0 && d_b == 0.0 {
                            continue;
                        }
                        let fj = wfn0[j];
                        let gj = wfn1[j];
                        let fij = fi * fj;
                        let gx = fi * gj[0] + fj * gi[0];
                        let gy = fi * gj[1] + fj * gi[1];
                        let gz = fi * gj[2] + fj * gi[2];
                        if d_a != 0.0 {
                            rho_a += d_a * fij;
                            grad_a[0] += d_a * gx;
                            grad_a[1] += d_a * gy;
                            grad_a[2] += d_a * gz;
                        }
                        if d_b != 0.0 {
                            rho_b += d_b * fij;
                            grad_b[0] += d_b * gx;
                            grad_b[1] += d_b * gy;
                            grad_b[2] += d_b * gz;
                        }
                    }
                }
                (rho_a, grad_a, rho_b, grad_b)
            })
            .collect();

        let n = vals.len();
        let mut rho_a = Vec::with_capacity(n);
        let mut rho_b = Vec::with_capacity(n);
        let mut grad_a = Vec::with_capacity(n);
        let mut grad_b = Vec::with_capacity(n);
        let mut saa = Vec::with_capacity(n);
        let mut sbb = Vec::with_capacity(n);
        let mut sab = Vec::with_capacity(n);
        for (ra, ga, rb, gb) in vals {
            saa.push(ga[0] * ga[0] + ga[1] * ga[1] + ga[2] * ga[2]);
            sbb.push(gb[0] * gb[0] + gb[1] * gb[1] + gb[2] * gb[2]);
            sab.push(ga[0] * gb[0] + ga[1] * gb[1] + ga[2] * gb[2]);
            rho_a.push(ra);
            rho_b.push(rb);
            grad_a.push(ga);
            grad_b.push(gb);
        }
        SpinRho {
            rho_a,
            rho_b,
            grad_a,
            grad_b,
            saa,
            sbb,
            sab,
        }
    }

    // 一次性构造 α、β 的交换相关矩阵
    //
    // V^σ_ij = Σ_p w_p [ φ_i φ_j v_ρσ + G_ij·(2 v_σσ ∇ρ_σ + v_αβ ∇ρ_σ' ) ]
    // 其中 G_ij = φ_i ∇φ_j + φ_j ∇φ_i。按格点并行，每个任务累加两个局部矩阵。
    fn calc_xcmat_spin(&self, rho: &SpinRho, vr_a: &[f64], vr_b: &[f64], vs_aa: &[f64], vs_bb: &[f64], vs_ab: &[f64]) -> (Array2<f64>, Array2<f64>) {
        let nato = self.base.wfn0s[0].len();
        let zero = || (vec![0.0f64; nato * nato], vec![0.0f64; nato * nato]);
        let mut acc = (0..self.base.grids.len())
            .into_par_iter()
            .with_min_len(grid_min_len(self.base.grids.len()))
            .fold(zero, |(mut xa, mut xb), p| {
                let weit = self.base.weits[p];
                if weit == 0.0 {
                    return (xa, xb);
                }
                let wfn0 = &self.base.wfn0s[p];
                let wfn1 = &self.base.wfn1s[p];
                let ga = rho.grad_a[p];
                let gb = rho.grad_b[p];
                let (vra, vrb) = (vr_a[p], vr_b[p]);
                let (vaa, vbb, vab) = (vs_aa[p], vs_bb[p], vs_ab[p]);
                for i in 0..nato {
                    let fi = wfn0[i];
                    let gi = wfn1[i];
                    if fi == 0.0 && gi == [0.0; 3] {
                        continue;
                    }
                    let irow = i * nato;
                    for j in 0..=i {
                        let fj = wfn0[j];
                        let gj = wfn1[j];
                        let gx = fi * gj[0] + fj * gi[0];
                        let gy = fi * gj[1] + fj * gi[1];
                        let gz = fi * gj[2] + fj * gi[2];
                        let term_a = gx * ga[0] + gy * ga[1] + gz * ga[2];
                        let term_b = gx * gb[0] + gy * gb[1] + gz * gb[2];
                        let fij = fi * fj;
                        let val_a = fij * vra + 2.0 * vaa * term_a + vab * term_b;
                        let val_b = fij * vrb + 2.0 * vbb * term_b + vab * term_a;
                        xa[irow + j] += val_a * weit;
                        xb[irow + j] += val_b * weit;
                    }
                }
                (xa, xb)
            })
            .reduce(zero, |(mut xa, mut xb), (ya, yb)| {
                for (x, y) in xa.iter_mut().zip(ya.iter()) {
                    *x += y;
                }
                for (x, y) in xb.iter_mut().zip(yb.iter()) {
                    *x += y;
                }
                (xa, xb)
            });
        // 只累加了下三角，这里补全对称部分
        for i in 0..nato {
            for j in 0..i {
                acc.0[j * nato + i] = acc.0[i * nato + j];
                acc.1[j * nato + i] = acc.1[i * nato + j];
            }
        }
        (
            Array2::from_shape_vec((nato, nato), acc.0).expect("xcmat(α) 维度不匹配"),
            Array2::from_shape_vec((nato, nato), acc.1).expect("xcmat(β) 维度不匹配"),
        )
    }

    /// 计算 ⟨S²⟩ = s_z(s_z+1) + n_β − Tr(D_α S D_β S)
    fn spin_sq(&self, dmat_a: &Array2<f64>, dmat_b: &Array2<f64>, smat: &Array2<f64>) -> f64 {
        let n = dmat_a.nrows();
        let sz = (self.nela as f64 - self.nelb as f64) / 2.0;
        let m = dmat_a.dot(smat).dot(dmat_b).dot(smat);
        let mut tr = 0.0;
        for i in 0..n {
            tr += m[(i, i)];
        }
        sz * (sz + 1.0) + self.nelb as f64 - tr
    }

    /// 执行自旋极化 KS-DFT 自洽场计算（泛函由构造时指定）
    pub fn run(&self) -> UKSResult {
        let mole = self.base.mole;
        let (nela, nelb) = (self.nela, self.nelb);
        assert!(nela > 0, "分子没有 α 电子，无法进行 UKS 计算");

        // 电子积分
        let int_caler = RsInt::new(mole);
        let t_int1e = Instant::now();
        let [smat, tmat, vmat] = int_caler.mat1e([1, 1, 1], &self.base.int_form);
        let time_int1e = t_int1e.elapsed();
        let t_int2e = Instant::now();
        let emat = int_caler.mat2e(&self.base.int_form);
        let time_int2e = t_int2e.elapsed();
        let hmat = &tmat + &vmat;
        let nato = smat.shape()[0];
        assert!(nela <= nato && nelb <= nato, "基函数数量({nato})不足以容纳占据轨道");
        let xmat = maths::matrix_nhalf(&smat);

        let rsfun = RsXcs::with_func(self.base.func);
        let exact_x = self.base.func.exact_x();
        let nuc_eng = self.base.nuc_eng();

        // 核心哈密顿量作为初猜（α、β 初始相同）
        let (engs0, cmat_guess) = maths::eigh(&xmat.dot(&hmat).dot(&xmat));
        let cmat0 = xmat.dot(&cmat_guess);
        let mut cmat_a = cmat0.clone();
        let mut cmat_b = cmat0;
        let mut dmat_a = calc_spin_dmat(nela, &cmat_a);
        let mut dmat_b = calc_spin_dmat(nelb, &cmat_b);
        let mut engs_a = engs0.to_vec();
        let mut engs_b = engs0.to_vec();

        let mut ele_eng = f64::INFINITY;
        let mut xc_eng = 0.0;
        let mut delta_eng = f64::INFINITY;
        let mut converged = false;
        let mut niter = 0;
        let mut diis_a = Diis::new(8);
        let mut diis_b = Diis::new(8);
        let t_scf = Instant::now();

        for i in 1..=self.base.max_iter {
            niter = i;
            let dmat = &dmat_a + &dmat_b;
            let rho = self.calc_rho_spin(&dmat_a, &dmat_b);
            let (zk, vr_a, vr_b, vs_aa, vs_bb, vs_ab) = rsfun.run_spin(&rho.rho_a, &rho.rho_b, &rho.saa, &rho.sbb, &rho.sab);
            let (xc_a, xc_b) = self.calc_xcmat_spin(&rho, &vr_a, &vr_b, &vs_aa, &vs_bb, &vs_ab);
            xc_eng = (0..rho.rho_a.len()).map(|p| zk[p] * (rho.rho_a[p] + rho.rho_b[p]) * self.base.weits[p]).sum();

            // 库伦项用总密度，交换项按自旋分别处理；精确交换按 a_x 加入
            let (jmat, kmat_a, kmat_b) = crate::rsint::jk_mat_spin(&dmat, &dmat_a, &dmat_b, &emat);
            let fmat_a = &hmat + &jmat - &kmat_a * exact_x + &xc_a;
            let fmat_b = &hmat + &jmat - &kmat_b * exact_x + &xc_b;

            // 能量：Tr(D_tot H) + ½Tr(D_tot J) − ½ a_x [Tr(D_a K_a) + Tr(D_b K_b)] + E_xc
            let eng = uks_ele_eng(&dmat, &hmat, &jmat, &dmat_a, &dmat_b, &kmat_a, &kmat_b, exact_x, xc_eng);
            delta_eng = (eng - ele_eng).abs();
            log::info!("SCF 迭代 {i:>3}: 电子能量 = {eng:>18.10}  ΔE = {delta_eng:>9.2e}  E_xc = {xc_eng:>16.8}");
            ele_eng = eng;
            if delta_eng < self.base.tol {
                converged = true;
                break;
            }

            // α、β 各用一套 DIIS；历史不足时退回 0.5 阻尼
            let err_a = Diis::commutator(&fmat_a, &dmat_a, &smat);
            let err_b = Diis::commutator(&fmat_b, &dmat_b, &smat);
            let ext_a = diis_a.extrapolate(&fmat_a, &err_a);
            let ext_b = diis_b.extrapolate(&fmat_b, &err_b);
            let fock_a = ext_a.as_ref().unwrap_or(&fmat_a);
            let fock_b = ext_b.as_ref().unwrap_or(&fmat_b);

            let (e_a, c_a) = maths::eigh(&xmat.dot(fock_a).dot(&xmat));
            let (e_b, c_b) = maths::eigh(&xmat.dot(fock_b).dot(&xmat));
            let cmat_a_new = xmat.dot(&c_a);
            let cmat_b_new = xmat.dot(&c_b);
            let dmat_a_new = calc_spin_dmat(nela, &cmat_a_new);
            let dmat_b_new = calc_spin_dmat(nelb, &cmat_b_new);
            engs_a = e_a.to_vec();
            engs_b = e_b.to_vec();
            cmat_a = cmat_a_new;
            cmat_b = cmat_b_new;
            if ext_a.is_some() {
                dmat_a = dmat_a_new;
            } else {
                dmat_a = &dmat_a * 0.5 + &dmat_a_new * 0.5;
            }
            if ext_b.is_some() {
                dmat_b = dmat_b_new;
            } else {
                dmat_b = &dmat_b * 0.5 + &dmat_b_new * 0.5;
            }
        }
        let time_scf = t_scf.elapsed();

        // 检查电子数
        let dmat_a_fin = calc_spin_dmat(nela, &cmat_a);
        let dmat_b_fin = calc_spin_dmat(nelb, &cmat_b);
        let na_fin = (&dmat_a_fin * &smat).sum();
        let nb_fin = (&dmat_b_fin * &smat).sum();
        if (na_fin - nela as f64).abs() > 1.0e-4 || (nb_fin - nelb as f64).abs() > 1.0e-4 {
            eprintln!("警告: α/β 电子数 {na_fin}/{nb_fin} 与预期 {nela}/{nelb} 不符");
        }
        let spin_sq = self.spin_sq(&dmat_a_fin, &dmat_b_fin, &smat);

        // 构造与积分形式一致的基组与系数（α 块在前、β 块在后）
        let basis = mole.basis().unify_form(&self.base.int_form);
        let atos = basis.atos();
        let ato_list = if self.base.int_form == "car" { &atos.cars } else { &atos.sphs };
        let ato_atms: Vec<usize> = ato_list.iter().map(|ato| ato.atm).collect();
        let ato_shls: Vec<usize> = ato_list.iter().map(|ato| ato.shl).collect();
        let ato_syms: Vec<String> = ato_list.iter().map(|ato| ato.sym.clone()).collect();

        let mut cmat = Array2::zeros((nato, 2 * nato));
        cmat.slice_mut(s![.., ..nato]).assign(&cmat_a);
        cmat.slice_mut(s![.., nato..]).assign(&cmat_b);
        let mut obt_engs = engs_a;
        obt_engs.extend(engs_b);
        let mut obt_occs = vec![0.0; 2 * nato];
        for occ in obt_occs.iter_mut().take(nela) {
            *occ = 1.0;
        }
        for occ in obt_occs.iter_mut().skip(nato).take(nelb) {
            *occ = 1.0;
        }
        let coefs = Coefs::new(ato_atms, ato_shls, ato_syms, obt_engs, obt_occs, cmat);

        UKSResult {
            ele_eng,
            xc_eng,
            nuc_eng,
            func: self.base.func,
            converged,
            niter,
            delta_eng,
            coefs,
            basis,
            form: self.base.int_form.clone(),
            ngrid: self.base.grids.len(),
            nela,
            nelb,
            spin_sq,
            time_int1e,
            time_int2e,
            time_scf,
        }
    }
}
