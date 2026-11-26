use crate::base::{Ato, Basis, Coefs, Mole};
use crate::maths;
use crate::props::gridprop;
use crate::props::gridprop::wfnfunc;
use crate::rsint::RsInt;
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
    /// SCF 中 J/K 构造的累计耗时
    pub time_jk: Duration,
    /// SCF 中格点相关（密度、XC 矩阵）的累计耗时
    pub time_grid: Duration,
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
    /// 格点总数
    pub ngrid: usize,
    /// 原子轨道数（与积分形式一致）
    pub nato: usize,
    /// 各格点处所有原子轨道的取值，按 `[p * nato + i]` 展开
    pub chi0: Vec<f64>,
    /// 各格点处所有原子轨道的一阶导数，按 `[p * 3 * nato + r * nato + i]` 展开
    pub chi1: Vec<f64>,
    /// Schwarz 筛选阈值（双电子积分）
    pub screen: f64,
    /// 双电子积分块允许占用的最大内存（字节），超过则退化为直接 SCF
    pub max_mem: usize,
    /// 显著原子轨道的距离截断（Bohr）；0 表示关闭（每格点用全部 AO）
    ///
    /// 打开后，每个格点块只保留与块内格点距离小于该截断的原子轨道，
    /// 大幅减少大分子格点积分的浮点量与带宽（见 `calc_rho_ao`）。
    pub ao_cut: f64,
    /// 原子间距矩阵，按 `[i * natm + j]` 存放
    atm_dist: Vec<f64>,
    /// 每个原子轨道所属的原子索引
    ao_atm: Vec<usize>,
    /// 每个原子轨道的有效半径（Bohr），由基函数径向部分衰减到 1e-7 决定
    ao_r: Vec<f64>,
    /// 按原子分组的原子轨道索引（升序）
    ao_by_atm: Vec<Vec<usize>>,
    /// 每个格点的所属原子
    grid_atm: Vec<usize>,
    /// 每个格点到所属原子中心的距离（Bohr）
    grid_r: Vec<f64>,
    /// 格点分块（与原子边界对齐），供显著 AO 路径使用
    blocks: Vec<GridBlock>,
}

/// 一个格点块：块内格点同属一个原子，半径落在 `[rmin, rmax]`
#[derive(Clone, Copy, Debug)]
pub struct GridBlock {
    /// 起始格点下标
    pub p0: usize,
    /// 结束格点下标（不含）
    pub p1: usize,
    /// 所属原子
    pub atm: usize,
    /// 块内最小径向半径
    pub rmin: f64,
    /// 块内最大径向半径
    pub rmax: f64,
}

/// 计算某个格点处所有原子轨道（指定形式）的取值与一阶导数
///
/// `wfnfunc::ato_wfn` 固定使用笛卡尔形式，而 KS 计算需要与积分形式
/// （`sph` 或 `car`）一致，否则系数矩阵与格点上的轨道值会错位。
fn ato_wfn_form(mole: &Mole, grid: &[f64; 3], form: &str) -> (Vec<f64>, Vec<f64>) {
    let atos = mole.basis().atos();
    let atos = if form == "car" { &atos.cars } else { &atos.sphs };
    let atoms = mole.atoms();
    let nato = atos.len();
    let mut wfn0 = vec![0.0; nato];
    // 一阶导数按 [r * nato + i] 存放，便于按方向连续读取
    let mut wfn1 = vec![0.0; nato * 3];
    for i in 0..nato {
        let atm = atos[i].atm;
        let [ax, ay, az] = atoms.get(atm).xyz;
        let xyz = [grid[0] - ax, grid[1] - ay, grid[2] - az];
        let coes = atos[i].coes_norm();
        let alps = &atos[i].alps;
        let lmn = atos[i].lmn();
        let (val0, val1, _val2) = wfnfunc::cgf(&xyz, lmn, alps, coes, 1);
        wfn0[i] = val0;
        for r in 0..3 {
            wfn1[r * nato + i] = val1[r];
        }
    }
    (wfn0, wfn1)
}

/// 格点分块的批大小（控制单个分块的 χ 体积，约 4 MB）
fn grid_batch(nato: usize) -> usize {
    (1_000_000 / (8 * nato.max(1))).clamp(128, 4096)
}

/// 估算一个原子轨道的有效半径（Bohr）
///
/// 高斯型轨道的径向部分为 `r^l Σ c_i exp(-a_i r²)`。这里用归一化收缩系数
/// 的上界（角向因子 ≤ 1）扫描到 `|χ(r)| < tol` 的位置，作为"该轨道只在
/// 这个半径内显著"的判据；超过该半径的格点不再计入该轨道的贡献。
fn ato_extent(ato: &Ato, tol: f64) -> f64 {
    let coes = ato.coes_norm();
    let l = ato.ang;
    let mut last = 0.0f64;
    let mut r = 0.0f64;
    while r <= 40.0 {
        let rl = match l {
            0 => 1.0,
            1 => r,
            2 => r * r,
            3 => r * r * r,
            _ => r.powi(l as i32),
        };
        let mut s = 0.0f64;
        for (a, c) in ato.alps.iter().zip(coes.iter()) {
            s += c.abs() * (-a * r * r).exp();
        }
        if s * rl >= tol {
            last = r;
        }
        r += 0.1;
    }
    last + 0.1
}

/// 并行遍历格点时的最小分块长度
///
/// `fold` 会为每个分块分配一份 `nato²` 累加器；分块太小会制造大量
/// 临时分配（每轮数百 MB），所以按线程数限制分块数量。
fn grid_min_len(ngrid: usize) -> usize {
    let chunks = rayon::current_num_threads() * 8;
    (ngrid / chunks).max(64)
}

/// 由格点上的原子轨道值/梯度与密度矩阵计算 ρ、∇ρ 与 |∇ρ|²
///
/// `chi0` 按 `[p * n + i]`、`chi1` 按 `[p * 3 * n + r * n + i]` 展开，其中 `n`
/// 是参与计算的原子轨道数（[`RKS`] 用全基组，[`super::lks::LKS`] 用 π 约化基）。
/// 分块做矩阵乘法：每块先算 `T = χ · D`，再用 `T` 与 `χ`、`∇χ` 逐行点积。
pub(crate) fn grid_rho(chi0: &[f64], chi1: &[f64], ngrid: usize, n: usize, dmat: &Array2<f64>) -> (Vec<f64>, Vec<[f64; 3]>, Vec<f64>) {
    let nb = grid_batch(n);
    let nchunk = ngrid.div_ceil(nb);
    let parts: Vec<(Vec<f64>, Vec<[f64; 3]>)> = (0..nchunk)
        .into_par_iter()
        .map(|c| {
            let p0 = c * nb;
            let p1 = ((c + 1) * nb).min(ngrid);
            let m = p1 - p0;
            let chi = ndarray::ArrayView2::from_shape((m, n), &chi0[p0 * n..p1 * n]).expect("chi0 形状错误");
            let t = chi.dot(dmat);
            let mut r0 = vec![0.0f64; m];
            let mut r1 = vec![[0.0f64; 3]; m];
            for q in 0..m {
                let trow = t.row(q);
                let crow = chi.row(q);
                let i0 = (p0 + q) * n * 3;
                let mut a0 = 0.0f64;
                let mut a1 = 0.0f64;
                let mut a2 = 0.0f64;
                let mut a3 = 0.0f64;
                let gx = &chi1[i0..i0 + n];
                let gy = &chi1[i0 + n..i0 + 2 * n];
                let gz = &chi1[i0 + 2 * n..i0 + 3 * n];
                for i in 0..n {
                    let ti = trow[i];
                    let ci = crow[i];
                    a0 += ti * ci;
                    a1 += ti * gx[i];
                    a2 += ti * gy[i];
                    a3 += ti * gz[i];
                }
                r0[q] = a0;
                r1[q] = [2.0 * a1, 2.0 * a2, 2.0 * a3];
            }
            (r0, r1)
        })
        .collect();
    let mut rho0s = Vec::with_capacity(ngrid);
    let mut rho1s = Vec::with_capacity(ngrid);
    let mut sigmas = Vec::with_capacity(ngrid);
    for (r0, r1) in parts {
        for (q, v) in r0.into_iter().enumerate() {
            let [dx, dy, dz] = r1[q];
            rho0s.push(v);
            rho1s.push([dx, dy, dz]);
            sigmas.push(dx * dx + dy * dy + dz * dz);
        }
    }
    (rho0s, rho1s, sigmas)
}

/// 由格点势分量与原子轨道值/梯度构造交换相关矩阵
///
/// `V_ij = Σ_p w_p [ φ_i φ_j v_ρ + 2 v_σ ∇ρ·(φ_i ∇φ_j + φ_j ∇φ_i) ]`，
/// 分块写成矩阵乘法：
/// `V = χᵀ diag(w v_ρ) χ + 2 Σ_r [ (χ·w v_σ g_r)ᵀ ∂_rχ + 转置 ]`。
pub(crate) fn grid_xcmat(chi0: &[f64], chi1: &[f64], ngrid: usize, n: usize, weits: &[f64], vrhos: &[f64], vsigmas: &[f64], rho1s: &[[f64; 3]]) -> Array2<f64> {
    let nb = grid_batch(n);
    let nchunk = ngrid.div_ceil(nb);
    let parts: Vec<Array2<f64>> = (0..nchunk)
        .into_par_iter()
        .map(|c| {
            let p0 = c * nb;
            let p1 = ((c + 1) * nb).min(ngrid);
            let m = p1 - p0;
            let chi = ndarray::ArrayView2::from_shape((m, n), &chi0[p0 * n..p1 * n]).expect("chi0 形状错误");
            let mut v = Array2::<f64>::zeros((n, n));
            // LDA 部分
            let a: Vec<f64> = (p0..p1).map(|p| weits[p] * vrhos[p]).collect();
            let av = ndarray::ArrayView1::from(&a[..]).insert_axis(ndarray::Axis(1));
            let xa = &chi * &av;
            v += &xa.t().dot(&chi);
            // 梯度（GGA）部分
            for r in 0..3 {
                let b: Vec<f64> = (p0..p1).map(|p| weits[p] * vsigmas[p] * rho1s[p][r]).collect();
                let bv = ndarray::ArrayView1::from(&b[..]).insert_axis(ndarray::Axis(1));
                let xb = &chi * &bv;
                // 导数在 chi1 中按点存放（步长 3*n），先整理成连续矩阵
                let mut dbuf = vec![0.0f64; m * n];
                for q in 0..m {
                    let src = (p0 + q) * 3 * n + r * n;
                    dbuf[q * n..(q + 1) * n].copy_from_slice(&chi1[src..src + n]);
                }
                let dchi = ndarray::ArrayView2::from_shape((m, n), &dbuf[..]).expect("chi1 形状错误");
                let mm = xb.t().dot(&dchi);
                v.scaled_add(2.0, &mm);
                v.scaled_add(2.0, &mm.t());
            }
            v
        })
        .collect();
    let mut out = Array2::<f64>::zeros((n, n));
    for v in parts {
        out += &v;
    }
    out
}

/// 由占据轨道系数构造闭壳层密度矩阵 `D = 2 C_occ C_occᵀ`
pub(crate) fn rdmat(nele: usize, cmat: &Array2<f64>) -> Array2<f64> {
    let coef = cmat.slice(s![.., ..nele / 2]).to_owned();
    coef.dot(&coef.t()) * 2.0
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
        Self::with_grid(mole, 75, 194, func)
    }

    /// 按泛函名创建（如 `"b3lyp"`、`"pbe0"`）；未知名称返回 `None`
    pub fn from_name(mole: &'m Mole, name: &str) -> Option<Self> {
        XcFunc::from_name(name).map(|func| Self::new(mole, func))
    }

    /// 用指定的径向/球面格点数与泛函创建 RKS 计算器
    pub fn with_grid(mole: &'m Mole, radial: usize, angular: usize, func: XcFunc) -> Self {
        Self::with_grid_prune(mole, radial, angular, func, false)
    }

    /// 同上，但可以关闭径向-角向剪枝（用于收敛性检查）
    pub fn with_grid_prune(mole: &'m Mole, radial: usize, angular: usize, func: XcFunc, prune: bool) -> Self {
        let int_form = mole.basis().int_form().to_string();
        let grid_caler = gridprop::dftgrid::Calculator::new(mole);
        let (grids, weits, grid_atm, grid_r) = grid_caler.mol_grids_meta(radial, angular, prune);
        let atos = mole.basis().atos();
        let nato = if int_form == "car" { atos.cars.len() } else { atos.sphs.len() };
        let ngrid = grids.len();
        let mut chi0 = vec![0.0f64; ngrid * nato];
        let mut chi1 = vec![0.0f64; ngrid * nato * 3];
        chi0.par_chunks_mut(nato.max(1)).zip(chi1.par_chunks_mut(nato.max(1) * 3)).zip(grids.par_iter()).for_each(|((c0, c1), grid)| {
            let (v0, v1) = ato_wfn_form(mole, grid, &int_form);
            c0.copy_from_slice(&v0);
            c1.copy_from_slice(&v1);
        });
        // 每个原子轨道的归属原子与有效半径（用于格点积分按块筛选）
        let ato_list = if int_form == "car" { &atos.cars } else { &atos.sphs };
        let ao_atm: Vec<usize> = ato_list.iter().map(|a| a.atm).collect();
        let ao_r: Vec<f64> = ato_list.iter().map(|a| ato_extent(a, 1.0e-7)).collect();
        let natm = mole.atoms.len();
        let mut ao_by_atm: Vec<Vec<usize>> = vec![Vec::new(); natm];
        for (i, &a) in ao_atm.iter().enumerate() {
            ao_by_atm[a].push(i);
        }
        let mut atm_dist = vec![0.0f64; natm * natm];
        for i in 0..natm {
            for j in 0..natm {
                atm_dist[i * natm + j] = mole.atoms.dist(i, j);
            }
        }
        // 格点分块：与原子边界对齐，保证块内格点同属一个原子中心
        let nb = grid_batch(nato);
        let mut blocks: Vec<GridBlock> = Vec::new();
        let mut s = 0usize;
        while s < ngrid {
            let a = grid_atm[s];
            let mut e = s;
            while e < ngrid && grid_atm[e] == a {
                e += 1;
            }
            let mut p = s;
            while p < e {
                let q = (p + nb).min(e);
                let mut rmin = f64::INFINITY;
                let mut rmax = 0.0f64;
                for k in p..q {
                    let r = grid_r[k];
                    if r < rmin {
                        rmin = r;
                    }
                    if r > rmax {
                        rmax = r;
                    }
                }
                blocks.push(GridBlock { p0: p, p1: q, atm: a, rmin, rmax });
                p = q;
            }
            s = e;
        }
        Self {
            mole,
            func,
            int_form,
            radial,
            angular,
            ngrid,
            nato,
            max_iter: 60,
            tol: 1e-8,
            grids,
            weits,
            chi0,
            chi1,
            screen: 1.0e-11,
            max_mem: 8 * 1024 * 1024 * 1024,
            ao_cut: 0.0,
            atm_dist,
            ao_atm,
            ao_r,
            ao_by_atm,
            grid_atm,
            grid_r,
            blocks,
        }
    }

    /// 计算一个格点块中"显著"的原子轨道下标（升序）
    ///
    /// 判据：块内格点都在以 `blk.atm` 为中心、半径区间 `[rmin, rmax]` 的球壳上，
    /// 因此到原子 `b` 中心的最小距离 `min_d` 可以精确给出；若
    /// `min_d >= ao_cut`（或超过该轨道自身有效半径），该轨道在块内可忽略。
    fn block_sig(&self, blk: &GridBlock) -> Vec<usize> {
        let natm = self.ao_by_atm.len();
        let cut = self.ao_cut;
        let mut sig = Vec::with_capacity(self.nato);
        for b in 0..natm {
            let run = &self.ao_by_atm[b];
            if run.is_empty() {
                continue;
            }
            let d = self.atm_dist[blk.atm * natm + b];
            let min_d = if d < blk.rmin {
                blk.rmin - d
            } else if d > blk.rmax {
                d - blk.rmax
            } else {
                0.0
            };
            for &i in run {
                if min_d < self.ao_r[i].min(cut) {
                    sig.push(i);
                }
            }
        }
        sig.sort_unstable();
        sig
    }

    /// 显著 AO 路径的统计信息：`(块数, 平均每块 AO 数, 最大每块 AO 数)`
    pub fn ao_stats(&self) -> (usize, f64, usize) {
        let mut total = 0usize;
        let mut mx = 0usize;
        for blk in &self.blocks {
            let m = self.block_sig(blk).len();
            total += m;
            mx = mx.max(m);
        }
        let n = self.blocks.len().max(1);
        (n, total as f64 / n as f64, mx)
    }

    /// 计算核-核排斥能
    pub fn nuc_eng(&self) -> f64 {
        self.mole.atoms.nuc_eng()
    }

    // 计算闭壳层密度矩阵 D = 2 * C_occ C_occ^T
    fn calc_rdmat(&self, nele: usize, cmat: &Array2<f64>) -> Array2<f64> {
        rdmat(nele, cmat)
    }

    /// 计算所有格点处的电子密度、密度梯度与 |∇ρ|²
    ///
    /// 分块做矩阵乘法：每个分块先算 `T = χ · D`（一次 GEMM），再用
    /// `T` 与 `χ`、`∇χ` 逐行点积得到 ρ 与 ∇ρ。相比逐点双重循环，
    /// 这样把 O(ngrid·nato²) 的乘法交给 BLAS 级的内核，且访存连续。
    pub fn calc_rho(&self, dmat: &Array2<f64>) -> (Vec<f64>, Vec<[f64; 3]>, Vec<f64>) {
        grid_rho(&self.chi0, &self.chi1, self.ngrid, self.nato, dmat)
    }

    // 计算交换相关矩阵
    //
    // V_ij = Σ_p w_p [ φ_i φ_j v_ρ + 2 v_σ ∇ρ·(φ_i ∇φ_j + φ_j ∇φ_i) ]
    //
    // 分块写成矩阵乘法：
    //   V = χᵀ diag(w v_ρ) χ + 2 Σ_r [ (χ·w v_σ g_r)ᵀ ∂_rχ + 转置 ]
    fn calc_xcmat(&self, vrhos: &[f64], vsigmas: &[f64], rho1s: &[[f64; 3]]) -> Array2<f64> {
        grid_xcmat(&self.chi0, &self.chi1, self.ngrid, self.nato, &self.weits, vrhos, vsigmas, rho1s)
    }

    /// 显著 AO 路径下的电子密度与梯度
    ///
    /// 与 [`Self::calc_rho`] 相同的数学，但每个格点块只取 `block_sig` 给出的
    /// 显著轨道：`ρ = Σ_{i∈sig} T_i χ_i`，`T = χ_sig · D_sig`。大分子上把
    /// 每块的 `nato²` 计算量降到 `m²`（m ≪ nato）。
    fn calc_rho_ao(&self, dmat: &Array2<f64>) -> (Vec<f64>, Vec<[f64; 3]>, Vec<f64>) {
        let nato = self.nato;
        let ngrid = self.ngrid;
        let dstd = dmat.as_standard_layout();
        let d = dstd.as_slice().expect("密度矩阵不是连续存储");
        let parts: Vec<(usize, Vec<f64>, Vec<[f64; 3]>)> = self
            .blocks
            .par_iter()
            .map(|blk| {
                let sig = self.block_sig(blk);
                let m = sig.len();
                let n = blk.p1 - blk.p0;
                let mut chi = vec![0.0f64; n * m];
                for q in 0..n {
                    let src = (blk.p0 + q) * nato;
                    for (k, &i) in sig.iter().enumerate() {
                        chi[q * m + k] = self.chi0[src + i];
                    }
                }
                let mut dsig = vec![0.0f64; m * m];
                for (ki, &i) in sig.iter().enumerate() {
                    let row = i * nato;
                    for (kj, &j) in sig.iter().enumerate() {
                        dsig[ki * m + kj] = d[row + j];
                    }
                }
                let chi_v = ndarray::ArrayView2::from_shape((n, m), &chi[..]).expect("chi0 形状错误");
                let d_v = ndarray::ArrayView2::from_shape((m, m), &dsig[..]).expect("密度子块形状错误");
                let t = chi_v.dot(&d_v);
                let mut r0 = vec![0.0f64; n];
                let mut r1 = vec![[0.0f64; 3]; n];
                for q in 0..n {
                    let cbase = q * m;
                    let gbase = (blk.p0 + q) * 3 * nato;
                    let trow = t.row(q);
                    let (mut a0, mut d0, mut d1, mut d2) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
                    for (k, &i) in sig.iter().enumerate() {
                        let ti = trow[k];
                        a0 += ti * chi[cbase + k];
                        d0 += ti * self.chi1[gbase + i];
                        d1 += ti * self.chi1[gbase + nato + i];
                        d2 += ti * self.chi1[gbase + 2 * nato + i];
                    }
                    r0[q] = a0;
                    r1[q] = [2.0 * d0, 2.0 * d1, 2.0 * d2];
                }
                (blk.p0, r0, r1)
            })
            .collect();
        let mut rho0s = vec![0.0f64; ngrid];
        let mut rho1s = vec![[0.0f64; 3]; ngrid];
        let mut sigmas = vec![0.0f64; ngrid];
        for (p0, r0, r1) in parts {
            for (q, v) in r0.into_iter().enumerate() {
                let [dx, dy, dz] = r1[q];
                rho0s[p0 + q] = v;
                rho1s[p0 + q] = [dx, dy, dz];
                sigmas[p0 + q] = dx * dx + dy * dy + dz * dz;
            }
        }
        (rho0s, rho1s, sigmas)
    }

    /// 显著 AO 路径下的交换相关矩阵
    ///
    /// 每块只在显著轨道子空间上累加 `V_sig = χᵀ diag(w v_ρ) χ + …`，再散射回
    /// 完整矩阵；块内累加器按线程折叠（`fold`/`reduce`），避免把上千个
    /// `nato×nato` 临时矩阵同时留在内存里。
    fn calc_xcmat_ao(&self, vrhos: &[f64], vsigmas: &[f64], rho1s: &[[f64; 3]]) -> Array2<f64> {
        let nato = self.nato;
        let zero = || vec![0.0f64; nato * nato];
        let acc = self
            .blocks
            .par_iter()
            .fold(zero, |mut acc, blk| {
                let sig = self.block_sig(blk);
                let m = sig.len();
                let n = blk.p1 - blk.p0;
                let mut chi = vec![0.0f64; n * m];
                for q in 0..n {
                    let src = (blk.p0 + q) * nato;
                    for (k, &i) in sig.iter().enumerate() {
                        chi[q * m + k] = self.chi0[src + i];
                    }
                }
                let chi_v = ndarray::ArrayView2::from_shape((n, m), &chi[..]).expect("chi0 形状错误");
                let mut v = vec![0.0f64; m * m];
                // LDA 部分
                let a: Vec<f64> = (blk.p0..blk.p1).map(|p| self.weits[p] * vrhos[p]).collect();
                let av = ndarray::ArrayView1::from(&a[..]).insert_axis(ndarray::Axis(1));
                let xa = &chi_v * &av;
                let t = xa.t().dot(&chi_v);
                for k in 0..m * m {
                    v[k] += t.as_slice().expect("LDA 子块非连续")[k];
                }
                // 梯度（GGA）部分
                for r in 0..3 {
                    let b: Vec<f64> = (blk.p0..blk.p1).map(|p| self.weits[p] * vsigmas[p] * rho1s[p][r]).collect();
                    let bv = ndarray::ArrayView1::from(&b[..]).insert_axis(ndarray::Axis(1));
                    let xb = &chi_v * &bv;
                    let mut dbuf = vec![0.0f64; n * m];
                    for q in 0..n {
                        let src = (blk.p0 + q) * 3 * nato + r * nato;
                        for (k, &i) in sig.iter().enumerate() {
                            dbuf[q * m + k] = self.chi1[src + i];
                        }
                    }
                    let dchi = ndarray::ArrayView2::from_shape((n, m), &dbuf[..]).expect("chi1 形状错误");
                    let mm = xb.t().dot(&dchi);
                    for ki in 0..m {
                        let row = ki * m;
                        for kj in 0..m {
                            v[row + kj] += 2.0 * (mm[(ki, kj)] + mm[(kj, ki)]);
                        }
                    }
                }
                // 散射回完整矩阵（原子轨道的下标在 sig 中升序）
                for (ki, &i) in sig.iter().enumerate() {
                    let row = i * nato;
                    let brow = ki * m;
                    for (kj, &j) in sig.iter().enumerate() {
                        acc[row + j] += v[brow + kj];
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
        Array2::from_shape_vec((nato, nato), acc).expect("XC 矩阵形状错误")
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
        let engine = crate::rsint::Engine::new(mole, &self.int_form, self.screen, self.max_mem);
        let time_int2e = t_int2e.elapsed();
        log::info!(
            "双电子积分: {} 个壳层块, {:.2} GB{}",
            engine.nblock,
            engine.mem_bytes as f64 / 1e9,
            if engine.is_direct() { "（直接 SCF：不缓存积分）" } else { "" }
        );
        let hmat = &tmat + &vmat;
        let nato = smat.shape()[0];
        assert!(nocc <= nato, "基函数数量({nato})不足以容纳 {nocc} 个占据轨道");
        let xmat = maths::matrix_nhalf(&smat);

        let rsfun = RsXcs::with_func(self.func);
        let exact_x = self.func.exact_x();
        let nuc_eng = self.nuc_eng();

        // 初猜：叠加原子密度（SAD）；不可用时退回核心哈密顿量
        let (_engs0, cmat_guess) = maths::eigh(&xmat.dot(&hmat).dot(&xmat));
        let mut cmat = xmat.dot(&cmat_guess);
        let mut dmat = self.calc_rdmat(nele, &cmat);
        if crate::qmcal::guess::use_sad_enabled() {
            if let Some(d) = crate::qmcal::guess::sad_dmat(&smat, &hmat, &crate::qmcal::guess::ao_owners(mole, &self.int_form), mole.atoms.nucs(), nele) {
                dmat = d;
            }
        }
        let mut obt_engs = Vec::new();

        let mut ele_eng = f64::INFINITY;
        let mut xc_eng = 0.0;
        let mut delta_eng = f64::INFINITY;
        let mut converged = false;
        let mut niter = 0;
        let mut diis = Diis::new(8);
        let use_ao = self.ao_cut > 0.0;
        let t_scf = Instant::now();
        let mut time_jk = Duration::ZERO;
        let mut time_grid = Duration::ZERO;
        for i in 1..=self.max_iter {
            niter = i;
            let t_grid = Instant::now();
            // 格点上的密度与梯度
            let t_prof = Instant::now();
            let (rho0s, rho1s, sigmas) = if use_ao { self.calc_rho_ao(&dmat) } else { self.calc_rho(&dmat) };
            let t_rho = t_prof.elapsed();
            let t_prof = Instant::now();
            let (zk, vrho, vsigma) = rsfun.run(&rho0s, &sigmas);
            let t_xcfun = t_prof.elapsed();
            let t_prof = Instant::now();
            let xcmat = if use_ao { self.calc_xcmat_ao(&vrho, &vsigma, &rho1s) } else { self.calc_xcmat(&vrho, &vsigma, &rho1s) };
            if std::env::var("RSWFN_PROF").is_ok() && i == 1 {
                eprintln!(
                    "PROF grid: rho {:.2} s ({} pts), xcfun {:.2} s, xcmat {:.2} s",
                    t_rho.as_secs_f64(),
                    rho0s.len(),
                    t_xcfun.as_secs_f64(),
                    t_prof.elapsed().as_secs_f64()
                );
            }
            xc_eng = (0..rho0s.len()).map(|p| zk[p] * rho0s[p] * self.weits[p]).sum();
            time_grid += t_grid.elapsed();

            // 库伦/交换矩阵与 Fock 矩阵
            // 精确交换按泛函系数 a_x 加入：F += -a_x * K/2
            let t_jk = Instant::now();
            let (jmat, kmat) = engine.jk(&dmat);
            time_jk += t_jk.elapsed();
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
            time_jk,
            time_grid,
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
    /// SCF 中 J/K 构造的累计耗时
    pub time_jk: Duration,
    /// SCF 中格点相关（密度、XC 矩阵）的累计耗时
    pub time_grid: Duration,
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
        Self::with_grid(mole, 75, 194, func)
    }

    /// 按泛函名创建（如 `"b3lyp"`、`"pbe0"`）；未知名称返回 `None`
    pub fn from_name(mole: &'m Mole, name: &str) -> Option<Self> {
        XcFunc::from_name(name).map(|func| Self::new(mole, func))
    }

    /// 用指定的径向/球面格点数与泛函创建 UKS 计算器
    pub fn with_grid(mole: &'m Mole, radial: usize, angular: usize, func: XcFunc) -> Self {
        Self::with_grid_prune(mole, radial, angular, func, false)
    }

    /// 用指定的径向/球面格点数创建 UKS 计算器，`prune` 控制径向-角向剪枝
    pub fn with_grid_prune(mole: &'m Mole, radial: usize, angular: usize, func: XcFunc, prune: bool) -> Self {
        let [nela, nelb] = mole.atoms.neles;
        assert!(nela >= nelb, "UKS 要求 α 电子数 ≥ β 电子数（{nela} < {nelb}）");
        Self {
            base: RKS::with_grid_prune(mole, radial, angular, func, prune),
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
                let wfn0 = &base.chi0[p * nato..(p + 1) * nato];
                let gb0 = p * nato * 3;
                let gx_ = &base.chi1[gb0..gb0 + nato];
                let gy_ = &base.chi1[gb0 + nato..gb0 + 2 * nato];
                let gz_ = &base.chi1[gb0 + 2 * nato..gb0 + 3 * nato];
                let mut rho_a = 0.0;
                let mut rho_b = 0.0;
                let mut grad_a = [0.0; 3];
                let mut grad_b = [0.0; 3];
                for i in 0..nato {
                    let fi = wfn0[i];
                    let gi = [gx_[i], gy_[i], gz_[i]];
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
                        let gj = [gx_[j], gy_[j], gz_[j]];
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
        let nato = self.base.nato;
        let zero = || (vec![0.0f64; nato * nato], vec![0.0f64; nato * nato]);
        let mut acc = (0..self.base.grids.len())
            .into_par_iter()
            .with_min_len(grid_min_len(self.base.grids.len()))
            .fold(zero, |(mut xa, mut xb), p| {
                let weit = self.base.weits[p];
                if weit == 0.0 {
                    return (xa, xb);
                }
                let wfn0 = &self.base.chi0[p * nato..(p + 1) * nato];
                let gb0 = p * nato * 3;
                let gx_ = &self.base.chi1[gb0..gb0 + nato];
                let gy_ = &self.base.chi1[gb0 + nato..gb0 + 2 * nato];
                let gz_ = &self.base.chi1[gb0 + 2 * nato..gb0 + 3 * nato];
                let ga = rho.grad_a[p];
                let gb = rho.grad_b[p];
                let (vra, vrb) = (vr_a[p], vr_b[p]);
                let (vaa, vbb, vab) = (vs_aa[p], vs_bb[p], vs_ab[p]);
                for i in 0..nato {
                    let fi = wfn0[i];
                    let gi = [gx_[i], gy_[i], gz_[i]];
                    if fi == 0.0 && gi == [0.0; 3] {
                        continue;
                    }
                    let irow = i * nato;
                    for j in 0..=i {
                        let fj = wfn0[j];
                        let gj = [gx_[j], gy_[j], gz_[j]];
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

    /// 显著 AO 路径下的 α/β 密度与梯度（每块只取 `block_sig` 给出的轨道）
    fn calc_rho_spin_ao(&self, dmat_a: &Array2<f64>, dmat_b: &Array2<f64>) -> SpinRho {
        let nato = self.base.nato;
        let ngrid = self.base.ngrid;
        let da_std = dmat_a.as_standard_layout();
        let da = da_std.as_slice().expect("dmat_a 不是连续存储");
        let db_std = dmat_b.as_standard_layout();
        let db = db_std.as_slice().expect("dmat_b 不是连续存储");
        let parts: Vec<(usize, Vec<f64>, Vec<[f64; 3]>, Vec<f64>, Vec<[f64; 3]>)> = self
            .base
            .blocks
            .par_iter()
            .map(|blk| {
                let sig = self.base.block_sig(blk);
                let m = sig.len();
                let n = blk.p1 - blk.p0;
                let mut chi = vec![0.0f64; n * m];
                for q in 0..n {
                    let src = (blk.p0 + q) * nato;
                    for (k, &i) in sig.iter().enumerate() {
                        chi[q * m + k] = self.base.chi0[src + i];
                    }
                }
                let mut dsa = vec![0.0f64; m * m];
                let mut dsb = vec![0.0f64; m * m];
                for (ki, &i) in sig.iter().enumerate() {
                    let row = i * nato;
                    for (kj, &j) in sig.iter().enumerate() {
                        dsa[ki * m + kj] = da[row + j];
                        dsb[ki * m + kj] = db[row + j];
                    }
                }
                let chi_v = ndarray::ArrayView2::from_shape((n, m), &chi[..]).expect("chi0 形状错误");
                let ta = chi_v.dot(&ndarray::ArrayView2::from_shape((m, m), &dsa[..]).expect("密度子块形状错误"));
                let tb = chi_v.dot(&ndarray::ArrayView2::from_shape((m, m), &dsb[..]).expect("密度子块形状错误"));
                let mut ra = vec![0.0f64; n];
                let mut ga = vec![[0.0f64; 3]; n];
                let mut rb = vec![0.0f64; n];
                let mut gb = vec![[0.0f64; 3]; n];
                for q in 0..n {
                    let cbase = q * m;
                    let gbase = (blk.p0 + q) * 3 * nato;
                    let ta_row = ta.row(q);
                    let tb_row = tb.row(q);
                    let (mut a0, mut ax, mut ay, mut az) = (0.0, 0.0, 0.0, 0.0);
                    let (mut b0, mut bx, mut by, mut bz) = (0.0, 0.0, 0.0, 0.0);
                    for (k, &i) in sig.iter().enumerate() {
                        let ca = ta_row[k];
                        a0 += ca * chi[cbase + k];
                        ax += ca * self.base.chi1[gbase + i];
                        ay += ca * self.base.chi1[gbase + nato + i];
                        az += ca * self.base.chi1[gbase + 2 * nato + i];
                        let cb = tb_row[k];
                        b0 += cb * chi[cbase + k];
                        bx += cb * self.base.chi1[gbase + i];
                        by += cb * self.base.chi1[gbase + nato + i];
                        bz += cb * self.base.chi1[gbase + 2 * nato + i];
                    }
                    ra[q] = a0;
                    ga[q] = [2.0 * ax, 2.0 * ay, 2.0 * az];
                    rb[q] = b0;
                    gb[q] = [2.0 * bx, 2.0 * by, 2.0 * bz];
                }
                (blk.p0, ra, ga, rb, gb)
            })
            .collect();
        let mut rho_a = vec![0.0f64; ngrid];
        let mut rho_b = vec![0.0f64; ngrid];
        let mut grad_a = vec![[0.0f64; 3]; ngrid];
        let mut grad_b = vec![[0.0f64; 3]; ngrid];
        let mut saa = vec![0.0f64; ngrid];
        let mut sbb = vec![0.0f64; ngrid];
        let mut sab = vec![0.0f64; ngrid];
        for (p0, ra, ga, rb, gb) in parts {
            for (q, v) in ra.into_iter().enumerate() {
                let p = p0 + q;
                let ga_q = ga[q];
                let gb_q = gb[q];
                rho_a[p] = v;
                rho_b[p] = rb[q];
                grad_a[p] = ga_q;
                grad_b[p] = gb_q;
                saa[p] = ga_q[0] * ga_q[0] + ga_q[1] * ga_q[1] + ga_q[2] * ga_q[2];
                sbb[p] = gb_q[0] * gb_q[0] + gb_q[1] * gb_q[1] + gb_q[2] * gb_q[2];
                sab[p] = ga_q[0] * gb_q[0] + ga_q[1] * gb_q[1] + ga_q[2] * gb_q[2];
            }
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

    /// 显著 AO 路径下的 α/β 交换相关矩阵
    fn calc_xcmat_spin_ao(&self, rho: &SpinRho, vr_a: &[f64], vr_b: &[f64], vs_aa: &[f64], vs_bb: &[f64], vs_ab: &[f64]) -> (Array2<f64>, Array2<f64>) {
        let nato = self.base.nato;
        let zero = || (vec![0.0f64; nato * nato], vec![0.0f64; nato * nato]);
        let acc = self
            .base
            .blocks
            .par_iter()
            .fold(zero, |(mut acca, mut accb), blk| {
                let sig = self.base.block_sig(blk);
                let m = sig.len();
                let n = blk.p1 - blk.p0;
                let mut chi = vec![0.0f64; n * m];
                for q in 0..n {
                    let src = (blk.p0 + q) * nato;
                    for (k, &i) in sig.iter().enumerate() {
                        chi[q * m + k] = self.base.chi0[src + i];
                    }
                }
                let chi_v = ndarray::ArrayView2::from_shape((n, m), &chi[..]).expect("chi0 形状错误");
                let mut va = vec![0.0f64; m * m];
                let mut vb = vec![0.0f64; m * m];
                // LDA 部分
                let wa: Vec<f64> = (blk.p0..blk.p1).map(|p| self.base.weits[p] * vr_a[p]).collect();
                let wb: Vec<f64> = (blk.p0..blk.p1).map(|p| self.base.weits[p] * vr_b[p]).collect();
                let wa_axis = ndarray::ArrayView1::from(&wa[..]).insert_axis(ndarray::Axis(1));
                let wb_axis = ndarray::ArrayView1::from(&wb[..]).insert_axis(ndarray::Axis(1));
                let ma = (&chi_v * &wa_axis).t().dot(&chi_v);
                let mb = (&chi_v * &wb_axis).t().dot(&chi_v);
                for k in 0..m * m {
                    va[k] += ma.as_slice().expect("LDA 子块非连续")[k];
                    vb[k] += mb.as_slice().expect("LDA 子块非连续")[k];
                }
                // 梯度（GGA）部分
                for r in 0..3 {
                    let wga: Vec<f64> = (blk.p0..blk.p1).map(|p| self.base.weits[p] * (2.0 * vs_aa[p] * rho.grad_a[p][r] + vs_ab[p] * rho.grad_b[p][r])).collect();
                    let wgb: Vec<f64> = (blk.p0..blk.p1).map(|p| self.base.weits[p] * (2.0 * vs_bb[p] * rho.grad_b[p][r] + vs_ab[p] * rho.grad_a[p][r])).collect();
                    let mut dbuf = vec![0.0f64; n * m];
                    for q in 0..n {
                        let src = (blk.p0 + q) * 3 * nato + r * nato;
                        for (k, &i) in sig.iter().enumerate() {
                            dbuf[q * m + k] = self.base.chi1[src + i];
                        }
                    }
                    let dchi = ndarray::ArrayView2::from_shape((n, m), &dbuf[..]).expect("chi1 形状错误");
                    let ga = (&chi_v * &ndarray::ArrayView1::from(&wga[..]).insert_axis(ndarray::Axis(1))).t().dot(&dchi);
                    let gbm = (&chi_v * &ndarray::ArrayView1::from(&wgb[..]).insert_axis(ndarray::Axis(1))).t().dot(&dchi);
                    for ki in 0..m {
                        for kj in 0..m {
                            va[ki * m + kj] += ga[(ki, kj)] + ga[(kj, ki)];
                            vb[ki * m + kj] += gbm[(ki, kj)] + gbm[(kj, ki)];
                        }
                    }
                }
                // 散射回完整矩阵
                for (ki, &i) in sig.iter().enumerate() {
                    let row = i * nato;
                    let brow = ki * m;
                    for (kj, &j) in sig.iter().enumerate() {
                        acca[row + j] += va[brow + kj];
                        accb[row + j] += vb[brow + kj];
                    }
                }
                (acca, accb)
            })
            .reduce(zero, |(mut aa, mut ba), (ab, bb)| {
                for (x, y) in aa.iter_mut().zip(ab.iter()) {
                    *x += y;
                }
                for (x, y) in ba.iter_mut().zip(bb.iter()) {
                    *x += y;
                }
                (aa, ba)
            });
        (
            Array2::from_shape_vec((nato, nato), acc.0).expect("xcmat(α) 形状错误"),
            Array2::from_shape_vec((nato, nato), acc.1).expect("xcmat(β) 形状错误"),
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
        let engine = crate::rsint::Engine::new(mole, &self.base.int_form, self.base.screen, self.base.max_mem);
        let time_int2e = t_int2e.elapsed();
        log::info!(
            "双电子积分: {} 个壳层块, {:.2} GB{}",
            engine.nblock,
            engine.mem_bytes as f64 / 1e9,
            if engine.is_direct() { "（直接 SCF：不缓存积分）" } else { "" }
        );
        let hmat = &tmat + &vmat;
        let nato = smat.shape()[0];
        assert!(nela <= nato && nelb <= nato, "基函数数量({nato})不足以容纳占据轨道");
        let xmat = maths::matrix_nhalf(&smat);

        let rsfun = RsXcs::with_func(self.base.func);
        let exact_x = self.base.func.exact_x();
        let nuc_eng = self.base.nuc_eng();

        // 初猜：叠加原子密度（α、β 初始相同）；不可用时退回核心哈密顿量
        let (engs0, cmat_guess) = maths::eigh(&xmat.dot(&hmat).dot(&xmat));
        let cmat0 = xmat.dot(&cmat_guess);
        let mut cmat_a = cmat0.clone();
        let mut cmat_b = cmat0;
        let mut dmat_a = calc_spin_dmat(nela, &cmat_a);
        let mut dmat_b = calc_spin_dmat(nelb, &cmat_b);
        if crate::qmcal::guess::use_sad_enabled() {
            let nele = nela + nelb;
            let owners = crate::qmcal::guess::ao_owners(mole, &self.base.int_form);
            if let Some(d) = crate::qmcal::guess::sad_dmat(&smat, &hmat, &owners, mole.atoms.nucs(), nele) {
                // 闭壳层极限下 α、β 各占一半；开壳层则按各自电子数等比缩放
                if nela == nelb {
                    dmat_a = &d * 0.5;
                    dmat_b = &d * 0.5;
                } else {
                    dmat_a = &d * (nela as f64 / nele as f64);
                    dmat_b = &d * (nelb as f64 / nele as f64);
                }
            }
        }
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
        let mut time_jk = Duration::ZERO;
        let mut time_grid = Duration::ZERO;
        let use_ao = self.base.ao_cut > 0.0;

        for i in 1..=self.base.max_iter {
            niter = i;
            let dmat = &dmat_a + &dmat_b;
            let t_grid = Instant::now();
            let rho = if use_ao { self.calc_rho_spin_ao(&dmat_a, &dmat_b) } else { self.calc_rho_spin(&dmat_a, &dmat_b) };
            let (zk, vr_a, vr_b, vs_aa, vs_bb, vs_ab) = rsfun.run_spin(&rho.rho_a, &rho.rho_b, &rho.saa, &rho.sbb, &rho.sab);
            let (xc_a, xc_b) = if use_ao {
                self.calc_xcmat_spin_ao(&rho, &vr_a, &vr_b, &vs_aa, &vs_bb, &vs_ab)
            } else {
                self.calc_xcmat_spin(&rho, &vr_a, &vr_b, &vs_aa, &vs_bb, &vs_ab)
            };
            xc_eng = (0..rho.rho_a.len()).map(|p| zk[p] * (rho.rho_a[p] + rho.rho_b[p]) * self.base.weits[p]).sum();
            time_grid += t_grid.elapsed();

            // 库伦项用总密度，交换项按自旋分别处理；精确交换按 a_x 加入
            let t_jk = Instant::now();
            let (jmat, kmat_a, kmat_b) = engine.jk_spin(&dmat, &dmat_a, &dmat_b);
            time_jk += t_jk.elapsed();
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
            time_jk,
            time_grid,
        }
    }
}
