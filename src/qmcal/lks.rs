//! 局部坐标系下、仅保留 π 成分的 Kohn-Sham（DFT）计算。
//!
//! 思路（MOCV：Molecular Orbital Coefficient Vector）：
//!
//! 1. 为每个原子构造局部坐标系（LCS），使局部 `x` 轴指向分子平面外法向
//!    （[`crate::props::atomprop::direction::Calculator::LCS`]）；
//! 2. 用 [`crate::rsint::RsInt::mat1e_lcs`] 得到局部坐标系的单电子积分；双电子
//!    积分则复用 [`crate::rsint::Engine`] 的壳层块（Schwarz 筛选 + 并行），
//!    逐块旋到局部坐标系后只保留 π 分量（[`crate::rsint::Engine::mat2e_local_pi`]）；
//! 3. 只保留局部 `x` 指数为奇数的原子轨道（π 成分：`p_x`、`d_xy`、`d_xz`、
//!    `f_xxx`、`f_xyz`、…），得到维度远小于全基组的 π 约化基；
//! 4. 在该约化基上做 KS 自洽场——Fock 矩阵与系数矩阵都只有 `nπ × nπ`。
//!
//! # 注意事项
//!
//! 这里得到的是 **π 子系统模型** 的能量，并 **不是** 分子总能量：
//! π 约化基丢掉了 σ 电子与核内层，若把 π 密度直接代入完整 KS 泛函会得到
//! 无意义的大负数（详见 `docs/pi电子能量.md`）。本实现的用途是以很小的
//! 矩阵代价获得 π 分子轨道及其轨道能，能量项仅供参考。

use std::collections::HashMap;
use std::time::{Duration, Instant};

use ndarray::{Array2, s};
use rayon::prelude::*;

use crate::base::{Coefs, Mole, Stm};
use crate::maths;
use crate::props::atomprop;
use crate::props::gridprop;
use crate::props::gridprop::wfnfunc;
use crate::props::orbtprop::mocv;
use crate::rsint::RsInt;
use crate::rsxcs::{RsXcs, XcFunc};

use super::diis::Diis;

/// 取方阵在给定行/列索引上的主子矩阵（π 约化）
fn submat(mat: &Array2<f64>, idx: &[usize]) -> Array2<f64> {
    let n = idx.len();
    let mut out = Array2::<f64>::zeros((n, n));
    for (a, &i) in idx.iter().enumerate() {
        for (b, &j) in idx.iter().enumerate() {
            out[(a, b)] = mat[(i, j)];
        }
    }
    out
}

/// 局部坐标系下 π-only Kohn-Sham 计算结果
pub struct LKSResult {
    /// π 子系统模型电子能量（不含核排斥）
    pub ele_eng: f64,
    /// 交换相关能
    pub xc_eng: f64,
    /// 核-核排斥能（完整分子的，仅供参考，未计入模型能量）
    pub nuc_eng: f64,
    /// 交换相关泛函
    pub func: XcFunc,
    /// 是否自洽收敛
    pub converged: bool,
    /// 实际迭代次数
    pub niter: usize,
    /// 最后一次迭代的能量变化（Hartree）
    pub delta_eng: f64,
    /// π 电子数
    pub npi_elec: usize,
    /// π 约化基大小（Fock / 系数矩阵维度）
    pub npi: usize,
    /// 完整笛卡尔基大小（对比用）
    pub nato: usize,
    /// π 约化基下的分子轨道能量（长度 `npi`）
    pub obt_engs: Vec<f64>,
    /// π 约化基下的分子轨道系数（`npi × npi`，行为局部 π 原子轨道）
    pub cmat: Array2<f64>,
    /// π 原子轨道在完整笛卡尔基中的索引
    pub pi_atos: Vec<usize>,
    /// 各原子的局部坐标系
    pub stms: HashMap<usize, Stm>,
    /// 数值积分格点数量
    pub ngrid: usize,
    /// 单电子积分耗时
    pub time_int1e: Duration,
    /// 双电子积分耗时
    pub time_int2e: Duration,
    /// SCF 迭代耗时
    pub time_scf: Duration,
}

impl LKSResult {
    /// π 子系统模型总能量（电子能量 + 完整核排斥能，仅作参考）
    pub fn total_eng(&self) -> f64 {
        self.ele_eng + self.nuc_eng
    }

    /// 占据轨道数（闭壳层，每轨道 2 电子）
    pub fn nocc(&self) -> usize {
        self.npi_elec / 2
    }

    /// 将 π 约化基下的系数展开到完整笛卡尔基
    ///
    /// 返回 `nato_car × npi` 的矩阵：第 `k` 列是第 `k` 个 π 分子轨道在全局
    /// 笛卡尔原子轨道上的系数（`C_glob = P · C_π`，其中 `P` 把局部 π 原子轨道
    /// 展开到全局基）。用于与全基组结果对照或后续分析。
    pub fn expand_coefs(&self, mole: &Mole) -> Array2<f64> {
        let shls = &mole.basis().shls;
        let uplo = shls.uplo("car");
        let nato = shls.nato("car");
        let npi = self.pi_atos.len();
        let cars = &mole.basis().atos().cars;
        // AO 索引 → 全局壳层索引（Atos 的 shl 字段是原子内壳层序号，不是全局序号）
        let mut ao2shl = vec![0usize; nato];
        for (s, [up, lo]) in uplo.iter().enumerate() {
            for idx in *up..*lo {
                ao2shl[idx] = s;
            }
        }
        // P：把每个局部 π 原子轨道展开到全局笛卡尔基（nato × npi）
        let mut proj = Array2::<f64>::zeros((nato, npi));
        let mut tcache: HashMap<usize, Array2<f64>> = HashMap::new();
        for (k, &g) in self.pi_atos.iter().enumerate() {
            let shl = ao2shl[g];
            let atm = cars[g].atm;
            let [up, lo] = uplo[shl];
            let t = tcache.entry(shl).or_insert_with(|| mocv::get_tmat_car(lo - up, &self.stms[&atm].mat));
            // 局部轨道 μ 在全局基上的展开为 T 的第 μ 行
            let mu = g - up;
            for nu in 0..(lo - up) {
                proj[(up + nu, k)] = t[(mu, nu)];
            }
        }
        proj.dot(&self.cmat)
    }

    /// 局部坐标系下 π 原子轨道的行标签（原子 + 局部轨道符号）
    ///
    /// 每个 π 原子轨道对应一个形如 `C1:PX` 的标签：`C1` 为原子（元素+序号），
    /// `PX` 为 *局部坐标系* 下的轨道符号（`PX` 即局部 x 方向，即面外 π 分量；
    /// d 壳层为 `XY`/`XZ`）。返回顺序与 [`Self::cmat`] 的行一致。
    pub fn local_labels(&self, mole: &Mole) -> Vec<String> {
        let cars = &mole.basis().atos().cars;
        self.pi_atos
            .iter()
            .map(|&k| {
                let ato = &cars[k];
                let elem = &mole.atoms.get(ato.atm).sym;
                format!("{}{}:{}", elem, ato.atm + 1, ato.sym)
            })
            .collect()
    }

    /// 打印（返回字符串）局部坐标系下的分子轨道系数表
    ///
    /// 行为 π 原子轨道（标签只含 π 成分，见 [`Self::local_labels`]），
    /// 列为 π 分子轨道，并给出轨道能。按每 8 列分块以便阅读。
    pub fn mo_table_local(&self, mole: &Mole) -> String {
        let labels = self.local_labels(mole);
        let npi = self.npi;
        let chunk = 8;
        let mut out = String::new();
        for start in (0..npi).step_by(chunk) {
            let end = (start + chunk).min(npi);
            out.push_str(&format!("{:<12}", "AO\\MO"));
            for j in start..end {
                out.push_str(&format!("{:>12}", j + 1));
            }
            out.push('\n');
            out.push_str(&format!("{:<12}", "E(Hartree)"));
            for j in start..end {
                out.push_str(&format!("{:>12.5}", self.obt_engs[j]));
            }
            out.push('\n');
            for (i, label) in labels.iter().enumerate() {
                out.push_str(&format!("{:<12}", label));
                for j in start..end {
                    out.push_str(&format!("{:>12.5}", self.cmat[(i, j)]));
                }
                out.push('\n');
            }
            out.push('\n');
        }
        out
    }

    /// 构造一个“正常（全局笛卡尔）基”的分子对象，用于导出 fch/mwfn
    ///
    /// 把 π 约化基下的系数展开回完整笛卡尔基（[`Self::expand_coefs`]），
    /// 配以 π 占据数，并使用统一为 `car` 形式的基组（与本实现使用的积分一致）。
    /// 得到的 `Mole` 可直接交给 [`crate::writer::FchWriter`] /
    /// [`crate::writer::MfnWriter`] 保存；其中分子轨道数为 `nπ`（小于基函数数），
    /// 电子数为 π 电子数。
    pub fn to_mole(&self, mole: &Mole) -> Mole {
        let basis = mole.basis().unify_form("car");
        let cars = &basis.atos().cars;
        let ato_atms: Vec<usize> = cars.iter().map(|a| a.atm).collect();
        let ato_shls: Vec<usize> = cars.iter().map(|a| a.shl).collect();
        let ato_syms: Vec<String> = cars.iter().map(|a| a.sym.clone()).collect();
        let cmat = self.expand_coefs(mole);
        let nocc = self.nocc();
        let obt_occs: Vec<f64> = (0..self.npi).map(|i| if i < nocc { 2.0 } else { 0.0 }).collect();
        let coefs = Coefs::new(ato_atms, ato_shls, ato_syms, self.obt_engs.clone(), obt_occs, cmat);
        Mole::new(mole.atoms.clone(), Some(basis), Some(coefs))
    }
}

/// 局部坐标系下 π-only Kohn-Sham 计算器
///
/// 在局部坐标系的 π 约化基上做闭壳层 KS 自洽场。相比 [`super::ks::RKS`]，
/// 电子积分只在局部坐标系计算一次后截取 π 子块，后续 Fock 构造、DIIS、
/// 对角化与格点积分都在 `nπ` 维上进行，矩阵规模显著变小。
pub struct LKS<'m> {
    /// 分子结构引用
    pub mole: &'m Mole,
    /// 交换相关泛函
    pub func: XcFunc,
    /// 径向格点数
    pub radial: usize,
    /// 球面格点数
    pub angular: usize,
    /// 最大 SCF 迭代次数
    pub max_iter: usize,
    /// 收敛阈值（相邻两次迭代的能量变化，Hartree）
    pub tol: f64,
    /// π 电子数；`None` 表示按 p 型 π 轨道数量自动估计
    pub npi_elec: Option<usize>,
    /// 各原子的局部坐标系
    pub stms: HashMap<usize, Stm>,
    /// π 原子轨道在完整笛卡尔基中的索引
    pub pi_atos: Vec<usize>,
    /// 所有格点坐标
    pub grids: Vec<[f64; 3]>,
    /// 所有格点权重
    pub weits: Vec<f64>,
    /// π 约化基大小
    pub npi: usize,
    /// Schwarz 筛选阈值（双电子积分）
    pub screen: f64,
    /// 双电子积分块允许占用的最大内存（字节），超过则退化为直接模式
    pub max_mem: usize,
    /// 各壳层的局部坐标变换矩阵（局部 = `T ·` 全局）
    pub ts: Vec<Array2<f64>>,
    /// 每个壳层内部的 π 分量：`(局部分量下标, π 约化基下标)`
    pub pi_comp: Vec<Vec<(usize, usize)>>,
    /// 各格点处 π 原子轨道的取值，按 `[p * npi + i]` 展开
    pub chi0: Vec<f64>,
    /// 各格点处 π 原子轨道的一阶导数，按 `[p * 3 * npi + r * npi + i]` 展开
    pub chi1: Vec<f64>,
}

/// 计算某个格点处局部坐标系下 π 原子轨道的取值与一阶导数
///
/// 先按笛卡尔形式求出全局原子轨道的值/梯度，再用各壳层的局部坐标变换
/// 矩阵 `T` 旋到局部坐标系（`v^l = T · v^g`），最后取出 π 成分。
fn pi_wfn_on_grid(mole: &Mole, uplo: &[[usize; 2]], ts: &[Array2<f64>], pi_atos: &[usize], grid: &[f64; 3]) -> (Vec<f64>, Vec<f64>) {
    let atos = mole.basis().atos();
    let cars = &atos.cars;
    let atoms = mole.atoms();
    let ncar = cars.len();

    let mut gv = vec![0.0; ncar];
    let mut gg = vec![[0.0; 3]; ncar];
    for k in 0..ncar {
        let atm = cars[k].atm;
        let [ax, ay, az] = atoms.get(atm).xyz;
        let xyz = [grid[0] - ax, grid[1] - ay, grid[2] - az];
        let coes = cars[k].coes_norm();
        let alps = &cars[k].alps;
        let lmn = cars[k].lmn();
        let (val0, val1, _val2) = wfnfunc::cgf(&xyz, lmn, alps, coes, 1);
        gv[k] = val0;
        gg[k] = val1;
    }

    // 逐壳层把全局值旋转到局部坐标系
    let mut lv = vec![0.0; ncar];
    let mut lg = vec![[0.0; 3]; ncar];
    for (s, [up, lo]) in uplo.iter().enumerate() {
        let t = &ts[s];
        let dim = lo - up;
        for a in 0..dim {
            let mut v = 0.0;
            let mut g = [0.0; 3];
            for b in 0..dim {
                let tab = t[(a, b)];
                if tab == 0.0 {
                    continue;
                }
                v += tab * gv[up + b];
                g[0] += tab * gg[up + b][0];
                g[1] += tab * gg[up + b][1];
                g[2] += tab * gg[up + b][2];
            }
            lv[up + a] = v;
            lg[up + a] = g;
        }
    }

    let n = pi_atos.len();
    // 一阶导数按 [r * n + i] 存放，与共享的格点积分函数保持一致
    let mut wfn0 = vec![0.0f64; n];
    let mut wfn1 = vec![0.0f64; n * 3];
    for (i, &k) in pi_atos.iter().enumerate() {
        wfn0[i] = lv[k];
        let [gx, gy, gz] = lg[k];
        wfn1[i] = gx;
        wfn1[n + i] = gy;
        wfn1[2 * n + i] = gz;
    }
    (wfn0, wfn1)
}

impl<'m> LKS<'m> {
    /// 创建新的 LKS 计算器（默认格点：100 径向 × 74 球面）
    pub fn new(mole: &'m Mole, func: XcFunc) -> Self {
        Self::with_grid(mole, 100, 74, func, None)
    }

    /// 按泛函名创建（如 `"b3lyp"`）；未知名称返回 `None`
    pub fn from_name(mole: &'m Mole, name: &str) -> Option<Self> {
        XcFunc::from_name(name).map(|func| Self::new(mole, func))
    }

    /// 用指定的径向/球面格点数、泛函与 π 电子数创建 LKS 计算器
    ///
    /// 初始化局部坐标系、π 原子轨道索引，以及所有 π 原子轨道在格点上的
    /// 取值和一阶导数。`npi_elec` 为 `None` 时在 [`Self::run`] 中自动估计。
    pub fn with_grid(mole: &'m Mole, radial: usize, angular: usize, func: XcFunc, npi_elec: Option<usize>) -> Self {
        // 1. 每个原子的局部坐标系（无坐标系时退回随机坐标系）
        let mut dir_caler = atomprop::direction::Calculator::new(mole);
        let natm = mole.atoms.len();
        let mut stms = HashMap::new();
        for atm in 0..natm {
            let stm = match dir_caler.LCS(atm, None) {
                Some(stm) => stm,
                None => {
                    log::warn!("原子 {} 没有局部坐标系，使用随机坐标系", atm + 1);
                    dir_caler.RLCS(atm)
                }
            };
            stms.insert(atm, stm);
        }

        // 2. 选出 π 对称的原子轨道：局部 x 指数为奇数
        let cars = &mole.basis().atos().cars;
        let pi_atos: Vec<usize> = (0..cars.len()).filter(|&k| cars[k].lmn()[0] % 2 == 1).collect();

        // 3. 分子数值积分格点
        let grid_caler = gridprop::dftgrid::Calculator::new(mole);
        let (grids, weits) = grid_caler.mol_grids(radial, angular);

        // 4. 各壳层的局部坐标变换矩阵（笛卡尔形式，与 RsInt::mat1e_lcs 一致）
        let uplo = mole.basis().shls.uplo("car").clone();
        let ts: Vec<Array2<f64>> = uplo
            .iter()
            .enumerate()
            .map(|(s, [up, lo])| {
                let atm = mole.basis().shls.get(s).atm;
                mocv::get_tmat_car(lo - up, &stms[&atm].mat)
            })
            .collect();

        // 5. 每个壳层内部的 π 分量：(局部分量下标, π 约化基下标)
        let nato = mole.basis().shls.nato("car");
        let mut ao2shl = vec![0usize; nato];
        for (s, [up, lo]) in uplo.iter().enumerate() {
            for idx in *up..*lo {
                ao2shl[idx] = s;
            }
        }
        let mut pi_comp: Vec<Vec<(usize, usize)>> = vec![Vec::new(); uplo.len()];
        for (pidx, &g) in pi_atos.iter().enumerate() {
            let s = ao2shl[g];
            pi_comp[s].push((g - uplo[s][0], pidx));
        }
        let npi = pi_atos.len();

        // 6. 各格点处 π 原子轨道的取值/梯度（连续存储，供共享的格点积分函数使用）
        let ngrid = grids.len();
        let nmax = npi.max(1);
        let mut chi0 = vec![0.0f64; ngrid * nmax];
        let mut chi1 = vec![0.0f64; ngrid * nmax * 3];
        chi0.par_chunks_mut(nmax).zip(chi1.par_chunks_mut(nmax * 3)).zip(grids.par_iter()).for_each(|((c0, c1), grid)| {
            let (v0, v1) = pi_wfn_on_grid(mole, &uplo, &ts, &pi_atos, grid);
            c0.copy_from_slice(&v0);
            c1.copy_from_slice(&v1);
        });

        Self {
            mole,
            func,
            radial,
            angular,
            max_iter: 60,
            tol: 1e-8,
            npi_elec,
            stms,
            pi_atos,
            grids,
            weits,
            npi,
            screen: 1.0e-11,
            max_mem: 8 * 1024 * 1024 * 1024,
            ts,
            pi_comp,
            chi0,
            chi1,
        }
    }

    /// 计算核-核排斥能（完整分子）
    pub fn nuc_eng(&self) -> f64 {
        self.mole.atoms.nuc_eng()
    }

    /// 自动估计 π 电子数
    ///
    /// 以"参与共轭的原子数量"为初猜：每个拥有 p 壳层的原子贡献一个 π 电子
    /// （如乙烯 2 个碳 → 2 个 π 电子，苯 6 个碳 → 6 个 π 电子）。注意
    /// 6-31G 等分裂价基组每个原子有多条 p 壳层，因此必须按 *原子* 去重，
    /// 而不能按 p 壳层/π 轨道个数统计。
    ///
    /// 对含杂原子孤对电子（如吡咯 N 贡献 2 个）、带电体系等情形可能不准，
    /// 此时应通过 `npi_elec` 显式指定。
    fn guess_npi_elec(&self) -> usize {
        let cars = &self.mole.basis().atos().cars;
        let mut atms: Vec<usize> = Vec::new();
        for &k in &self.pi_atos {
            if cars[k].ang == 1 && !atms.contains(&cars[k].atm) {
                atms.push(cars[k].atm);
            }
        }
        let npi = self.pi_atos.len();
        let mut nele = atms.len().clamp(1, 2 * npi);
        if nele % 2 == 1 {
            nele = nele.saturating_sub(1).max(2).min(2 * npi);
        }
        nele.max(2)
    }

    /// 计算所有格点处 π 电子密度、密度梯度与 |∇ρ|²
    fn calc_rho(&self, dmat: &Array2<f64>) -> (Vec<f64>, Vec<[f64; 3]>, Vec<f64>) {
        super::ks::grid_rho(&self.chi0, &self.chi1, self.grids.len(), self.npi, dmat)
    }

    /// 由格点上的交换相关势构造交换相关矩阵（π 约化维度）
    ///
    /// `V_ij = Σ_p w_p [ φ_i φ_j v_ρ + 2 v_σ ∇ρ·(φ_i ∇φ_j + φ_j ∇φ_i) ]`
    fn calc_xcmat(&self, vrhos: &[f64], vsigmas: &[f64], rho1s: &[[f64; 3]]) -> Array2<f64> {
        super::ks::grid_xcmat(&self.chi0, &self.chi1, self.grids.len(), self.npi, &self.weits, vrhos, vsigmas, rho1s)
    }

    /// 执行局部坐标系 π-only KS 自洽场计算
    pub fn run(&self) -> LKSResult {
        let mole = self.mole;
        let [nela, nelb] = mole.atoms.neles;
        assert_eq!(nela, nelb, "LKS 仅支持闭壳层体系：α电子数({nela}) ≠ β电子数({nelb})");
        let npi = self.pi_atos.len();
        assert!(npi > 0, "分子中没有 π 成分的原子轨道（需要平面/共轭体系）");
        let nato = mole.basis().shls.nato("car");

        // 单电子积分：局部坐标系下（按壳层旋转），再截取 π 子块
        let int_caler = RsInt::new(mole);
        let t_int1e = Instant::now();
        let [smat_l, tmat_l, vmat_l] = int_caler.mat1e_lcs(&self.stms, [1, 1, 1]);
        let time_int1e = t_int1e.elapsed();
        let smat = submat(&smat_l, &self.pi_atos);
        let hmat = submat(&(&tmat_l + &vmat_l), &self.pi_atos);

        // 双电子积分：复用 KS 的壳层块引擎（Schwarz 筛选 + 并行构建），逐块旋到
        // 局部坐标系后只保留 π 分量，直接得到 nπ 维的双电子积分
        // （不再做全维 4 指标变换，也不再需要完整 Mat2e）。
        let t_int2e = Instant::now();
        let emat = {
            let engine = crate::rsint::Engine::new(mole, "car", self.screen, self.max_mem);
            log::info!(
                "双电子积分: {} 个壳层块, {:.2} GB{}",
                engine.nblock,
                engine.mem_bytes as f64 / 1e9,
                if engine.is_direct() { "（直接模式：不缓存积分）" } else { "" }
            );
            engine.mat2e_local_pi(&self.ts, &self.pi_comp, npi)
        };
        let time_int2e = t_int2e.elapsed();

        // π 电子数
        let npi_elec = self.npi_elec.unwrap_or_else(|| self.guess_npi_elec());
        assert!(npi_elec % 2 == 0, "LKS 需要偶数个 π 电子（当前 {npi_elec}）");
        let nocc = npi_elec / 2;
        assert!(nocc <= npi, "π 轨道数({npi})不足以容纳 {nocc} 个占据轨道（π 电子数 {npi_elec}）");

        let xmat = maths::matrix_nhalf(&smat);
        let (_engs0, cmat_guess) = maths::eigh(&xmat.dot(&hmat).dot(&xmat));
        let mut cmat = xmat.dot(&cmat_guess);
        let mut dmat = super::ks::rdmat(npi_elec, &cmat);
        let mut obt_engs: Vec<f64> = Vec::new();

        let rsfun = RsXcs::with_func(self.func);
        let exact_x = self.func.exact_x();
        let nuc_eng = self.nuc_eng();

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
            let (jmat, kmat) = crate::rsint::jk_mat(&dmat, &emat);
            let fmat = &hmat + &jmat - &kmat * (exact_x * 0.5) + &xcmat;

            let eng = super::ks::rks_ele_eng(&dmat, &hmat, &jmat, &kmat, exact_x, xc_eng);
            delta_eng = (eng - ele_eng).abs();
            log::info!("LKS 迭代 {i:>3}: π 电子能量 = {eng:>18.10}  ΔE = {delta_eng:>9.2e}  E_xc = {xc_eng:>16.8}");
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
            let dmat_step = super::ks::rdmat(npi_elec, &cmat_step);
            obt_engs = engs.to_vec();
            cmat = cmat_step;
            if extrap.is_some() {
                dmat = dmat_step;
            } else {
                dmat = &dmat * 0.5 + &dmat_step * 0.5;
            }
        }
        let time_scf = t_scf.elapsed();

        // 检查 π 电子数
        let dmat_fin = super::ks::rdmat(npi_elec, &cmat);
        let nele_fin = (&dmat_fin * &smat).sum();
        if (nele_fin - npi_elec as f64).abs() > 1.0e-4 {
            eprintln!("警告: π 密度矩阵电子数 {nele_fin} 与预期 {npi_elec} 不符");
        }

        LKSResult {
            ele_eng,
            xc_eng,
            nuc_eng,
            func: self.func,
            converged,
            niter,
            delta_eng,
            npi_elec,
            npi,
            nato,
            obt_engs,
            cmat,
            pi_atos: self.pi_atos.clone(),
            stms: self.stms.clone(),
            ngrid: self.grids.len(),
            time_int1e,
            time_int2e,
            time_scf,
        }
    }
}
