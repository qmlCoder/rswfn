//! 电子静电势的解析计算
//!
//! 参考 Multiwfn 3.8 (`function.f90` 中的 `eleesp1` / `cubesp`)。电子对静电势的
//! 贡献不做实空间格点积分 `∫ρ(r')/|r-r'|d³r'`，而是直接对密度矩阵做解析积分：
//!
//! ```text
//! V_ele(r) = -Σ_{μν} D_{μν} ⟨μ|1/|r-r'||ν⟩
//! ```
//!
//! 其中三中心积分 `⟨a|1/|r_C||b⟩` 用 McMurchie-Davidson 展开 + Boys 函数计算，
//! 与 Multiwfn 的 `eleesp1` 等价。相比格点积分有三个好处：
//!
//! - 没有格点积分误差，不需要很密的积分格点；
//! - 没有 `1/|r-r'|` 的自作用奇点，查询点可以任意取；
//! - 总电子数严格守恒，分子远端不会出现虚假的单极势。
//!
//! 返回值为原子单位（Hartree），且已含电子的负号。
//!
//! 参考: Helgaker, Jørgensen, Olsen, *Molecular Electronic-Structure Theory*,
//! 第 9 章; Cook, *Handbook of Computational Quantum Chemistry*, p237-p281.
//!
//! 完整推导见 `docs/静电势.md`。

use crate::base::Mole;
use ndarray::Array2;
use rayon::prelude::*;
use std::f64::consts::PI;

/// Boys 函数最大阶数，需覆盖 h 型 GTF 的两两叠加（5+5=10）
const MAX_FN: usize = 12;
/// 单个方向上 A 系数的最大个数（角动量之和为 10 时需 91 个）
const MAX_COF: usize = 96;

/// 未归一化的原生高斯函数
///
/// 归一化系数已经并入收缩系数，因此这里只保留指数、角动量和中心。
#[derive(Clone, Copy)]
pub struct Prim {
    pub lmn: [usize; 3],
    pub alp: f64,
    pub xyz: [f64; 3],
}

/// 一个 primitive 对上与查询点无关的量
struct PairData {
    /// 重叠中心 P
    p: [f64; 3],
    /// α_a + α_b
    ep: f64,
    /// 2π/ep * K_ab，其中 K_ab = exp(-α_a α_b/ep * |A-B|²)
    prefac: f64,
    /// 角动量之和，也是 Boys 函数的最高阶
    maxfn: usize,
    /// 三个方向上的 A 系数（与查询点无关的部分）
    axes: [AxisCoef; 3],
}

/// 单个方向上的 A 系数
///
/// 查询点确定后 `A_k = coef_k * PC^(exp_k)`，`map_k` 用于查询 Boys 函数 `F_{map_k}`。
struct AxisCoef {
    coef: Vec<f64>,
    exp: Vec<u8>,
    map: Vec<u8>,
}

impl AxisCoef {
    /// 预先算出与查询点无关的部分（Multiwfn 中 `Alri`/`Amsj`/`Antk` 的系数）
    fn new(l1: usize, l2: usize, pa: f64, pb: f64, gamma: f64) -> Self {
        let mut coef = Vec::new();
        let mut exp = Vec::new();
        let mut map = Vec::new();
        for l in 0..=(l1 + l2) {
            let tl = if l % 2 == 1 { -1.0 } else { 1.0 };
            let fjtmp = fj(l, l1, l2, pa, pb) * tl * fact(l);
            for r in 0..=(l / 2) {
                for i in 0..=((l - 2 * r) / 2) {
                    // 对应 Multiwfn 的 getAfac，但把 PC 的幂次单独拆出来
                    let ti = if i % 2 == 1 { -1.0 } else { 1.0 };
                    let pw = l - 2 * r - 2 * i;
                    let c = fjtmp * ti * (0.25 / gamma).powi((r + i) as i32) / (fact(r) * fact(i) * fact(pw));
                    coef.push(c);
                    exp.push(pw as u8);
                    map.push((l - 2 * r - i) as u8);
                }
            }
        }
        Self { coef, exp, map }
    }

    /// 查询点给定后求 A 系数
    #[inline]
    fn eval(&self, pc: f64, out: &mut [f64; MAX_COF]) {
        for k in 0..self.coef.len() {
            out[k] = self.coef[k] * pc.powi(self.exp[k] as i32);
        }
    }
}

impl PairData {
    fn new(ga: &Prim, gb: &Prim) -> Self {
        let ep = ga.alp + gb.alp;
        let mut p = [0.0; 3];
        let mut pa = [0.0; 3];
        let mut pb = [0.0; 3];
        let mut sq_ab = 0.0;
        for d in 0..3 {
            let (xa, xb) = (ga.xyz[d], gb.xyz[d]);
            p[d] = (ga.alp * xa + gb.alp * xb) / ep;
            pa[d] = p[d] - xa;
            pb[d] = p[d] - xb;
            sq_ab += (xa - xb) * (xa - xb);
        }
        let prefac = 2.0 * PI / ep * f64::exp(-ga.alp * gb.alp / ep * sq_ab);
        Self {
            p,
            ep,
            prefac,
            maxfn: ga.lmn.iter().sum::<usize>() + gb.lmn.iter().sum::<usize>(),
            axes: [
                AxisCoef::new(ga.lmn[0], gb.lmn[0], pa[0], pb[0], ep),
                AxisCoef::new(ga.lmn[1], gb.lmn[1], pa[1], pb[1], ep),
                AxisCoef::new(ga.lmn[2], gb.lmn[2], pa[2], pb[2], ep),
            ],
        }
    }

    /// 计算 primitive 对的核引积分 ⟨a|1/|r-r_C||b⟩
    fn attraction(&self, rc: &[f64; 3], al: &mut [f64; MAX_COF], am: &mut [f64; MAX_COF], an: &mut [f64; MAX_COF]) -> f64 {
        // 重叠中心到查询点的矢量 PC = P - C
        let mut pc = [0.0; 3];
        let mut sq_pc = 0.0;
        for d in 0..3 {
            pc[d] = self.p[d] - rc[d];
            sq_pc += pc[d] * pc[d];
        }

        let x = self.ep * sq_pc;
        let expnx = f64::exp(-x);
        let maxfn = self.maxfn;
        let mut fnv = [0.0f64; MAX_FN + 1];
        fnv[maxfn] = fmch(maxfn, x, expnx);
        // 由高阶向低阶递推，数值稳定
        let mut nu = maxfn;
        while nu > 0 {
            fnv[nu - 1] = (expnx + 2.0 * x * fnv[nu]) / (2 * nu - 1) as f64;
            nu -= 1;
        }

        self.axes[0].eval(pc[0], al);
        self.axes[1].eval(pc[1], am);
        self.axes[2].eval(pc[2], an);

        let mapl = &self.axes[0].map;
        let mapm = &self.axes[1].map;
        let mapn = &self.axes[2].map;

        let mut term = 0.0;
        for l in 0..mapl.len() {
            let av = al[l];
            if av == 0.0 {
                continue;
            }
            for m in 0..mapm.len() {
                let bv = am[m];
                if bv == 0.0 {
                    continue;
                }
                let ab = av * bv;
                for n in 0..mapn.len() {
                    term += ab * an[n] * fnv[mapl[l] as usize + mapm[m] as usize + mapn[n] as usize];
                }
            }
        }
        term * self.prefac
    }
}

/// 电子静电势的解析计算器
///
/// 构造时把密度矩阵展开到 primitive 基并预计算与查询点无关的量，
/// 之后可以对任意多个点反复求值。
pub struct EspCalculator {
    /// 上三角存储的 primitive 密度矩阵，与 `pairs` 同序
    dprim: Vec<f64>,
    /// 每个 primitive 对 (a >= b) 的预计算量
    pairs: Vec<(usize, usize, PairData)>,
}

impl EspCalculator {
    /// 根据分子与笛卡尔形式的密度矩阵构建计算器
    ///
    /// `dmat` 必须是 `mole.coefs().get_dmat("car", mole)` 返回的矩阵，
    /// 其指标与 `mole.basis().atos().cars` 一一对应。
    pub fn new(mole: &Mole, dmat: &Array2<f64>) -> Self {
        let atos = &mole.basis().atos().cars;
        let xyzs = mole.atoms().xyzs();
        let nato = dmat.nrows();

        // 展开成 primitive 基，并记录每个收缩轨道的起始下标
        let mut prims: Vec<Prim> = Vec::new();
        let mut offs = Vec::with_capacity(nato);
        for ato in atos.iter().take(nato) {
            offs.push(prims.len());
            let lmn = *ato.lmn();
            let xyz = xyzs[ato.atm];
            for &alp in ato.alps.iter() {
                prims.push(Prim { lmn, alp, xyz });
            }
        }
        let nprim = prims.len();

        // 上三角（含对角）的 primitive 密度
        let mut dprim = vec![0.0; nprim * (nprim + 1) / 2];
        let idx = |a: usize, b: usize| {
            if a >= b { a * (a + 1) / 2 + b } else { b * (b + 1) / 2 + a }
        };
        for i in 0..nato {
            let ci = atos[i].coes_norm();
            for j in i..nato {
                let dij = dmat[(i, j)];
                if dij == 0.0 {
                    continue;
                }
                let cj = atos[j].coes_norm();
                for (a, &ca) in ci.iter().enumerate() {
                    for (b, &cb) in cj.iter().enumerate() {
                        if i == j && b < a {
                            continue;
                        }
                        dprim[idx(offs[i] + a, offs[j] + b)] += dij * ca * cb;
                    }
                }
            }
        }

        // 只为密度非零的对预先算好几何量
        let mut pairs = Vec::with_capacity(nprim * (nprim + 1) / 2);
        for a in 0..nprim {
            for b in 0..=a {
                if dprim[a * (a + 1) / 2 + b] == 0.0 {
                    continue;
                }
                pairs.push((a, b, PairData::new(&prims[a], &prims[b])));
            }
        }

        Self { dprim, pairs }
    }

    /// 一点处电子对静电势的贡献（原子单位，含负号，通常 ≤ 0）
    pub fn val(&self, rc: &[f64; 3]) -> f64 {
        let mut al = [0.0f64; MAX_COF];
        let mut am = [0.0f64; MAX_COF];
        let mut an = [0.0f64; MAX_COF];
        let mut sum = 0.0;
        for (a, b, pair) in self.pairs.iter() {
            let d = self.dprim[a * (a + 1) / 2 + b];
            let v = pair.attraction(rc, &mut al, &mut am, &mut an);
            sum += if a == b { d * v } else { 2.0 * d * v };
        }
        -sum
    }

    /// 批量计算（按查询点并行）
    pub fn vals(&self, rcs: &[[f64; 3]]) -> Vec<f64> {
        rcs.par_iter().map(|rc| self.val(rc)).collect()
    }
}

/// `j!`（0 ≤ j ≤ 20）
#[inline]
fn fact(j: usize) -> f64 {
    const FACTS: [f64; 21] = [
        1.0,
        1.0,
        2.0,
        6.0,
        24.0,
        120.0,
        720.0,
        5040.0,
        40320.0,
        362880.0,
        3628800.0,
        39916800.0,
        479001600.0,
        6227020800.0,
        87178291200.0,
        1307674368000.0,
        20922789888000.0,
        355687428096000.0,
        6402373705728000.0,
        121645100408832000.0,
        2432902008176640000.0,
    ];
    FACTS[j]
}

/// 参考 Multiwfn 的 `fj`（Cook book p237）
fn fj(j: usize, l: usize, m: usize, aa: f64, bb: f64) -> f64 {
    let imax = j.min(l);
    let imin = j.saturating_sub(m);
    let mut sum = 0.0;
    for k in imin..=imax {
        let pre = fact(l) / fact(l - k) / fact(k) * fact(m) / fact(m + k - j) / fact(j - k);
        let at = if l - k != 0 { aa.powi((l - k) as i32) } else { 1.0 };
        let bt = if m + k - j != 0 { bb.powi((m + k - j) as i32) } else { 1.0 };
        sum += pre * at * bt;
    }
    sum
}

/// Boys 函数 F_m(x) = ∫₀¹ t^(2m) exp(-x t²) dt
///
/// 与 Multiwfn 的 `Fmch` 相同：x 较小时用带 exp(-x) 因子的级数展开，
/// x 较大时用渐近展开，避免相消误差。
fn fmch(m: usize, x: f64, expnx: f64) -> f64 {
    // Multiwfn 用 1e-8；这里收紧到 1e-12，多几次迭代即可换来接近机器精度的 Boys 函数
    const EPS: f64 = 1e-12;
    if x <= 10.0 {
        if expnx == 0.0 {
            return 0.0;
        }
        let mut a = m as f64 + 0.5;
        let mut term = 1.0 / a;
        let mut partsum = term;
        for _ in 2..=50 {
            a += 1.0;
            term *= x / a;
            partsum += term;
            if term / partsum < EPS {
                return 0.5 * partsum * expnx;
            }
        }
        return 0.5 * partsum * expnx;
    }

    // 渐近展开: F_m(x) = (2m-1)!!/(2x)^m * √(π/x)/2 - exp(-x)/(2x) * Σ...
    let mut a = m as f64 - 0.5;
    let mut b = m as f64 + 0.5;
    let xd = 1.0 / x;
    let mut approx = 0.5 * PI.sqrt() * xd.sqrt() * xd.powi(m as i32);
    for _ in 1..=m {
        b -= 1.0;
        approx *= b;
    }
    let fimult = 0.5 * expnx * xd;
    if fimult == 0.0 {
        return approx;
    }
    let fiprop = fimult / approx;
    let mut term = 1.0;
    let mut partsum = 1.0;
    let notrms = (x + m as f64) as usize;
    for _ in 2..=notrms {
        term *= a * xd;
        partsum += term;
        if (term * fiprop / partsum).abs() < EPS {
            return approx - fimult * partsum;
        }
        a -= 1.0;
    }
    approx - fimult * partsum
}

// TODO: 手写积分后端已存档到 docs/手写积分后端.md 并从代码移除；
// 原先用于对照手写核引/静电势的测试依赖 crate::rsint::hand，暂不再编译。
// 后续如需要，用 rsint 后端重写这类交叉验证测试。
