//! 交换相关泛函计算
//!
//! 泛函的具体实现由 `use-rsxcs` 特性选择:
//! - 启用: 使用外部 `rsxcs` crate (libxc 7.1.2 子集的 Rust 实现), 并可直接使用
//!   其支持的全部泛函（按 libxc 名称或编号选择，见 [`XcFunc::External`]）;
//! - 关闭: 使用内置的 `DFTxclib` 翻译（`hand`）。
//!
//! 两种后端通过本模块导出的统一接口（[`XcComponent`]、[`RksEval`]、[`UksEval`]、
//! [`XcFunc`]、[`RsXcs`]、[`eval_rks`]、[`eval_uks`]）对外提供，上层代码不感知
//! 后端差异。后端实现分别在 `hand.rs`（内部）与 `rsxcs.rs`（外部库）中。

use rayon::prelude::*;

/// 泛函逐点求值的并行分块大小（点）
const CHUNK_PTS: usize = 4096;

mod hand;
#[cfg(feature = "use-rsxcs")]
mod rsxcs;

/// 统一接口：单个分量的闭壳层/开壳层求值入口（由内置 `DFTxclib` 后端实现）
#[cfg(not(feature = "use-rsxcs"))]
pub use hand::{eval_rks, eval_uks};
/// 统一接口：单个分量的闭壳层/开壳层求值入口（由外部 `rsxcs` 后端实现）
#[cfg(feature = "use-rsxcs")]
pub use rsxcs::{eval_rks, eval_uks};

/// 本项目使用的交换/相关分量
///
/// 内部 `DFTxclib` 与外部 `rsxcs` 0.2 都通过同一组 [`eval_rks`] /
/// [`eval_uks`] 接口求值，因此上层代码不感知后端差异。
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum XcComponent {
    /// Slater 交换（LDA）
    LdaX,
    /// VWN5 相关（LDA）
    Vwn5,
    /// VWN3 相关（LDA，B3LYP 使用）
    Vwn3,
    /// Becke88 交换梯度修正（GGA）
    B88,
    /// LYP 相关（GGA）
    Lyp,
    /// Perdew86 相关梯度修正（GGA）
    P86,
    /// PBE 交换（GGA）
    PbeX,
    /// PBE 相关（GGA）
    PbeC,
}

impl XcComponent {
    /// 全部分量
    pub const ALL: [XcComponent; 8] = [
        XcComponent::LdaX,
        XcComponent::Vwn5,
        XcComponent::Vwn3,
        XcComponent::B88,
        XcComponent::Lyp,
        XcComponent::P86,
        XcComponent::PbeX,
        XcComponent::PbeC,
    ];

    /// 是否为 GGA 分量（需要 sigma）
    pub fn is_gga(self) -> bool {
        matches!(self, Self::B88 | Self::Lyp | Self::P86 | Self::PbeX | Self::PbeC)
    }
}

/// 闭壳层下单个分量的求值结果
pub struct RksEval {
    /// 每粒子能量
    pub exc: Vec<f64>,
    /// dE/dρ
    pub vrho: Vec<f64>,
    /// dE/dσ，σ = |∇ρ|²
    pub vsigma: Vec<f64>,
}

/// 开壳层下单个分量的求值结果
pub struct UksEval {
    /// 每粒子能量
    pub exc: Vec<f64>,
    /// dE/dρ_α
    pub vrho_a: Vec<f64>,
    /// dE/dρ_β
    pub vrho_b: Vec<f64>,
    /// dE/dσ_αα
    pub vsigma_aa: Vec<f64>,
    /// dE/dσ_ββ
    pub vsigma_bb: Vec<f64>,
    /// dE/dσ_αβ
    pub vsigma_ab: Vec<f64>,
}

/// 常用交换相关泛函
///
/// 每个泛函由若干 DFT 分量（交换 + 相关）与精确交换系数 `a_x` 组成；
/// `a_x` 由 KS 模块在构建 Fock 矩阵时单独加入（闭壳层按 `a_x/2` 缩放 K，
/// 开壳层按 `a_x` 直接缩放 K_α/K_β）。
///
/// 启用 `use-rsxcs` 时，还可通过 [`XcFunc::from_name`] 直接使用外部 `rsxcs`
/// 支持的任意 libxc 泛函（如 `"gga_x_pw91"`、`"hyb_gga_xc_b3lyp"`）；这类泛函
/// 由 [`XcFunc::External`] 承载，并作为整体交给外部库求值。
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum XcFunc {
    /// Slater 交换 + VWN5 相关（LDA）
    Svwn5,
    /// Becke88 交换 + LYP 相关
    BLYP,
    /// Becke88 交换 + Perdew86 相关
    BP86,
    /// PBE 交换 + PBE 相关
    PBE,
    /// PBE0（25% 精确交换 + 75% PBE 交换 + PBE 相关）
    PBE0,
    /// B3LYP（20% 精确交换；局部相关用 VWN3，与 Gaussian 一致）
    B3LYP,
    /// 外部 `rsxcs` 支持的任意泛函（存放 libxc 编号），仅启用 `use-rsxcs` 时可用
    External(u32),
}

impl XcFunc {
    /// 全部内置泛函
    pub const ALL: [XcFunc; 6] = [XcFunc::Svwn5, XcFunc::BLYP, XcFunc::BP86, XcFunc::PBE, XcFunc::PBE0, XcFunc::B3LYP];

    /// 从泛函名解析（忽略大小写、连字符与下划线）
    ///
    /// 先匹配内置名：`lda`/`svwn5`、`blyp`、`bp86`、`pbe`、`pbe0`（`pbe1pbe`）、
    /// `b3lyp`；未命中且启用 `use-rsxcs` 时，再按 libxc 短名（如 `gga_x_pw91`、
    /// `hyb_gga_xc_pbeh`）解析为 [`XcFunc::External`]。
    pub fn from_name(name: &str) -> Option<Self> {
        let key: String = name.chars().filter(|c| c.is_ascii_alphanumeric()).map(|c| c.to_ascii_lowercase()).collect();
        match key.as_str() {
            "lda" | "lsda" | "svwn" | "svwn5" => Some(Self::Svwn5),
            "blyp" => Some(Self::BLYP),
            "bp86" => Some(Self::BP86),
            "pbe" => Some(Self::PBE),
            "pbe0" | "pbeh" | "pbe1pbe" => Some(Self::PBE0),
            "b3lyp" => Some(Self::B3LYP),
            _ => from_external_name(name),
        }
    }

    /// 规范名称（内置为小写短名；外部泛函为 libxc 短名）
    pub fn name(&self) -> &'static str {
        match self {
            Self::Svwn5 => "svwn5",
            Self::BLYP => "blyp",
            Self::BP86 => "bp86",
            Self::PBE => "pbe",
            Self::PBE0 => "pbe0",
            Self::B3LYP => "b3lyp",
            Self::External(id) => external_name(*id),
        }
    }

    /// 精确（Hartree-Fock）交换的混合系数
    pub fn exact_x(&self) -> f64 {
        match self {
            Self::PBE0 => 0.25,
            Self::B3LYP => 0.20,
            Self::External(id) => external_exact_x(*id),
            _ => 0.0,
        }
    }

    /// 是否为杂化泛函（含精确交换）
    pub fn is_hybrid(&self) -> bool {
        self.exact_x() > 0.0
    }

    /// 若为外部 `rsxcs` 泛函，返回其 libxc 编号
    pub fn external_id(&self) -> Option<u32> {
        match self {
            Self::External(id) => Some(*id),
            _ => None,
        }
    }

    /// DFT 分量及其系数（闭壳层/开壳层共用）
    ///
    /// 外部泛函整体交给外部库求值，没有分量分解，返回空切片。
    fn parts(&self) -> &'static [(XcComponent, f64)] {
        match self {
            Self::Svwn5 => &SVWN5_PARTS,
            Self::BLYP => &BLYP_PARTS,
            Self::BP86 => &BP86_PARTS,
            Self::PBE => &PBE_PARTS,
            Self::PBE0 => &PBE0_PARTS,
            Self::B3LYP => &B3LYP_PARTS,
            Self::External(_) => &[],
        }
    }
}

/// 按 libxc 短名解析外部泛函（仅启用 `use-rsxcs` 时可用）
#[cfg(feature = "use-rsxcs")]
fn from_external_name(name: &str) -> Option<XcFunc> {
    ::rsxcs::Functional::from_name(name).map(|f| XcFunc::External(f.info().id))
}

#[cfg(not(feature = "use-rsxcs"))]
fn from_external_name(_name: &str) -> Option<XcFunc> {
    None
}

/// 外部泛函的名称（仅启用 `use-rsxcs` 时可解析）
#[cfg(feature = "use-rsxcs")]
fn external_name(id: u32) -> &'static str {
    ::rsxcs::Functional::from_id(id).map(|f| f.info().name).unwrap_or("external")
}

#[cfg(not(feature = "use-rsxcs"))]
fn external_name(_id: u32) -> &'static str {
    "external"
}

/// 杂化外部泛函的精确交换系数
///
/// `rsxcs::Info` 未提供该系数，这里按 libxc 编号给出（0.2 中仅这两个杂化泛函）。
fn external_exact_x(id: u32) -> f64 {
    match id {
        402 => 0.20, // hyb_gga_xc_b3lyp
        406 => 0.25, // hyb_gga_xc_pbeh (PBE0)
        _ => 0.0,
    }
}

/// 当前后端原生支持的泛函列表
///
/// - 内置后端（未启用 `use-rsxcs`）：返回内置的 6 个泛函（[`XcFunc::ALL`]）；
/// - 外部后端（启用 `use-rsxcs`）：返回外部 `rsxcs` 库支持的全部 libxc 泛函
///   （以 [`XcFunc::External`] 表示）。
///
/// 即列表只反映当前使用的后端，不与另一种后端混合。
#[cfg(feature = "use-rsxcs")]
pub fn supported() -> Vec<XcFunc> {
    external_ids().map(XcFunc::External).collect()
}

/// 当前后端原生支持的泛函列表（内置后端）
#[cfg(not(feature = "use-rsxcs"))]
pub fn supported() -> Vec<XcFunc> {
    XcFunc::ALL.to_vec()
}

/// 当前后端支持的泛函名称列表（顺序与 [`supported`] 一致）
///
/// 可直接用于提示、校验或命令行帮助。
pub fn supported_names() -> Vec<&'static str> {
    supported().iter().map(|f| f.name()).collect()
}

/// 外部 `rsxcs` 库支持的全部泛函的 libxc 编号
#[cfg(feature = "use-rsxcs")]
fn external_ids() -> impl Iterator<Item = u32> {
    ::rsxcs::Functional::ALL.iter().map(|f| f.info().id)
}

impl std::fmt::Display for XcFunc {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.name())
    }
}

const SVWN5_PARTS: [(XcComponent, f64); 2] = [(XcComponent::LdaX, 1.0), (XcComponent::Vwn5, 1.0)];

const BLYP_PARTS: [(XcComponent, f64); 2] = [(XcComponent::B88, 1.0), (XcComponent::Lyp, 1.0)];

const BP86_PARTS: [(XcComponent, f64); 2] = [(XcComponent::B88, 1.0), (XcComponent::P86, 1.0)];

const PBE_PARTS: [(XcComponent, f64); 2] = [(XcComponent::PbeX, 1.0), (XcComponent::PbeC, 1.0)];

const PBE0_PARTS: [(XcComponent, f64); 2] = [(XcComponent::PbeX, 0.75), (XcComponent::PbeC, 1.0)];

const B3LYP_PARTS: [(XcComponent, f64); 4] = [(XcComponent::LdaX, 0.08), (XcComponent::B88, 0.72), (XcComponent::Lyp, 0.81), (XcComponent::Vwn3, 0.19)];

/// 泛函计算器
///
/// 由 [`RsXcs::with_func`] 指定具体泛函；[`RsXcs::new`] 默认使用 B3LYP。
/// DFT 分量由 [`XcFunc`] 定义，精确交换系数通过 [`RsXcs::exact_x`]
/// 提供给 KS 模块。
pub struct RsXcs {
    func: XcFunc,
}

impl RsXcs {
    /// 创建默认（B3LYP）泛函计算器
    pub fn new() -> Self {
        Self::with_func(XcFunc::B3LYP)
    }

    /// 创建指定泛函的计算器
    pub fn with_func(func: XcFunc) -> Self {
        Self { func }
    }

    /// 当前使用的泛函
    pub fn func(&self) -> XcFunc {
        self.func
    }

    /// 当前泛函的精确交换系数
    pub fn exact_x(&self) -> f64 {
        self.func.exact_x()
    }

    /// 计算交换相关泛函的能量及其一阶导数
    ///
    /// 泛函由构造时指定的 [`XcFunc`] 决定（默认 B3LYP）。对于外部泛函
    /// （[`XcFunc::External`]），整份泛函由外部 `rsxcs` 库求值。
    ///
    /// 返回 `(zk, vrho, vsigma)`:
    /// - `zk`: 每个粒子的交换相关能 (exc), 与 libxc 的 `exc` 语义一致
    /// - `vrho`: 能量密度对电子密度的导数
    /// - `vsigma`: 能量密度对 `sigma = |∇ρ|²` 的导数
    pub fn run(&self, rho: &[f64], sigma: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        if let Some(id) = self.func.external_id() {
            return external_run_rks(id, rho, sigma);
        }

        let npt = rho.len();
        assert_eq!(sigma.len(), npt, "rho 与 sigma 长度不一致");
        if npt == 0 {
            return (Vec::new(), Vec::new(), Vec::new());
        }

        // 泛函的 DFT 部分（精确交换由 KS 代码按 exact_x 单独处理）。
        // 每个分量都是逐点独立的标量运算（含 exp/pow，属延迟受限），
        // 按格点分块并行可把大分子格点部分的泛函求值时间降到 1/线程数量级。
        let parts = self.func.parts();
        let chunks: Vec<(Vec<f64>, Vec<f64>, Vec<f64>)> = (0..npt.div_ceil(CHUNK_PTS))
            .into_par_iter()
            .map(|c| {
                let p0 = c * CHUNK_PTS;
                let p1 = (p0 + CHUNK_PTS).min(npt);
                let n = p1 - p0;
                let r = &rho[p0..p1];
                let s = &sigma[p0..p1];
                let mut zk = vec![0.0; n];
                let mut vrho = vec![0.0; n];
                let mut vsigma = vec![0.0; n];
                for &(comp, coef) in parts {
                    let out = eval_rks(comp, r, s);
                    for i in 0..n {
                        zk[i] += coef * out.exc[i];
                        vrho[i] += coef * out.vrho[i];
                        vsigma[i] += coef * out.vsigma[i];
                    }
                }
                (zk, vrho, vsigma)
            })
            .collect();
        let mut zk = Vec::with_capacity(npt);
        let mut vrho = Vec::with_capacity(npt);
        let mut vsigma = Vec::with_capacity(npt);
        for (a, b, c) in chunks {
            zk.extend(a);
            vrho.extend(b);
            vsigma.extend(c);
        }

        (zk, vrho, vsigma)
    }

    /// 计算自旋极化交换相关泛函的能量及其一阶导数
    ///
    /// 使用与闭壳层 [`RsXcs::run`] 相同的泛函；精确交换由调用方
    /// （UKS 模块）按 [`RsXcs::exact_x`] 额外加入。外部泛函由 `rsxcs` 库整体求值。
    ///
    /// 返回 `(zk, vrho_a, vrho_b, vsigma_aa, vsigma_bb, vsigma_ab)`：
    /// - `zk`: 每个粒子的交换相关能（exc）
    /// - `vrho_a`/`vrho_b`: 能量密度对 α/β 密度的导数
    /// - `vsigma_aa`/`vsigma_bb`/`vsigma_ab`: 能量密度对
    ///   `σ_αα = |∇ρ_α|²`、`σ_ββ = |∇ρ_β|²`、`σ_αβ = ∇ρ_α·∇ρ_β` 的导数
    #[allow(clippy::too_many_arguments)]
    pub fn run_spin(&self, rho_a: &[f64], rho_b: &[f64], saa: &[f64], sbb: &[f64], sab: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
        if let Some(id) = self.func.external_id() {
            return external_run_uks(id, rho_a, rho_b, saa, sbb, sab);
        }

        let npt = rho_a.len();
        assert_eq!(rho_b.len(), npt, "rho_a 与 rho_b 长度不一致");
        assert_eq!(saa.len(), npt, "sigma_aa 长度不一致");
        assert_eq!(sbb.len(), npt, "sigma_bb 长度不一致");
        assert_eq!(sab.len(), npt, "sigma_ab 长度不一致");

        if npt == 0 {
            return (Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new());
        }

        let parts = self.func.parts();
        type SpinChunk = (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>);
        let chunks: Vec<SpinChunk> = (0..npt.div_ceil(CHUNK_PTS))
            .into_par_iter()
            .map(|c| {
                let p0 = c * CHUNK_PTS;
                let p1 = (p0 + CHUNK_PTS).min(npt);
                let n = p1 - p0;
                let (ra, rb) = (&rho_a[p0..p1], &rho_b[p0..p1]);
                let (saa_c, sbb_c, sab_c) = (&saa[p0..p1], &sbb[p0..p1], &sab[p0..p1]);
                let mut zk = vec![0.0; n];
                let mut va = vec![0.0; n];
                let mut vb = vec![0.0; n];
                let mut vaa = vec![0.0; n];
                let mut vbb = vec![0.0; n];
                let mut vab = vec![0.0; n];
                for &(comp, coef) in parts {
                    let out = eval_uks(comp, ra, rb, saa_c, sbb_c, sab_c);
                    for i in 0..n {
                        zk[i] += coef * out.exc[i];
                        va[i] += coef * out.vrho_a[i];
                        vb[i] += coef * out.vrho_b[i];
                        vaa[i] += coef * out.vsigma_aa[i];
                        vbb[i] += coef * out.vsigma_bb[i];
                        vab[i] += coef * out.vsigma_ab[i];
                    }
                }
                (zk, va, vb, vaa, vbb, vab)
            })
            .collect();
        let mut zk = Vec::with_capacity(npt);
        let mut vrho_a = Vec::with_capacity(npt);
        let mut vrho_b = Vec::with_capacity(npt);
        let mut vs_aa = Vec::with_capacity(npt);
        let mut vs_bb = Vec::with_capacity(npt);
        let mut vs_ab = Vec::with_capacity(npt);
        for (a, b, c, d, e, f) in chunks {
            zk.extend(a);
            vrho_a.extend(b);
            vrho_b.extend(c);
            vs_aa.extend(d);
            vs_bb.extend(e);
            vs_ab.extend(f);
        }

        (zk, vrho_a, vrho_b, vs_aa, vs_bb, vs_ab)
    }
}

/// 闭壳层下整份外部泛函的求值（仅启用 `use-rsxcs` 时可用）
#[cfg(feature = "use-rsxcs")]
fn external_run_rks(id: u32, rho: &[f64], sigma: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    rsxcs::run_rks(id, rho, sigma)
}

/// 开壳层下整份外部泛函的求值（仅启用 `use-rsxcs` 时可用）
#[cfg(feature = "use-rsxcs")]
#[allow(clippy::too_many_arguments)]
fn external_run_uks(id: u32, rho_a: &[f64], rho_b: &[f64], saa: &[f64], sbb: &[f64], sab: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    rsxcs::run_uks(id, rho_a, rho_b, saa, sbb, sab)
}

/// 未启用 `use-rsxcs` 时无法求值外部泛函
#[cfg(not(feature = "use-rsxcs"))]
fn external_run_rks(_id: u32, _rho: &[f64], _sigma: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    panic!("外部泛函需要启用 use-rsxcs 特性");
}

/// 未启用 `use-rsxcs` 时无法求值外部泛函
#[cfg(not(feature = "use-rsxcs"))]
#[allow(clippy::too_many_arguments)]
fn external_run_uks(_id: u32, _rho_a: &[f64], _rho_b: &[f64], _saa: &[f64], _sbb: &[f64], _sab: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    panic!("外部泛函需要启用 use-rsxcs 特性");
}

impl Default for RsXcs {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::{RsXcs, XcComponent, XcFunc};

    /// 能量密度 = 每粒子能量 * 密度
    fn energy_density(exc: &[f64], rho: &[f64], i: usize) -> f64 {
        exc[i] * rho[i]
    }

    #[test]
    fn derivatives_match_finite_difference() {
        let rho0 = [0.5_f64, 1.0, 2.0];
        let sigma0 = [0.2_f64, 0.0, 0.8];
        let fun = RsXcs::new();
        let (_zk, vrho, vsigma) = fun.run(&rho0, &sigma0);

        // ∂(ρ·exc)/∂ρ = vrho
        for i in 0..rho0.len() {
            let h = 1.0e-6 * rho0[i].max(1.0);
            let mut rho_p = rho0;
            let mut rho_m = rho0;
            rho_p[i] += h;
            rho_m[i] -= h;
            let (zk_p, _, _) = fun.run(&rho_p, &sigma0);
            let (zk_m, _, _) = fun.run(&rho_m, &sigma0);
            let num = (energy_density(&zk_p, &rho_p, i) - energy_density(&zk_m, &rho_m, i)) / (2.0 * h);
            assert!((num - vrho[i]).abs() < 1.0e-4 * vrho[i].abs().max(1.0), "vrho[{i}]: finite-difference={num}, analytic={}", vrho[i]);
        }

        // ∂(ρ·exc)/∂σ = vsigma
        for i in 0..rho0.len() {
            if sigma0[i] == 0.0 {
                continue; // 边界点, 避免负 sigma
            }
            let h = 1.0e-6 * sigma0[i].max(1.0);
            let mut sigma_p = sigma0;
            let mut sigma_m = sigma0;
            sigma_p[i] += h;
            sigma_m[i] -= h;
            let (zk_p, _, _) = fun.run(&rho0, &sigma_p);
            let (zk_m, _, _) = fun.run(&rho0, &sigma_m);
            let num = (energy_density(&zk_p, &rho0, i) - energy_density(&zk_m, &rho0, i)) / (2.0 * h);
            assert!((num - vsigma[i]).abs() < 1.0e-4 * vsigma[i].abs().max(1.0), "vsigma[{i}]: finite-difference={num}, analytic={}", vsigma[i]);
        }
    }

    /// 自旋极化泛函在闭壳层极限（ρ_α = ρ_β）下应与未极化泛函一致
    #[test]
    fn spin_closed_shell_matches_unpolarized() {
        let rho = [0.1_f64, 0.5, 1.0, 2.0];
        let sigma = [0.05_f64, 0.2, 0.0, 0.8];
        let fun = RsXcs::new();
        let (zk0, vrho0, vsig0) = fun.run(&rho, &sigma);

        // 闭壳层: ρ_α = ρ_β = ρ/2, σ_αα = σ_ββ = σ_αβ = |∇ρ|²/4
        let rho_a: Vec<f64> = rho.iter().map(|r| 0.5 * r).collect();
        let rho_b = rho_a.clone();
        let saa: Vec<f64> = sigma.iter().map(|s| 0.25 * s).collect();
        let sbb = saa.clone();
        let sab = saa.clone();
        let (zk, vra, vrb, vaa, vbb, vab) = fun.run_spin(&rho_a, &rho_b, &saa, &sbb, &sab);

        for i in 0..rho.len() {
            assert!((zk[i] - zk0[i]).abs() < 1e-9, "zk[{i}]: {} vs {}", zk[i], zk0[i]);
            // ρ_α = ρ/2 ⇒ ∂E/∂ρ = (vrho_a + vrho_b)/2
            assert!((vra[i] + vrb[i] - 2.0 * vrho0[i]).abs() < 1e-7 * vrho0[i].abs().max(1.0), "vrho[{i}]: {} vs {}", vra[i] + vrb[i], 2.0 * vrho0[i]);
            // ∂E/∂σ = (v_αα + v_ββ + v_αβ)/4
            if sigma[i] != 0.0 {
                assert!(
                    (vaa[i] + vbb[i] + vab[i] - 4.0 * vsig0[i]).abs() < 1e-7 * vsig0[i].abs().max(1.0),
                    "vsigma[{i}]: {} vs {}",
                    vaa[i] + vbb[i] + vab[i],
                    4.0 * vsig0[i]
                );
            }
        }
    }

    /// 自旋极化泛函的解析导数应与有限差分一致
    #[test]
    fn spin_derivatives_match_finite_difference() {
        let rho_a0 = [0.4_f64, 0.8, 1.2];
        let rho_b0 = [0.1_f64, 0.4, 0.6];
        let saa0 = [0.2_f64, 0.1, 0.8];
        let sbb0 = [0.05_f64, 0.3, 0.2];
        let sab0 = [0.1_f64, 0.15, 0.4];
        let fun = RsXcs::new();

        let energy_density = |ra: &[f64], rb: &[f64], aa: &[f64], bb: &[f64], ab: &[f64], i: usize| {
            let (zk, _, _, _, _, _) = fun.run_spin(ra, rb, aa, bb, ab);
            zk[i] * (ra[i] + rb[i])
        };

        let (_, vra, vrb, vaa, vbb, vab) = fun.run_spin(&rho_a0, &rho_b0, &saa0, &sbb0, &sab0);

        for i in 0..rho_a0.len() {
            // ∂/∂ρ_α
            let h = 1e-6 * rho_a0[i].max(1.0);
            let (mut ra_p, mut ra_m) = (rho_a0, rho_a0);
            ra_p[i] += h;
            ra_m[i] -= h;
            let num = (energy_density(&ra_p, &rho_b0, &saa0, &sbb0, &sab0, i) - energy_density(&ra_m, &rho_b0, &saa0, &sbb0, &sab0, i)) / (2.0 * h);
            assert!((num - vra[i]).abs() < 1e-4 * vra[i].abs().max(1.0), "vrho_a[{i}]: fd={num}, analytic={}", vra[i]);

            // ∂/∂ρ_β
            let h = 1e-6 * rho_b0[i].max(1.0);
            let (mut rb_p, mut rb_m) = (rho_b0, rho_b0);
            rb_p[i] += h;
            rb_m[i] -= h;
            let num = (energy_density(&rho_a0, &rb_p, &saa0, &sbb0, &sab0, i) - energy_density(&rho_a0, &rb_m, &saa0, &sbb0, &sab0, i)) / (2.0 * h);
            assert!((num - vrb[i]).abs() < 1e-4 * vrb[i].abs().max(1.0), "vrho_b[{i}]: fd={num}, analytic={}", vrb[i]);

            // ∂/∂σ_αα
            let h = 1e-6 * saa0[i].max(1.0);
            let (mut s_p, mut s_m) = (saa0, saa0);
            s_p[i] += h;
            s_m[i] -= h;
            let num = (energy_density(&rho_a0, &rho_b0, &s_p, &sbb0, &sab0, i) - energy_density(&rho_a0, &rho_b0, &s_m, &sbb0, &sab0, i)) / (2.0 * h);
            assert!((num - vaa[i]).abs() < 1e-4 * vaa[i].abs().max(1.0), "vsigma_aa[{i}]: fd={num}, analytic={}", vaa[i]);

            // ∂/∂σ_ββ
            let h = 1e-6 * sbb0[i].max(1.0);
            let (mut s_p, mut s_m) = (sbb0, sbb0);
            s_p[i] += h;
            s_m[i] -= h;
            let num = (energy_density(&rho_a0, &rho_b0, &saa0, &s_p, &sab0, i) - energy_density(&rho_a0, &rho_b0, &saa0, &s_m, &sab0, i)) / (2.0 * h);
            assert!((num - vbb[i]).abs() < 1e-4 * vbb[i].abs().max(1.0), "vsigma_bb[{i}]: fd={num}, analytic={}", vbb[i]);

            // ∂/∂σ_αβ
            let h = 1e-6 * sab0[i].max(1.0);
            let (mut s_p, mut s_m) = (sab0, sab0);
            s_p[i] += h;
            s_m[i] -= h;
            let num = (energy_density(&rho_a0, &rho_b0, &saa0, &sbb0, &s_p, i) - energy_density(&rho_a0, &rho_b0, &saa0, &sbb0, &s_m, i)) / (2.0 * h);
            assert!((num - vab[i]).abs() < 1e-4 * vab[i].abs().max(1.0), "vsigma_ab[{i}]: fd={num}, analytic={}", vab[i]);
        }
    }

    /// 内部 `DFTxclib` 与外部 `rsxcs` 0.2 必须实现同一组分量接口
    #[test]
    fn components_export_identical_interface() {
        let rho = [0.5_f64, 1.0];
        let sigma = [0.2_f64, 0.8];
        let rho_a = [0.3_f64, 0.6];
        let rho_b = [0.2_f64, 0.4];
        let saa = [0.1_f64, 0.3];
        let sbb = [0.05_f64, 0.2];
        let sab = [0.02_f64, 0.1];
        for comp in XcComponent::ALL {
            let rks = super::eval_rks(comp, &rho, &sigma);
            assert_eq!(rks.exc.len(), rho.len());
            assert_eq!(rks.vrho.len(), rho.len());
            assert_eq!(rks.vsigma.len(), rho.len());

            let uks = super::eval_uks(comp, &rho_a, &rho_b, &saa, &sbb, &sab);
            assert_eq!(uks.exc.len(), rho.len());
            assert_eq!(uks.vrho_a.len(), rho.len());
            assert_eq!(uks.vrho_b.len(), rho.len());
            assert_eq!(uks.vsigma_aa.len(), rho.len());
            assert_eq!(uks.vsigma_bb.len(), rho.len());
            assert_eq!(uks.vsigma_ab.len(), rho.len());
        }
    }

    /// 所有注册的泛函都能按名字解析，且能量与导数均为有限值
    #[test]
    fn all_functionals_are_usable() {
        assert_eq!(XcFunc::from_name("B3LYP"), Some(XcFunc::B3LYP));
        assert_eq!(XcFunc::from_name("pbe-0"), Some(XcFunc::PBE0));
        assert_eq!(XcFunc::from_name("LDA"), Some(XcFunc::Svwn5));
        assert_eq!(XcFunc::from_name("not-a-func"), None);

        assert_eq!(XcFunc::B3LYP.exact_x(), 0.20);
        assert_eq!(XcFunc::PBE0.exact_x(), 0.25);
        assert!(!XcFunc::PBE.is_hybrid());
        assert_eq!(XcFunc::ALL.len(), 6);

        let rho = [0.001_f64, 0.1, 0.5, 1.0, 2.0];
        let sigma = [0.0001_f64, 0.05, 0.2, 0.8, 1.5];
        let rho_a: Vec<f64> = rho.iter().map(|r| 0.6 * r).collect();
        let rho_b: Vec<f64> = rho.iter().map(|r| 0.4 * r).collect();
        let saa: Vec<f64> = sigma.iter().map(|s| 0.3 * s).collect();
        let sbb: Vec<f64> = sigma.iter().map(|s| 0.2 * s).collect();
        let sab: Vec<f64> = sigma.iter().map(|s| 0.1 * s).collect();

        for func in XcFunc::ALL {
            let fun = RsXcs::with_func(func);
            let (zk, vrho, vsigma) = fun.run(&rho, &sigma);
            let (zs, vra, vrb, vaa, vbb, vab) = fun.run_spin(&rho_a, &rho_b, &saa, &sbb, &sab);
            for i in 0..rho.len() {
                assert!(zk[i].is_finite(), "{func}: zk[{i}] 非有限值");
                assert!(vrho[i].is_finite(), "{func}: vrho[{i}] 非有限值");
                assert!(vsigma[i].is_finite(), "{func}: vsigma[{i}] 非有限值");
                assert!(zs[i].is_finite(), "{func}: spin zk[{i}] 非有限值");
                assert!(vra[i].is_finite() && vrb[i].is_finite(), "{func}: spin vrho 非有限值");
                assert!(vaa[i].is_finite() && vbb[i].is_finite() && vab[i].is_finite(), "{func}: spin vsigma 非有限值");
            }
        }
    }

    /// 启用 `use-rsxcs` 时，`from_name` 可解析外部库支持的全部泛函
    #[cfg(feature = "use-rsxcs")]
    #[test]
    fn external_functionals_by_name() {
        // 内置名字仍映射到内置泛函（不受 External 影响）
        assert_eq!(XcFunc::from_name("b3lyp"), Some(XcFunc::B3LYP));

        // libxc 短名映射到 External
        let pw91 = XcFunc::from_name("gga_x_pw91").expect("gga_x_pw91 应可解析");
        assert_eq!(pw91.external_id(), Some(109));
        assert_eq!(pw91.name(), "gga_x_pw91");

        // 外部杂化泛函按 libxc 编号取精确交换系数
        let b3lyp = XcFunc::from_name("hyb_gga_xc_b3lyp").expect("hyb_gga_xc_b3lyp 应可解析");
        assert_eq!(b3lyp.exact_x(), 0.20);

        // 整份外部泛函求值应得到有限值
        let rho = [0.001_f64, 0.1, 0.5, 1.0];
        let sigma = [0.0001_f64, 0.05, 0.2, 0.8];
        let fun = RsXcs::with_func(pw91);
        let (zk, vrho, vsigma) = fun.run(&rho, &sigma);
        for i in 0..rho.len() {
            assert!(zk[i].is_finite() && vrho[i].is_finite() && vsigma[i].is_finite());
        }
    }

    /// 支持列表：只反映当前后端的原生泛函，且每个名字都能被 `from_name` 解析回来
    #[test]
    fn supported_lists_all_parseable_functionals() {
        let supported = super::supported();
        let names = super::supported_names();
        assert_eq!(names.len(), supported.len());
        for (func, name) in supported.iter().zip(names.iter()) {
            assert_eq!(XcFunc::from_name(name), Some(*func), "{name} 无法解析回 {func}");
        }

        // 内置后端：列表就是内置的 6 个泛函
        #[cfg(not(feature = "use-rsxcs"))]
        for func in XcFunc::ALL {
            assert!(supported.contains(&func), "缺少内置泛函 {func}");
        }

        // 外部后端：列表全部是外部泛函
        #[cfg(feature = "use-rsxcs")]
        {
            assert!(!supported.is_empty());
            assert!(supported.iter().all(|f| f.external_id().is_some()));
        }
    }
}
