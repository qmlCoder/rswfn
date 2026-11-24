//! 内置（手写）交换相关泛函后端。
//!
//! - [`DFTxclib`] 是 `DFTxclib.F` 的 Rust 直译，仅在未启用 `use-rsxcs` 时编译；
//! - 内部 VWN3 由两种后端共用：外部 `rsxcs` 的 `lda_c_vwn_3` 与 Gaussian 采用
//!   的自旋插值不一致，为保证 B3LYP 能量不随特性开关变化，VWN3 一律用本实现。

use std::f64::consts::PI;

use super::{RksEval, UksEval};

#[cfg(not(feature = "use-rsxcs"))]
#[allow(non_snake_case)]
#[path = "DFTxclib.rs"]
mod DFTxclib;

// ---------------------------------------------------------------------------
// 内部 VWN3（两种后端共用）
// ---------------------------------------------------------------------------

/// VWN 参数集（四个参数）
#[derive(Clone, Copy)]
struct VwnParams {
    a: f64,
    b: f64,
    c: f64,
    x0: f64,
}

/// VWN3 未极化（paramagnetic，ζ = 0）参数，取自 libxc 的 `LDA_C_VWN_3`
const VWN3_P: VwnParams = VwnParams {
    a: 0.0310907,
    b: 13.0720,
    c: 42.7198,
    x0: -0.409286,
};

/// VWN3 完全极化（ferromagnetic，ζ = 1）参数
const VWN3_F: VwnParams = VwnParams {
    a: 0.01554535,
    b: 20.1231,
    c: 101.578,
    x0: -0.743294,
};

/// 由 `rs` 计算 VWN 的相关能 `eps` 及其对 `rs` 的导数
fn vwn_eps(rs: f64, p: &VwnParams) -> (f64, f64) {
    let q = (4.0 * p.c - p.b * p.b).sqrt();
    // x = sqrt(rs)
    let x = rs.sqrt();
    let xx = x * x + p.b * x + p.c; // X(x)
    let dx = 2.0 * x + p.b;
    let at = (q / dx).atan();
    let ln1 = (x * x / xx).ln();
    let ln2 = ((x - p.x0) * (x - p.x0) / xx).ln();
    let bx0 = p.b * p.x0 / (p.x0 * p.x0 + p.b * p.x0 + p.c);
    let b2x0 = p.b + 2.0 * p.x0;
    let eps = p.a * (ln1 + 2.0 * p.b / q * at - bx0 * (ln2 + 2.0 * b2x0 / q * at));

    // d eps / dx
    let dat = -2.0 * q / (dx * dx) / (1.0 + (q / dx) * (q / dx));
    let deps_dx =
        p.a * (2.0 / x - dx / xx + 2.0 * p.b / q * dat - bx0 * (2.0 / (x - p.x0) - dx / xx + 2.0 * b2x0 / q * dat));
    let deps_drs = deps_dx / (2.0 * x);
    (eps, deps_drs)
}

/// VWN3（Vosko-Wilk-Nusair 参量化 III）的未极化 LDA 相关能
///
/// Gaussian 的 B3LYP 使用 VWN3。外部 `rsxcs` 的 `lda_c_vwn_3` 与这里的
/// 参数化不一致（其未极化极限退化为 VWN5），因此两种后端都使用本实现，
/// 以保证 B3LYP 能量不随特性开关变化。
///
/// `zk` 返回能量密度（单位体积），`vrho` 返回其对密度的导数。
fn vwn3(rho: &[f64], zk: &mut [f64], vrho: &mut [f64]) {
    for i in 0..rho.len() {
        let r = rho[i];
        zk[i] = 0.0;
        vrho[i] = 0.0;
        if r <= 1.0e-20 {
            continue;
        }
        let rs = (3.0 / (4.0 * PI * r)).powf(1.0 / 3.0);
        let (eps, deps_drs) = vwn_eps(rs, &VWN3_P);
        zk[i] = r * eps;
        // vrho = d(ρ·eps)/dρ = eps - (rs/3)·d eps/d rs
        vrho[i] = eps - rs / 3.0 * deps_drs;
    }
}

/// VWN3 的自旋极化 LDA 相关能
///
/// 采用自旋插值 `ε_c(rs, ζ) = ε_P + (ε_F − ε_P)·f(ζ)`，其中
/// `f(ζ) = [(1+ζ)^(4/3) + (1−ζ)^(4/3) − 2] / (2^(4/3) − 2)`。
/// 两种后端都使用本实现，原因见 [`vwn3`]。
///
/// `zk` 返回能量密度（单位体积），`vrho_a`/`vrho_b` 返回其对 α/β 密度的导数。
fn vwn3_spin(rho_a: &[f64], rho_b: &[f64], zk: &mut [f64], vrho_a: &mut [f64], vrho_b: &mut [f64]) {
    let denom = 2.0_f64.powf(4.0 / 3.0) - 2.0;
    for i in 0..rho_a.len() {
        let ra = rho_a[i].max(0.0);
        let rb = rho_b[i].max(0.0);
        let rho = ra + rb;
        zk[i] = 0.0;
        vrho_a[i] = 0.0;
        vrho_b[i] = 0.0;
        if rho <= 1.0e-20 {
            continue;
        }
        let rs = (3.0 / (4.0 * PI * rho)).powf(1.0 / 3.0);
        let (eps_p, deps_p) = vwn_eps(rs, &VWN3_P);
        let (eps_f, deps_f) = vwn_eps(rs, &VWN3_F);

        let zeta = (ra - rb) / rho;
        let pa = 1.0 + zeta;
        let pb = 1.0 - zeta;
        let f = (pa.powf(4.0 / 3.0) + pb.powf(4.0 / 3.0) - 2.0) / denom;
        // df/dζ
        let df = 4.0 / 3.0 * (pa.powf(1.0 / 3.0) - pb.powf(1.0 / 3.0)) / denom;

        let eps = eps_p + (eps_f - eps_p) * f;
        let deps_drs = deps_p + (deps_f - deps_p) * f;
        let deps_dzeta = (eps_f - eps_p) * df;

        zk[i] = rho * eps;
        vrho_a[i] = eps - rs / 3.0 * deps_drs + (1.0 - zeta) * deps_dzeta;
        vrho_b[i] = eps - rs / 3.0 * deps_drs - (1.0 + zeta) * deps_dzeta;
    }
}

/// 内置 VWN3 的闭壳层求值结果（两种后端共用）
pub(super) fn vwn3_rks(rho: &[f64]) -> RksEval {
    let n = rho.len();
    let mut zk = vec![0.0; n];
    let mut vrho = vec![0.0; n];
    vwn3(rho, &mut zk, &mut vrho);
    for i in 0..n {
        zk[i] = if rho[i] > 1.0e-20 { zk[i] / rho[i] } else { 0.0 };
    }
    RksEval {
        exc: zk,
        vrho,
        vsigma: vec![0.0; n],
    }
}

/// 内置 VWN3 的开壳层求值结果（两种后端共用）
pub(super) fn vwn3_uks(rho_a: &[f64], rho_b: &[f64]) -> UksEval {
    let n = rho_a.len();
    let mut zk = vec![0.0; n];
    let mut vrho_a = vec![0.0; n];
    let mut vrho_b = vec![0.0; n];
    vwn3_spin(rho_a, rho_b, &mut zk, &mut vrho_a, &mut vrho_b);
    for i in 0..n {
        let rho = rho_a[i] + rho_b[i];
        zk[i] = if rho > 1.0e-20 { zk[i] / rho } else { 0.0 };
    }
    UksEval {
        exc: zk,
        vrho_a,
        vrho_b,
        vsigma_aa: vec![0.0; n],
        vsigma_bb: vec![0.0; n],
        vsigma_ab: vec![0.0; n],
    }
}

// ---------------------------------------------------------------------------
// DFTxclib 后端（仅未启用 use-rsxcs 时编译）
// ---------------------------------------------------------------------------

#[cfg(not(feature = "use-rsxcs"))]
mod backend {
    use super::super::{RksEval, UksEval, XcComponent};
    use super::{DFTxclib, vwn3_rks, vwn3_uks};

    /// 内置 `DFTxclib` 中 restricted (闭壳层) 泛函共用的签名
    type RksFn = fn(i32, usize, &[f64], &[f64], &mut [f64], &mut [f64], &mut [f64], &mut [f64], &mut [f64], &mut [f64]);

    /// `DFTxclib` 中 unrestricted (自旋极化) 泛函共用的签名
    ///
    /// 输入 5 个切片（ρ_α、ρ_β、σ_αα、σ_ββ、σ_αβ），
    /// 输出 6 个一阶导数切片，其余 15 个二阶导数切片在本仓库里不使用（传空切片）。
    #[allow(clippy::type_complexity)]
    type UksFn = fn(
        i32,
        usize,
        &[f64],
        &[f64],
        &[f64],
        &[f64],
        &[f64],
        &mut [f64],
        &mut [f64],
        &mut [f64],
        &mut [f64],
        &mut [f64],
        &mut [f64],
        &mut [f64],
        &mut [f64],
        &mut [f64],
        &mut [f64],
        &mut [f64],
        &mut [f64],
        &mut [f64],
        &mut [f64],
        &mut [f64],
        &mut [f64],
        &mut [f64],
        &mut [f64],
        &mut [f64],
        &mut [f64],
        &mut [f64],
    );

    fn rks_func(comp: XcComponent) -> Option<RksFn> {
        Some(match comp {
            XcComponent::LdaX => DFTxclib::rks_x_lda,
            XcComponent::Vwn5 => DFTxclib::rks_c_vwn5,
            XcComponent::B88 => DFTxclib::rks_x_b88,
            XcComponent::Lyp => DFTxclib::rks_c_lyp,
            XcComponent::P86 => DFTxclib::rks_c_p86,
            XcComponent::PbeX => DFTxclib::rks_x_pbe,
            XcComponent::PbeC => DFTxclib::rks_c_pbe,
            XcComponent::Vwn3 => return None,
        })
    }

    fn uks_func(comp: XcComponent) -> Option<UksFn> {
        Some(match comp {
            XcComponent::LdaX => DFTxclib::uks_x_lda,
            XcComponent::Vwn5 => DFTxclib::uks_c_vwn5,
            XcComponent::B88 => DFTxclib::uks_x_b88,
            XcComponent::Lyp => DFTxclib::uks_c_lyp,
            XcComponent::P86 => DFTxclib::uks_c_p86,
            XcComponent::PbeX => DFTxclib::uks_x_pbe,
            XcComponent::PbeC => DFTxclib::uks_c_pbe,
            XcComponent::Vwn3 => return None,
        })
    }

    pub fn eval_rks(comp: XcComponent, rho: &[f64], sigma: &[f64]) -> RksEval {
        if comp == XcComponent::Vwn3 {
            return vwn3_rks(rho);
        }
        let n = rho.len();
        let mut zk = vec![0.0; n];
        let mut vrho = vec![0.0; n];
        let mut vsigma = vec![0.0; n];
        let func = rks_func(comp).expect("未知的泛函分量");
        let mut vsigma_aa = vec![0.0; n];
        func(
            1,
            n,
            rho,
            sigma,
            &mut zk,
            &mut vrho,
            &mut vsigma_aa,
            &mut [],
            &mut [],
            &mut [],
        );
        // DFTxclib 的 restricted vsigma 对应 sigma_aa = |∇ρ|²/4
        for i in 0..n {
            vsigma[i] = 0.25 * vsigma_aa[i];
        }
        for i in 0..n {
            zk[i] = if rho[i] > 1.0e-20 { zk[i] / rho[i] } else { 0.0 };
        }
        RksEval { exc: zk, vrho, vsigma }
    }

    pub fn eval_uks(comp: XcComponent, rho_a: &[f64], rho_b: &[f64], saa: &[f64], sbb: &[f64], sab: &[f64]) -> UksEval {
        if comp == XcComponent::Vwn3 {
            return vwn3_uks(rho_a, rho_b);
        }
        let n = rho_a.len();
        let mut zk = vec![0.0; n];
        let mut vrho_a = vec![0.0; n];
        let mut vrho_b = vec![0.0; n];
        let mut vsigma_aa = vec![0.0; n];
        let mut vsigma_bb = vec![0.0; n];
        let mut vsigma_ab = vec![0.0; n];
        let func = uks_func(comp).expect("未知的泛函分量");
        func(
            1,
            n,
            rho_a,
            rho_b,
            saa,
            sbb,
            sab,
            &mut zk,
            &mut vrho_a,
            &mut vrho_b,
            &mut vsigma_aa,
            &mut vsigma_bb,
            &mut vsigma_ab,
            &mut [],
            &mut [],
            &mut [],
            &mut [],
            &mut [],
            &mut [],
            &mut [],
            &mut [],
            &mut [],
            &mut [],
            &mut [],
            &mut [],
            &mut [],
            &mut [],
            &mut [],
        );
        for i in 0..n {
            let rho = rho_a[i] + rho_b[i];
            zk[i] = if rho > 1.0e-20 { zk[i] / rho } else { 0.0 };
        }
        UksEval {
            exc: zk,
            vrho_a,
            vrho_b,
            vsigma_aa,
            vsigma_bb,
            vsigma_ab,
        }
    }
}

#[cfg(not(feature = "use-rsxcs"))]
pub use backend::{eval_rks, eval_uks};
