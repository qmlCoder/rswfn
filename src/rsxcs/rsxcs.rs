//! 外部 `rsxcs` 0.2 后端：直接调用 libxc 7.1.2 子集的 Rust 实现。
//!
//! 本模块负责两件事：
//! - 实现统一接口 [`eval_rks`] / [`eval_uks`]，把项目内的 [`XcComponent`] 映射到
//!   `rsxcs` 的泛函枚举（其中 VWN3 仍用内置实现，见 [`crate::rsxcs::hand`]）；
//! - 支持直接求值外部库提供的**任意**泛函（[`run_rks`] / [`run_uks`]），供
//!   [`crate::rsxcs::XcFunc::External`] 使用。

use super::hand::{vwn3_rks, vwn3_uks};
use super::{RksEval, UksEval, XcComponent};
use ::rsxcs::{Family, Functional, Order, Spin, Xc};

/// 把本项目使用的分量映射到 rsxcs 的泛函枚举（VWN3 走内部实现）
fn functional(comp: XcComponent) -> Functional {
    match comp {
        XcComponent::LdaX => Functional::LdaX,
        XcComponent::Vwn5 => Functional::LdaCVwn,
        XcComponent::Vwn3 => unreachable!("VWN3 使用内置实现"),
        XcComponent::B88 => Functional::GgaXB88,
        XcComponent::Lyp => Functional::GgaCLyp,
        XcComponent::P86 => Functional::GgaCP86,
        XcComponent::PbeX => Functional::GgaXPbe,
        XcComponent::PbeC => Functional::GgaCPbe,
    }
}

fn eval(comp: XcComponent, spin: Spin, rho: &[f64], sigma: &[f64]) -> ::rsxcs::Evaluation {
    let xc = Xc::from_functional(functional(comp));
    let sigma = if comp.is_gga() { Some(sigma) } else { None };
    xc.eval(spin, rho, sigma, None, None, Order::Vxc)
        .expect("rsxcs 泛函求值失败")
}

pub fn eval_rks(comp: XcComponent, rho: &[f64], sigma: &[f64]) -> RksEval {
    if comp == XcComponent::Vwn3 {
        return vwn3_rks(rho);
    }
    let out = eval(comp, Spin::Unpolarized, rho, sigma);
    let vsigma = if comp.is_gga() {
        out.vsigma
    } else {
        vec![0.0; rho.len()]
    };
    RksEval {
        exc: out.exc,
        vrho: out.vxc,
        vsigma,
    }
}

pub fn eval_uks(comp: XcComponent, rho_a: &[f64], rho_b: &[f64], saa: &[f64], sbb: &[f64], sab: &[f64]) -> UksEval {
    if comp == XcComponent::Vwn3 {
        return vwn3_uks(rho_a, rho_b);
    }
    let n = rho_a.len();
    let mut rho = Vec::with_capacity(2 * n);
    for i in 0..n {
        rho.push(rho_a[i]);
        rho.push(rho_b[i]);
    }
    let mut sigma = Vec::new();
    if comp.is_gga() {
        sigma.reserve(3 * n);
        for i in 0..n {
            sigma.push(saa[i]);
            sigma.push(sab[i]);
            sigma.push(sbb[i]);
        }
    }
    let out = eval(comp, Spin::Polarized, &rho, &sigma);
    let mut vrho_a = vec![0.0; n];
    let mut vrho_b = vec![0.0; n];
    let mut vsigma_aa = vec![0.0; n];
    let mut vsigma_bb = vec![0.0; n];
    let mut vsigma_ab = vec![0.0; n];
    for i in 0..n {
        vrho_a[i] = out.vxc[2 * i];
        vrho_b[i] = out.vxc[2 * i + 1];
        if comp.is_gga() {
            vsigma_aa[i] = out.vsigma[3 * i];
            vsigma_ab[i] = out.vsigma[3 * i + 1];
            vsigma_bb[i] = out.vsigma[3 * i + 2];
        }
    }
    UksEval {
        exc: out.exc,
        vrho_a,
        vrho_b,
        vsigma_aa,
        vsigma_bb,
        vsigma_ab,
    }
}

// ---------------------------------------------------------------------------
// 整份外部泛函求值（供 `XcFunc::External` 使用，支持 rsxcs 的全部泛函）
// ---------------------------------------------------------------------------

/// LDA 泛函不需要 sigma；其余（GGA/mGGA）需要。
fn needs_sigma(func: Functional) -> bool {
    !matches!(func.family().semilocal(), Family::Lda)
}

pub(super) fn run_rks(id: u32, rho: &[f64], sigma: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let func = Functional::from_id(id).unwrap_or_else(|| panic!("未知的 libxc 泛函编号 {id}"));
    let sigma = if needs_sigma(func) { Some(sigma) } else { None };
    let out = Xc::from_functional(func)
        .eval(Spin::Unpolarized, rho, sigma, None, None, Order::Vxc)
        .unwrap_or_else(|e| panic!("rsxcs 泛函 {} 求值失败: {e}", func.info().name));
    let vsigma = if out.vsigma.is_empty() {
        vec![0.0; rho.len()]
    } else {
        out.vsigma
    };
    (out.exc, out.vxc, vsigma)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn run_uks(
    id: u32,
    rho_a: &[f64],
    rho_b: &[f64],
    saa: &[f64],
    sbb: &[f64],
    sab: &[f64],
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let func = Functional::from_id(id).unwrap_or_else(|| panic!("未知的 libxc 泛函编号 {id}"));
    let n = rho_a.len();
    let mut rho = Vec::with_capacity(2 * n);
    for i in 0..n {
        rho.push(rho_a[i]);
        rho.push(rho_b[i]);
    }
    let mut sigma = Vec::new();
    if needs_sigma(func) {
        sigma.reserve(3 * n);
        for i in 0..n {
            sigma.push(saa[i]);
            sigma.push(sab[i]);
            sigma.push(sbb[i]);
        }
    }
    let sigma = if sigma.is_empty() { None } else { Some(sigma.as_slice()) };
    let out = Xc::from_functional(func)
        .eval(Spin::Polarized, &rho, sigma, None, None, Order::Vxc)
        .unwrap_or_else(|e| panic!("rsxcs 泛函 {} 求值失败: {e}", func.info().name));

    let mut vrho_a = vec![0.0; n];
    let mut vrho_b = vec![0.0; n];
    let mut vsigma_aa = vec![0.0; n];
    let mut vsigma_bb = vec![0.0; n];
    let mut vsigma_ab = vec![0.0; n];
    let gga = !out.vsigma.is_empty();
    for i in 0..n {
        vrho_a[i] = out.vxc[2 * i];
        vrho_b[i] = out.vxc[2 * i + 1];
        if gga {
            vsigma_aa[i] = out.vsigma[3 * i];
            vsigma_ab[i] = out.vsigma[3 * i + 1];
            vsigma_bb[i] = out.vsigma[3 * i + 2];
        }
    }
    (out.exc, vrho_a, vrho_b, vsigma_aa, vsigma_bb, vsigma_ab)
}
