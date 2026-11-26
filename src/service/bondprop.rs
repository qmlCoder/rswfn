//! 键属性分析任务。

use ndarray::Array2;

use crate::base::Mole;
use crate::props::bondprop;
use crate::report::{Column, Report, Table};

use super::{ServiceResult, check_atom, report_for, require_coefs};

/// 键级表：只保留 |键级| > 0.01 的原子对
fn bond_order_table(mole: &Mole, omat: &Array2<f64>) -> Table {
    let syms = mole.atoms().syms();
    let natm = mole.atoms().len();
    let mut table = Table::new(vec![Column::left("原子对"), Column::right("序号"), Column::right("键级")]);
    for i in 0..natm {
        for j in (i + 1)..natm {
            let order = omat[(i, j)];
            if order.abs() > 0.01 {
                let pair = format!("{}-{}", syms[i], syms[j]);
                table.push_row(vec![pair, format!("{}-{}", i + 1, j + 1), format!("{order:.6}")]);
            }
        }
    }
    table
}

fn bond_order_report(mole: &Mole, omat: &Array2<f64>, title: &str) -> Report {
    let mut report = report_for(title, mole);
    report.tables.push(bond_order_table(mole, omat));
    report.notes.push("仅显示 |键级| > 0.01 的原子对".to_string());
    report
}

/// 校验原子对（1 基，且不能是同一个原子）
fn check_pair(mole: &Mole, atm1: usize, atm2: usize) -> ServiceResult<()> {
    check_atom(mole, atm1)?;
    check_atom(mole, atm2)?;
    if atm1 == atm2 {
        return Err("需要两个不同的原子".to_string());
    }
    Ok(())
}

/// Mayer 键级
pub fn mayer(mole: &Mole) -> ServiceResult<Report> {
    require_coefs(mole)?;
    let caler = bondprop::order::Calculator::new(mole);
    let omat = caler.mayer();
    Ok(bond_order_report(mole, &omat, "Mayer 键级分析"))
}

/// Wiberg 键级
pub fn wiberg(mole: &Mole) -> ServiceResult<Report> {
    require_coefs(mole)?;
    let caler = bondprop::order::Calculator::new(mole);
    let omat = caler.wiberg();
    Ok(bond_order_report(mole, &omat, "Wiberg 键级分析"))
}

/// 指定分子轨道对键级的贡献
pub fn order_obt(mole: &Mole, obt: usize, otype: &str) -> ServiceResult<Report> {
    require_coefs(mole)?;
    let nobt = mole.coefs().obt_occs.len();
    if obt == 0 || obt > nobt {
        return Err(format!("轨道序号越界: {obt}（共 {nobt} 个轨道，序号从 1 开始）"));
    }
    let caler = bondprop::order::Calculator::new(mole);
    let omat = caler.order_obt(otype, obt - 1);
    Ok(bond_order_report(mole, &omat, &format!("分子轨道 {obt} 对{otype}键级的贡献")))
}

/// π 键级（POCV / MOCV）
pub fn pi_order(mole: &Mole, method: &str, otype: &str) -> ServiceResult<Report> {
    require_coefs(mole)?;
    let caler = bondprop::order::Calculator::new(mole);
    let omat = match method {
        "pocv" => caler.pi_pocv(otype)?.1,
        "mocv" => caler.pi_mocv(otype).2,
        other => return Err(format!("未知的轨道改造方法: {other}（可用: pocv | mocv）")),
    };
    Ok(bond_order_report(mole, &omat, &format!("π键级分析（{}，{otype}）", method.to_uppercase())))
}

/// 键级分解为 σ、π、δ、φ
pub fn decom(mole: &Mole, atm1: usize, atm2: usize) -> ServiceResult<Report> {
    require_coefs(mole)?;
    check_pair(mole, atm1, atm2)?;
    let caler = bondprop::order::Calculator::new(mole);
    let vals = caler.decom(atm1 - 1, atm2 - 1, "mayer");
    let syms = mole.atoms().syms();

    let mut report = report_for(format!("Mayer 键级分解: 原子{}({}) - 原子{}({})", atm1, syms[atm1 - 1], atm2, syms[atm2 - 1]), mole);
    let mut table = Table::new(vec![Column::left("分量"), Column::right("键级")]);
    for (name, value) in ["σ键级", "π键级", "δ键级", "φ键级"].iter().zip(vals) {
        table.push_row(vec![name.to_string(), format!("{value:.6}")]);
    }
    report.tables.push(table);
    Ok(report)
}

/// 键能量分解（bond / rest / cross）
pub fn bond_energy(mole: &Mole, atm1: usize, atm2: usize) -> ServiceResult<Report> {
    require_coefs(mole)?;
    check_pair(mole, atm1, atm2)?;
    let caler = bondprop::energy::Calculator::new(mole);
    let e = caler.bond_parts(atm1 - 1, atm2 - 1)?;
    let syms = mole.atoms().syms();

    let mut report = report_for(format!("键能量分解: 原子{}({}) - 原子{}({})（单位: Hartree）", atm1, syms[atm1 - 1], atm2, syms[atm2 - 1]), mole);
    let mut table = Table::new(vec![Column::left("项目"), Column::left("表达式"), Column::right("数值")]);
    let rows = [
        ("键自能量", "E_bond", e.e_bond),
        ("其余部分", "E_rest", e.e_rest),
        ("键-其余交叉项", "E_cross", e.e_cross),
        ("  · 库仑", "Tr(D_bond J_rest)", e.e_cross_coul),
        ("  · 精确交换", "-a_x/2 Tr(D_bond K_rest)", e.e_cross_exch),
        ("  · XC 残差", "ΔE_xc", e.e_cross_xc),
        ("完整电子能量", "E_elec", e.e_elec),
        ("核排斥能", "E_nuc", e.e_nuc),
        ("总能量", "E_elec + E_nuc", e.e_total),
        ("键可加贡献", "E_bond + E_cross", e.bond()),
    ];
    for (label, expr, value) in rows {
        table.push_row(vec![label.to_string(), expr.to_string(), format!("{value:.10}")]);
    }
    report.tables.push(table);
    Ok(report)
}
