//! 轨道属性分析任务。

use crate::base::Mole;
use crate::props::atomprop;
use crate::props::orbtprop;
use crate::report::{Column, Report, Table};

use super::{ServiceResult, require_coefs};

/// 每个分子轨道的 π 电子数表
fn pi_eles_table(pi_eles: &[f64]) -> Table {
    let mut table = Table::new(vec![Column::right("轨道"), Column::right("π 电子数")]);
    for (i, &val) in pi_eles.iter().enumerate() {
        table.push_row(vec![(i + 1).to_string(), format!("{val:.8}")]);
    }
    table
}

/// 每个分子轨道的 π 电子数（POCV 方法）
pub fn pi_eles_pocv(mole: &Mole) -> ServiceResult<Report> {
    require_coefs(mole)?;
    let caler = orbtprop::obtval::Calculator::new(mole);
    let pi_eles = caler.pi_eles_pocv();
    let mut report = super::report_for("分子轨道 π 电子数（POCV 方法）", mole);
    report.tables.push(pi_eles_table(&pi_eles));
    Ok(report)
}

/// 每个分子轨道的 π 电子数（MOCV 方法）
pub fn pi_eles_mocv(mole: &Mole) -> ServiceResult<Report> {
    require_coefs(mole)?;
    let caler = orbtprop::obtval::Calculator::new(mole);
    let pi_eles = caler.pi_eles_mocv();
    let mut report = super::report_for("分子轨道 π 电子数（MOCV 方法）", mole);
    report.tables.push(pi_eles_table(&pi_eles));
    Ok(report)
}

/// 指定分子轨道的原子电子布居
pub fn composition(mole: &Mole, obt: usize, ctype: &str) -> ServiceResult<Report> {
    require_coefs(mole)?;
    let nobt = mole.coefs().obt_occs.len();
    if obt == 0 || obt > nobt {
        return Err(format!("轨道序号越界: {obt}（共 {nobt} 个轨道，序号从 1 开始）"));
    }
    let caler = atomprop::charge::Calculator::new(mole);
    let pops = caler.eles_obt(ctype, obt - 1);
    Ok(super::atomprop::population_report(mole, &pops, &format!("第 {obt} 个分子轨道的 {ctype} 电子布居")))
}
