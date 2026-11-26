//! 原子属性分析任务。

use std::collections::HashMap;

use crate::base::{Mole, Vec3};
use crate::props::atomprop;
use crate::props::atomprop::activity;
use crate::report::{Column, Field, Report, Table};

use super::{ServiceResult, report_for, require_coefs};

/// 原子布居表：原子 / 符号 / 数值
fn population_table(mole: &Mole, pops: &[f64], value_name: &str) -> Table {
    let syms = mole.atoms().syms();
    let mut table = Table::new(vec![Column::right("原子"), Column::left("符号"), Column::right(value_name)]);
    for (i, &pop) in pops.iter().enumerate() {
        let sym = syms.get(i).cloned().unwrap_or_default();
        table.push_row(vec![(i + 1).to_string(), sym, format!("{pop:.8}")]);
    }
    table
}

/// 原子布居结果（含总电子数）
pub(super) fn population_report(mole: &Mole, pops: &[f64], title: &str) -> Report {
    let mut report = report_for(title, mole);
    report.fields.push(Field::new("总电子数", format!("{:.8}", pops.iter().sum::<f64>())));
    report.tables.push(population_table(mole, pops, "电子布居"));
    report
}

/// 福井函数表（7 列）
fn fukui_table(mole: &Mole, vals: &[[f64; 7]]) -> Table {
    let syms = mole.atoms().syms();
    let mut table = Table::new(vec![
        Column::right("原子"),
        Column::left("符号"),
        Column::right("e0(N)"),
        Column::right("e(N+1)"),
        Column::right("e(N-1)"),
        Column::right("f-"),
        Column::right("f+"),
        Column::right("f0"),
        Column::right("Δf"),
    ]);
    for (i, v) in vals.iter().enumerate() {
        let sym = syms.get(i).cloned().unwrap_or_default();
        table.push_row(vec![
            (i + 1).to_string(),
            sym,
            format!("{:.8}", v[0]),
            format!("{:.8}", v[1]),
            format!("{:.8}", v[2]),
            format!("{:.6}", v[3]),
            format!("{:.6}", v[4]),
            format!("{:.6}", v[5]),
            format!("{:.6}", v[6]),
        ]);
    }
    table
}

/// π 方向向量表
fn dirs_table(mole: &Mole, dirs: &HashMap<usize, Vec3>) -> Table {
    let syms = mole.atoms().syms();
    let mut idxs: Vec<usize> = dirs.keys().copied().collect();
    idxs.sort_unstable();
    let mut table = Table::new(vec![Column::right("原子"), Column::left("符号"), Column::right("X"), Column::right("Y"), Column::right("Z")]);
    for idx in idxs {
        let [x, y, z] = dirs[&idx].into_arr();
        let sym = syms.get(idx).cloned().unwrap_or_default();
        table.push_row(vec![(idx + 1).to_string(), sym, format!("{x:.10}"), format!("{y:.10}"), format!("{z:.10}")]);
    }
    table
}

/// Mulliken 电子布居
pub fn mulliken(mole: &Mole) -> ServiceResult<Report> {
    require_coefs(mole)?;
    let caler = atomprop::charge::Calculator::new(mole);
    Ok(population_report(mole, &caler.mulliken(), "Mulliken 电子布居分析"))
}

/// Löwdin 电子布居
pub fn lowdin(mole: &Mole) -> ServiceResult<Report> {
    require_coefs(mole)?;
    let caler = atomprop::charge::Calculator::new(mole);
    Ok(population_report(mole, &caler.lowdin(), "Löwdin 电子布居分析"))
}

/// Hirshfeld 电子布居
pub fn hirshfeld(mole: &Mole) -> ServiceResult<Report> {
    require_coefs(mole)?;
    let caler = atomprop::charge::Calculator::new(mole);
    Ok(population_report(mole, &caler.hirshfeld(), "Hirshfeld 电子布居分析"))
}

/// 自旋密度布居
pub fn spin(mole: &Mole, ctype: &str) -> ServiceResult<Report> {
    require_coefs(mole)?;
    if !mole.is_open() {
        return Err("非开壳层分子（多重度=1）无法计算自旋布居".to_string());
    }
    let caler = atomprop::charge::Calculator::new(mole);
    let spin = caler.spin(ctype);
    let syms = mole.atoms().syms();

    let mut report = report_for(format!("自旋密度布居分析（{ctype}）"), mole);
    let mut table = Table::new(vec![Column::right("原子"), Column::left("符号"), Column::right("α电子"), Column::right("β电子"), Column::right("α - β")]);
    for (i, &[a, b, d]) in spin.iter().enumerate() {
        let sym = syms.get(i).cloned().unwrap_or_default();
        table.push_row(vec![(i + 1).to_string(), sym, format!("{a:.8}"), format!("{b:.8}"), format!("{d:.8}")]);
    }
    report.tables.push(table);
    Ok(report)
}

/// RESP 电荷
pub fn resp(mole: &Mole, pt_density: f64, one_stage: bool) -> ServiceResult<Report> {
    require_coefs(mole)?;
    let caler = atomprop::resp::Calculator::new(mole);
    let mut opts = atomprop::resp::RespOptions::default();
    opts.pt_density = pt_density;
    opts.two_stage = !one_stage;
    let res = caler.resp_with(&opts)?;

    let title = if res.stage2 { "RESP 电荷分析（标准两阶段）" } else { "RESP 电荷分析（第一阶段）" };
    let mut report = report_for(title, mole);
    let syms = mole.atoms().syms();
    report.fields.push(Field::new("MK 拟合点数", res.n_points.to_string()));
    report.fields.push(Field::new("总电荷", format!("{:.8}", res.charges.iter().sum::<f64>())));
    report.fields.push(Field::new("ESP RMSE", format!("{:.6e}", res.rmse)));
    report.fields.push(Field::new("ESP RRMSE", format!("{:.6e}", res.rrmse)));
    let iterations = if res.stage2 {
        format!("阶段1 {} 次，阶段2 {} 次", res.stage1_iter, res.stage2_iter)
    } else {
        format!("阶段1 {} 次", res.stage1_iter)
    };
    report.fields.push(Field::new("迭代", iterations));

    let mut table = if res.stage2 {
        Table::new(vec![Column::right("原子"), Column::left("符号"), Column::right("阶段1电荷"), Column::right("RESP电荷")])
    } else {
        Table::new(vec![Column::right("原子"), Column::left("符号"), Column::right("RESP电荷")])
    };
    for (i, &charge) in res.charges.iter().enumerate() {
        let sym = syms.get(i).cloned().unwrap_or_default();
        if res.stage2 {
            table.push_row(vec![(i + 1).to_string(), sym, format!("{:.8}", res.stage1[i]), format!("{charge:.8}")]);
        } else {
            table.push_row(vec![(i + 1).to_string(), sym, format!("{charge:.8}")]);
        }
    }
    report.tables.push(table);

    if !res.converged {
        report.notes.push("警告: RESP 拟合未收敛".to_string());
    }
    Ok(report)
}

/// 原子 π 电子布居（POCV / MOCV）
pub fn pi_pop(mole: &Mole, ctype: &str, method: &str) -> ServiceResult<Report> {
    require_coefs(mole)?;
    let caler = atomprop::charge::Calculator::new(mole);
    let eles = match method {
        "pocv" => caler.pi_pocv(ctype)?.1,
        "mocv" => caler.pi_mocv(ctype)?.2,
        other => return Err(format!("未知的轨道改造方法: {other}（可用: pocv | mocv）")),
    };
    Ok(population_report(mole, &eles, &format!("原子π电子布居（{}，{ctype}）", method.to_uppercase())))
}

/// 福井函数与双描述符
pub fn fukui(mole: &Mole, mole_n: &Mole, mole_p: &Mole, ctype: &str) -> ServiceResult<Report> {
    require_coefs(mole)?;
    require_coefs(mole_n)?;
    require_coefs(mole_p)?;
    let caler = activity::Calculator::new(mole);
    let vals = caler.fukui(mole_n, mole_p, ctype);

    let mut report = report_for(format!("福井函数与双描述符分析（{ctype}）"), mole);
    report.tables.push(fukui_table(mole, &vals));
    Ok(report)
}

/// π 轨道的福井函数与双描述符
pub fn fukui_pi(mole: &Mole, mole_n: &Mole, mole_p: &Mole, ctype: &str, method: &str) -> ServiceResult<Report> {
    require_coefs(mole)?;
    require_coefs(mole_n)?;
    require_coefs(mole_p)?;
    let caler = activity::Calculator::new(mole);

    let (mut report, dirs, vals) = match method {
        "pocv" => {
            let (dirs, vals) = caler.fukui_pi_pocv(mole_n, mole_p, ctype)?;
            let report = report_for(format!("π福井函数与双描述符分析（POCV，{ctype}）"), mole);
            (report, Some(dirs), vals)
        }
        "mocv" => {
            let (_stms, _atos, vals) = caler.fukui_pi_mocv(mole_n, mole_p, ctype)?;
            let report = report_for(format!("π福井函数与双描述符分析（MOCV，{ctype}）"), mole);
            (report, None, vals)
        }
        other => return Err(format!("未知的轨道改造方法: {other}（可用: pocv | mocv）")),
    };

    report.tables.push(fukui_table(mole, &vals));
    if let Some(dirs) = dirs {
        report.tables.push(dirs_table(mole, &dirs));
    }
    Ok(report)
}

/// 自由价
pub fn freev(mole: &Mole, otype: &str) -> ServiceResult<Report> {
    require_coefs(mole)?;
    let caler = activity::Calculator::new(mole);
    let vals = caler.freev(otype);
    Ok(freev_report(mole, &vals, &format!("自由价分析（{otype}）")))
}

/// π 自由价
pub fn freev_pi(mole: &Mole, otype: &str) -> ServiceResult<Report> {
    require_coefs(mole)?;
    let caler = activity::Calculator::new(mole);
    let (_dirs, vals) = caler.freev_pi(otype)?;
    Ok(freev_report(mole, &vals, &format!("π自由价分析（{otype}）")))
}

fn freev_report(mole: &Mole, vals: &[[f64; 3]], title: &str) -> Report {
    let syms = mole.atoms().syms();
    let mut report = report_for(title, mole);
    let mut table = Table::new(vec![Column::right("原子"), Column::left("符号"), Column::right("标准价"), Column::right("键级和"), Column::right("自由价")]);
    for (i, v) in vals.iter().enumerate() {
        let sym = syms.get(i).cloned().unwrap_or_default();
        table.push_row(vec![(i + 1).to_string(), sym, format!("{:.6}", v[0]), format!("{:.6}", v[1]), format!("{:.6}", v[2])]);
    }
    report.tables.push(table);
    report
}
