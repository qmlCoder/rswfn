//! 分子属性分析任务。

use crate::base::Mole;
use crate::props::moleprop;
use crate::report::{Column, Field, Report, Table};

use super::{ServiceResult, report_for, require_coefs};

/// 解析形如 "1-3,5,7-9" 的原子列表（1 基、含逗号/分号/空格分隔与区间），返回 0 基索引
pub fn parse_atom_list(spec: &str, natm: usize) -> ServiceResult<Vec<usize>> {
    let mut idxs = vec![];
    for token in spec.split(|c: char| c == ',' || c == ';' || c.is_whitespace()) {
        let token = token.trim();
        if token.is_empty() {
            continue;
        }
        let nums: Vec<usize> = token.split('-').map(|s| s.trim().parse().map_err(|_| format!("无法解析原子序号: {spec}"))).collect::<Result<_, _>>()?;
        let (a, b) = match nums[..] {
            [a] => (a, a),
            [a, b] => (a, b),
            _ => return Err(format!("无法解析原子序号区间: {spec}")),
        };
        if a == 0 || b == 0 || a > b || b > natm {
            return Err(format!("原子序号越界: {spec}（分子共有 {natm} 个原子，序号从 1 开始）"));
        }
        idxs.extend((a..=b).map(|i| i - 1));
    }
    Ok(idxs)
}

/// 芳香性分析（PISD）
pub fn aromacity(mole: &Mole, rings: &[String]) -> ServiceResult<Report> {
    require_coefs(mole)?;
    let natm = mole.atoms().len();
    let caler = moleprop::aromacity::Calculator::new(mole);

    // 未指定环时自动检测分子中的最小环（SSSR）
    let rings: Vec<Vec<usize>> = if rings.is_empty() {
        mole.rings()
    } else {
        rings.iter().map(|ring| parse_atom_list(ring, natm)).collect::<Result<_, _>>()?
    };

    if rings.is_empty() {
        // 无环结构：回退为计算全分子非氢键的π键级和
        let vals = caler.PISD(None)?;
        let mut report = report_for("芳香性分析（π键级方法）", mole);
        report.fields.push(Field::new("全分子π键级和", format!("{:.6}", vals[0])));
        report.notes.push("未检测到环，按全分子统计".to_string());
        return Ok(report);
    }

    let vals = caler.PISD(Some(&rings))?;
    let syms = mole.atoms().syms();

    let mut report = report_for("芳香性分析（π键级方法）", mole);
    let mut table = Table::new(vec![Column::right("环"), Column::left("原子列表"), Column::right("PISD 值")]);
    for (i, (ring, &val)) in rings.iter().zip(vals.iter()).enumerate() {
        let desc: Vec<String> = ring.iter().map(|&a| format!("{}({})", a + 1, syms[a])).collect();
        table.push_row(vec![(i + 1).to_string(), desc.join(" "), format!("{val:.6}")]);
    }
    report.tables.push(table);
    Ok(report)
}

/// π 电子能量分解（POCV / MOCV）
pub fn pi_energy(mole: &Mole, method: &str) -> ServiceResult<Report> {
    require_coefs(mole)?;
    let caler = moleprop::energy::Calculator::new(mole);
    let e = match method {
        "pocv" => caler.pi_pocv_parts(),
        "mocv" => caler.pi_mocv_parts(),
        other => return Err(format!("未知的轨道改造方法: {other}（可用: pocv | mocv）")),
    }?;

    let mut report = report_for(format!("π 电子能量分解（{} 方法，单位: Hartree）", method.to_uppercase()), mole);
    let mut table = Table::new(vec![Column::left("项目"), Column::left("表达式"), Column::right("数值")]);
    let rows = [
        ("π 自能量", "E_π", e.e_pi),
        ("σ 自能量", "E_σ", e.e_sigma),
        ("π-σ 交叉项", "E_πσ", e.e_cross),
        ("  · 库仑", "Tr(D_π J_σ)", e.e_cross_coul),
        ("  · 精确交换", "-a_x/2 Tr(D_π K_σ)", e.e_cross_exch),
        ("  · XC 残差", "ΔE_xc", e.e_cross_xc),
        ("完整电子能量", "E_elec", e.e_elec),
        ("核排斥能", "E_nuc", e.e_nuc),
        ("总能量", "E_elec + E_nuc", e.e_total),
        ("π 可加贡献", "E_π + E_πσ", e.pi()),
    ];
    for (label, expr, value) in rows {
        table.push_row(vec![label.to_string(), expr.to_string(), format!("{value:.10}")]);
    }
    report.tables.push(table);
    Ok(report)
}
