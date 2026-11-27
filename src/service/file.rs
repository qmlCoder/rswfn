//! 文件级任务：坐标提取、格式转换与支持信息。

use std::path::{Path, PathBuf, absolute};

use crate::base::Mole;
use crate::reader::{GjfReader, LogReader};
use crate::report::{Column, Field, Report, Table};
use crate::tools::log::Tool;
use crate::utils;
use crate::writer::{DmoWriter, GjfWriter, Writer, XyzWriter};

use super::{ServiceResult, absolute_str, ext_of};

/// 提取分子坐标，可导出 .xyz
pub fn coord(input: &Path, output: Option<&Path>, unit: &str) -> ServiceResult<Report> {
    let mole = super::load_any_mole(input)?;
    let atoms = mole.atoms();
    let natm = atoms.len();
    let (xyzs, unit_name) = match unit {
        "bohr" => (atoms.xyzs().clone(), "Bohr"),
        "ang" => (atoms.xyzs_ang(), "Å"),
        other => return Err(format!("未知的坐标单位: {other}（可用: ang | bohr）")),
    };

    let mut report = Report::new("分子坐标");
    report.fields = vec![
        Field::new("文件", mole.path.clone().unwrap_or_default()),
        Field::new("分子式", mole.formula()),
        Field::new("电荷", mole.charge().to_string()),
        Field::new("多重度", mole.multip().to_string()),
        Field::new("原子数", natm.to_string()),
    ];

    let syms = atoms.syms();
    let mut table = Table::new(vec![
        Column::right("序号"),
        Column::left("符号"),
        Column::right(format!("X ({unit_name})")),
        Column::right(format!("Y ({unit_name})")),
        Column::right(format!("Z ({unit_name})")),
    ]);
    for i in 0..natm {
        let [x, y, z] = xyzs[i];
        table.push_row(vec![(i + 1).to_string(), syms.get(i).cloned().unwrap_or_default(), format!("{x:.10}"), format!("{y:.10}"), format!("{z:.10}")]);
    }
    report.tables.push(table);

    // 可选：导出到 .xyz 文件（xyz 格式固定使用埃）
    if let Some(output) = output {
        let dest = absolute(output).map_err(|err| format!("无法解析输出路径 {}: {err}", output.display()))?;
        if ext_of(&dest).as_deref() != Some("xyz") {
            return Err(format!("坐标导出目前仅支持 .xyz 文件格式，收到: {}", dest.display()));
        }
        let dest_path = dest.to_str().ok_or_else(|| format!("输出路径不是有效的 UTF-8: {}", dest.display()))?;
        let mut writer = XyzWriter::new();
        writer.title = format!("{} extracted by rswfn", mole.formula());
        writer.read_mole(&mole);
        writer.save(dest_path);
        report.notes.push(format!("坐标已导出到: {dest_path}"));
    }
    Ok(report)
}

/// 文件格式转换（log/out → gjf，gjf → dmi）
pub fn trans(input: &Path, output: &Path) -> ServiceResult<Report> {
    let from_path = absolute_str(input, "输入文件")?;
    let dest_path = absolute_str(output, "输出文件")?;
    let from_type = ext_of(Path::new(&from_path)).ok_or_else(|| format!("输入文件缺少扩展名: {from_path}"))?;
    let dest_type = ext_of(Path::new(&dest_path)).ok_or_else(|| format!("输出文件缺少扩展名: {dest_path}"))?;

    match (from_type.as_str(), dest_type.as_str()) {
        ("log" | "out", "gjf") => {
            let reader = LogReader::new(&from_path);
            let mole = Mole::from_read(&reader);
            let job = reader.read_job();
            let mut writer = GjfWriter::new();
            writer.read_mole(&mole);
            writer.job = job.to_owned();
            writer.chk = utils::file_name(&dest_path);
            writer.charge = reader.charge as isize;
            writer.multip = reader.multip as usize;
            writer.save(&dest_path);
        }
        ("gjf", "dmi") => {
            let reader = GjfReader::new(&from_path);
            let mole = Mole::from_read(&reader);
            let mut writer = DmoWriter::blank();
            writer.read_mole(&mole);
            writer.save(&dest_path);
        }
        _ => {
            return Err(format!("不支持的转换类型: {from_type} -> {dest_type}（支持 log/out→gjf、gjf→dmi）"));
        }
    }

    let mut report = Report::new("文件格式转换");
    report.fields.push(Field::new("转换", format!("{from_type} -> {dest_type}")));
    report.fields.push(Field::new("输入", from_path));
    report.fields.push(Field::new("输出", dest_path));
    report.notes.push("文件转换成功".to_string());
    Ok(report)
}

/// 支持信息（Supporting Information）汇总
pub fn si(inputs: &[PathBuf], output: Option<&Path>) -> ServiceResult<Report> {
    if inputs.is_empty() {
        return Err("至少需要指定一个输入文件".to_string());
    }
    let paths: Vec<String> = inputs.iter().map(|path| absolute_str(path, "输入文件")).collect::<ServiceResult<_>>()?;
    let text = Tool::new().get_SI(&paths);

    let mut report = Report::new("支持信息（Supporting Information）");
    report.fields.push(Field::new("输入文件数", inputs.len().to_string()));
    match output {
        Some(output) => {
            let dest_path = absolute_str(output, "输出文件")?;
            std::fs::write(&dest_path, &text).map_err(|err| format!("写入支持信息文件失败: {err}"))?;
            report.fields.push(Field::new("输出文件", dest_path));
        }
        None => report.notes.push(text),
    }
    Ok(report)
}
