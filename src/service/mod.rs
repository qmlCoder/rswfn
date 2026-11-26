//! 分析任务层。
//!
//! CLI 与 TUI 都只负责收集参数，实际计算在这里完成并返回
//! [`Report`](crate::report::Report)。函数不再直接 `println!`，也不会在
//! 可预期的错误上 `panic!`，而是返回 `Err(String)` 交给前端展示。

use std::path::{Path, absolute};

use crate::base::Mole;
use crate::reader::{FchReader, GjfReader, LogReader, MfnReader};
use crate::report::{Field, Report};

pub mod atomprop;
pub mod bondprop;
pub mod file;
pub mod gridprop;
pub mod moleprop;
pub mod orbtprop;
pub mod qmcal;

/// 服务层统一错误类型：给用户看的中文错误信息。
pub type ServiceResult<T> = Result<T, String>;

/// 读取带波函数信息的文件（log/out/gjf/fch/fchk/mwfn）。
///
/// 与旧的 CLI 行为一致，会记录绝对路径，供能量分解等需要回读原文件的分析使用。
pub fn load_mole(filepath: &Path) -> ServiceResult<Mole> {
    let path = absolute(filepath).map_err(|err| format!("无法解析输入路径 {}: {err}", filepath.display()))?;
    if !path.is_file() {
        return Err(format!("文件不存在或不是普通文件: {}", path.display()));
    }
    let path_str = path.to_str().ok_or_else(|| format!("路径不是有效的 UTF-8: {}", path.display()))?.to_string();
    let ftype = ext_of(&path).unwrap_or_default();
    let mut mole = match ftype.as_str() {
        "log" | "out" => Mole::from_read(&LogReader::new(&path_str)),
        "gjf" => Mole::from_read(&GjfReader::new(&path_str)),
        "fch" | "fchk" => Mole::from_read(&FchReader::new(&path_str)),
        "mwfn" => Mole::from_read(&MfnReader::new(&path_str)),
        other => return Err(format!("不支持的波函数文件类型: .{other}，支持: log/out/gjf/fch/fchk/mwfn")),
    };
    mole.path = Some(path_str);
    Ok(mole)
}

/// 读取任意 `Mole::from_file` 支持的结构/波函数文件（用于坐标提取与格式转换）。
pub fn load_any_mole(filepath: &Path) -> ServiceResult<Mole> {
    let path = absolute(filepath).map_err(|err| format!("无法解析输入路径 {}: {err}", filepath.display()))?;
    if !path.is_file() {
        return Err(format!("文件不存在或不是普通文件: {}", path.display()));
    }
    let path_str = path.to_str().ok_or_else(|| format!("路径不是有效的 UTF-8: {}", path.display()))?;
    Ok(Mole::from_file(path_str))
}

/// 小写扩展名
pub(crate) fn ext_of(path: &Path) -> Option<String> {
    path.extension().map(|ext| ext.to_string_lossy().to_ascii_lowercase())
}

/// 检查分子是否含基组与轨道系数，避免计算器内部 panic。
pub(crate) fn require_coefs(mole: &Mole) -> ServiceResult<()> {
    if mole.basis.is_none() {
        return Err("该分子不含基组信息，请使用 fch/fchk/log/out/mwfn 文件".to_string());
    }
    if mole.coefs.is_none() {
        return Err("该分子不含轨道系数，请使用 fch/fchk/log/out/mwfn 文件".to_string());
    }
    Ok(())
}

/// 检查 1 基原子序号是否在范围内。
pub(crate) fn check_atom(mole: &Mole, idx: usize) -> ServiceResult<()> {
    let natm = mole.atoms().len();
    if idx == 0 || idx > natm {
        return Err(format!("原子序号越界: {idx}（分子共有 {natm} 个原子，序号从 1 开始）"));
    }
    Ok(())
}

/// 分子基本信息字段，所有分析结果共用。
pub(crate) fn molecule_fields(mole: &Mole) -> Vec<Field> {
    vec![
        Field::new("文件", mole.path.clone().unwrap_or_else(|| "<内存分子>".to_string())),
        Field::new("分子式", mole.formula()),
        Field::new("电荷", mole.charge().to_string()),
        Field::new("多重度", mole.multip().to_string()),
        Field::new("原子数", mole.atoms().len().to_string()),
    ]
}

/// 创建带分子基本信息的结果对象。
pub(crate) fn report_for(title: impl Into<String>, mole: &Mole) -> Report {
    let mut report = Report::new(title);
    report.fields = molecule_fields(mole);
    report
}
