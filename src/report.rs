//! 结构化分析结果。
//!
//! 各分析模块只负责构造 [`Report`]，不直接打印。CLI 用 [`Display`] 输出
//! 可读文本，用 [`Report::to_json`] 输出 JSON；TUI 复用同一结构渲染界面。
//!
//! 表格单元格统一用已格式化的字符串，这样文本对齐、JSON 序列化和 TUI 渲染
//! 共享同一份精度设置。

use std::fmt;

/// 表格列的对齐方式。
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Align {
    /// 左对齐（文字）
    Left,
    /// 右对齐（数值）
    Right,
    /// 居中
    Center,
}

/// 表格列定义。
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Column {
    /// 列标题
    pub name: String,
    /// 对齐方式
    pub align: Align,
}

impl Column {
    /// 创建左对齐列
    pub fn left(name: impl Into<String>) -> Self {
        Self { name: name.into(), align: Align::Left }
    }

    /// 创建右对齐列
    pub fn right(name: impl Into<String>) -> Self {
        Self { name: name.into(), align: Align::Right }
    }

    /// 创建居中列
    pub fn center(name: impl Into<String>) -> Self {
        Self { name: name.into(), align: Align::Center }
    }
}

/// 结果表格：列定义 + 已格式化的行。
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Table {
    /// 列定义
    pub columns: Vec<Column>,
    /// 数据行，每行的单元格数量应与列数一致
    pub rows: Vec<Vec<String>>,
}

impl Table {
    /// 按列定义创建空表
    pub fn new(columns: Vec<Column>) -> Self {
        Self { columns, rows: Vec::new() }
    }

    /// 追加一行
    pub fn push_row<R, S>(&mut self, cells: R)
    where
        R: IntoIterator<Item = S>,
        S: Into<String>,
    {
        self.rows.push(cells.into_iter().map(Into::into).collect());
    }

    /// 是否没有可显示的内容
    pub fn is_empty(&self) -> bool {
        self.columns.is_empty() || self.rows.is_empty()
    }
}

/// 标量字段（文件、电荷、能量等键值对）。
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Field {
    /// 字段名
    pub label: String,
    /// 字段值（已格式化）
    pub value: String,
}

impl Field {
    /// 创建字段
    pub fn new(label: impl Into<String>, value: impl Into<String>) -> Self {
        Self { label: label.into(), value: value.into() }
    }
}

/// 一次分析的完整结果。
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Report {
    /// 标题
    pub title: String,
    /// 标量字段
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub fields: Vec<Field>,
    /// 结果表格（通常一个，可能有多个）
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub tables: Vec<Table>,
    /// 备注（警告、导出路径、附加说明等）
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub notes: Vec<String>,
}

impl Report {
    /// 创建只带标题的空结果
    pub fn new(title: impl Into<String>) -> Self {
        Self { title: title.into(), ..Self::default() }
    }

    /// 追加一个标量字段
    pub fn field(mut self, label: impl Into<String>, value: impl Into<String>) -> Self {
        self.fields.push(Field::new(label, value));
        self
    }

    /// 追加一张表格
    pub fn table(mut self, table: Table) -> Self {
        self.tables.push(table);
        self
    }

    /// 追加一条备注
    pub fn note(mut self, note: impl Into<String>) -> Self {
        self.notes.push(note.into());
        self
    }

    /// 序列化为缩进后的 JSON
    pub fn to_json(&self) -> String {
        serde_json::to_string_pretty(self).expect("结果序列化为 JSON 失败")
    }

    /// 序列化为单行 JSON
    pub fn to_json_compact(&self) -> String {
        serde_json::to_string(self).expect("结果序列化为 JSON 失败")
    }
}

/// 单个字符在终端中占用的列数（CJK 汉字、假名、全角符号等占 2 列）
fn char_width(c: char) -> usize {
    let wide = matches!(
        c as u32,
        0x1100..=0x115F
            | 0x2E80..=0x303E
            | 0x3041..=0x33FF
            | 0x3400..=0x4DBF
            | 0x4E00..=0x9FFF
            | 0xA000..=0xA4CF
            | 0xAC00..=0xD7A3
            | 0xF900..=0xFAFF
            | 0xFE30..=0xFE6F
            | 0xFF00..=0xFF60
            | 0xFFE0..=0xFFE6
            | 0x1F300..=0x1F64F
            | 0x1F900..=0x1F9FF
            | 0x20000..=0x3FFFD
    );
    if wide { 2 } else { 1 }
}

/// 字符串在终端中的显示宽度
pub fn disp_width(s: &str) -> usize {
    s.chars().map(char_width).sum()
}

/// 按显示宽度补齐，用于表格单元格
fn pad_cell(s: &str, width: usize, align: Align) -> String {
    let w = disp_width(s);
    if w >= width {
        return s.to_string();
    }
    let gap = width - w;
    match align {
        Align::Left => format!("{s}{}", " ".repeat(gap)),
        Align::Right => format!("{}{s}", " ".repeat(gap)),
        Align::Center => {
            let left = gap / 2;
            format!("{}{s}{}", " ".repeat(left), " ".repeat(gap - left))
        }
    }
}

/// 每列的实际宽度
fn column_widths(table: &Table) -> Vec<usize> {
    let mut widths: Vec<usize> = table.columns.iter().map(|c| disp_width(&c.name)).collect();
    for row in &table.rows {
        for (i, width) in widths.iter_mut().enumerate() {
            if let Some(cell) = row.get(i) {
                *width = (*width).max(disp_width(cell));
            }
        }
    }
    widths
}

/// 表格总宽度（列宽之和 + 列间距）
fn table_width(widths: &[usize]) -> usize {
    if widths.is_empty() { 0 } else { widths.iter().sum::<usize>() + 2 * (widths.len() - 1) }
}

fn write_table(f: &mut fmt::Formatter<'_>, table: &Table, widths: &[usize]) -> fmt::Result {
    let headers: Vec<String> = table.columns.iter().zip(widths).map(|(col, &w)| pad_cell(&col.name, w, col.align)).collect();
    writeln!(f, "{}", headers.join("  "))?;

    let total = table_width(widths);
    writeln!(f, "{}", "-".repeat(total))?;

    for row in &table.rows {
        let cells: Vec<String> = table
            .columns
            .iter()
            .zip(widths)
            .enumerate()
            .map(|(i, (col, &w))| pad_cell(row.get(i).map(String::as_str).unwrap_or(""), w, col.align))
            .collect();
        writeln!(f, "{}", cells.join("  "))?;
    }
    Ok(())
}

impl fmt::Display for Report {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let label_width = self.fields.iter().map(|field| disp_width(&field.label)).max().unwrap_or(0);
        let layouts: Vec<Vec<usize>> = self.tables.iter().map(column_widths).collect();

        // 分隔线宽度取标题、字段和表格行宽的最大值
        let mut width = disp_width(&self.title).max(40);
        for field in &self.fields {
            width = width.max(label_width + 2 + disp_width(&field.value));
        }
        for layout in &layouts {
            width = width.max(table_width(layout));
        }

        writeln!(f)?;
        writeln!(f, "{}", self.title)?;
        writeln!(f, "{}", "=".repeat(width))?;

        for field in &self.fields {
            let gap = " ".repeat(label_width - disp_width(&field.label));
            writeln!(f, "{}{}: {}", field.label, gap, field.value)?;
        }

        for (i, (table, layout)) in self.tables.iter().zip(&layouts).enumerate() {
            if i > 0 || !self.fields.is_empty() {
                writeln!(f)?;
            }
            write_table(f, table, layout)?;
        }

        for note in &self.notes {
            writeln!(f, "{}", note)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> Report {
        let mut table = Table::new(vec![Column::right("原子"), Column::left("符号")]);
        table.push_row(vec!["1", "C"]);
        table.push_row(vec!["2", "H"]);
        Report::new("Mulliken 电子布居").field("电荷", "0").table(table).note("备注")
    }

    #[test]
    fn json_roundtrip() {
        let report = sample();
        let back: Report = serde_json::from_str(&report.to_json()).expect("反序列化失败");
        assert_eq!(report, back);
    }

    #[test]
    fn text_contains_content() {
        let text = sample().to_string();
        assert!(text.contains("Mulliken 电子布居"));
        assert!(text.contains("原子"));
        assert!(text.contains("备注"));
    }
}
