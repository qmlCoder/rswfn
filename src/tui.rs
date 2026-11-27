//! 全屏交互界面（基于 ratatui）。
//!
//! 设计要点：
//! - 菜单/表单由 [`TASKS`] 元数据驱动，CLI 的每个分析任务都在这里有对应入口；
//! - 计算在后台线程执行，结果通过通道回到界面，避免长任务卡死 UI；
//! - 进度信息走 `log`，在 TUI 中收集到日志面板，不会污染全屏画面；
//! - 所有可预期错误都变成弹窗，不会退出程序。

use std::collections::{HashSet, VecDeque};
use std::io::{self, Write};
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::path::{Path, PathBuf};
use std::sync::mpsc::{self, Receiver, TryRecvError};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::{Duration, Instant};

use ratatui::backend::CrosstermBackend;
use ratatui::crossterm::event::{self, DisableBracketedPaste, EnableBracketedPaste, Event, KeyCode, KeyEvent, KeyEventKind, KeyModifiers};
use ratatui::crossterm::execute;
use ratatui::crossterm::terminal::{EnterAlternateScreen, LeaveAlternateScreen, disable_raw_mode, enable_raw_mode};
use ratatui::layout::{Alignment, Constraint, Layout, Rect};
use ratatui::style::{Color, Modifier, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, Cell, Clear, List, ListItem, ListState, Paragraph, Row, Table, TableState, Wrap};
use ratatui::{Frame, Terminal};

use crate::report::Report;
use crate::service::{self, atomprop, bondprop, file, gridprop, moleprop, orbtprop, qmcal};
use crate::viewer::{self, View};

// ===== 选项常量 =====

const CTYPES: &[&str] = &["mulliken", "lowdin", "hirshfeld"];
const METHODS: &[&str] = &["pocv", "mocv"];
const OTYPES: &[&str] = &["mayer", "wiberg"];
const UNITS: &[&str] = &["ang", "bohr"];
const TOPO_FUNCS: &[&str] = &["mol_rho_dm", "pi_rho_pocv", "pi_rho_mocv"];
const DFT_FUNCS: &[&str] = &["b3lyp", "pbe0", "pbe", "blyp", "bp86", "svwn5"];
const DENSITY_FUNCS: &[&str] = &[
    "density",
    "pro_rho",
    "spin",
    "rho_a",
    "rho_b",
    "del_rho",
    "pi_pocv",
    "pi_mocv",
    "rdg",
    "iri",
    "laplace",
    "signl2rho",
    "ked",
    "lagkin",
    "hamkin",
    "elf",
    "lol",
    "infoentro",
    "activate",
    "curvature",
    "fukui_plus",
    "fukui_minus",
    "fukui_zero",
    "fukui_df",
];
const EXTENSIONS: &[&str] = &["log", "out", "fch", "fchk", "gjf", "mwfn", "xyz", "mol", "sdf", "mde", "dmi", "dmo", "inp", "cub"];

// ===== 任务元数据 =====

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TaskId {
    Mulliken,
    Lowdin,
    Hirshfeld,
    Spin,
    Resp,
    PiPop,
    Fukui,
    FukuiPi,
    Freev,
    FreevPi,
    Mayer,
    Wiberg,
    OrderObt,
    PiOrder,
    Decom,
    BondEnergy,
    Aromacity,
    PiEnergy,
    PiElesPocv,
    PiElesMocv,
    Composition,
    GridDensity,
    GridPotential,
    GridTopology,
    QmHf,
    QmDft,
    Coord,
    Trans,
    Dhmo,
    Si,
}

struct TaskInfo {
    id: TaskId,
    name: &'static str,
    desc: &'static str,
    /// 对应的 CLI 子命令路径
    cli: &'static str,
}

static TASKS: &[TaskInfo] = &[
    TaskInfo {
        id: TaskId::Mulliken,
        cli: "atomprop charge mulliken",
        name: "Mulliken 布居",
        desc: "Mulliken 电子布居分析",
    },
    TaskInfo {
        id: TaskId::Lowdin,
        cli: "atomprop charge lowdin",
        name: "Löwdin 布居",
        desc: "Löwdin 电子布居分析",
    },
    TaskInfo {
        id: TaskId::Hirshfeld,
        cli: "atomprop charge hirshfeld",
        name: "Hirshfeld 布居",
        desc: "Hirshfeld 电子布居分析",
    },
    TaskInfo {
        id: TaskId::Spin,
        cli: "atomprop charge spin",
        name: "自旋布居",
        desc: "自旋密度布居（仅开壳层）",
    },
    TaskInfo {
        id: TaskId::Resp,
        cli: "atomprop charge resp",
        name: "RESP 电荷",
        desc: "限制性静电势电荷拟合",
    },
    TaskInfo {
        id: TaskId::PiPop,
        cli: "atomprop charge pi-pop",
        name: "π 电子布居",
        desc: "POCV / MOCV 原子 π 电子布居",
    },
    TaskInfo {
        id: TaskId::Fukui,
        cli: "atomprop activity fukui",
        name: "福井函数",
        desc: "福井函数与双描述符（需 N/N+1/N-1）",
    },
    TaskInfo {
        id: TaskId::FukuiPi,
        cli: "atomprop activity fukui-pi",
        name: "π 福井函数",
        desc: "π 轨道福井函数与双描述符",
    },
    TaskInfo {
        id: TaskId::Freev,
        cli: "atomprop activity freev",
        name: "自由价",
        desc: "自由价（free valence）",
    },
    TaskInfo {
        id: TaskId::FreevPi,
        cli: "atomprop activity freev-pi",
        name: "π 自由价",
        desc: "π 自由价（POCV）",
    },
    TaskInfo {
        id: TaskId::Mayer,
        cli: "bondprop order mayer",
        name: "Mayer 键级",
        desc: "Mayer 键级分析",
    },
    TaskInfo {
        id: TaskId::Wiberg,
        cli: "bondprop order wiberg",
        name: "Wiberg 键级",
        desc: "Wiberg 键级分析",
    },
    TaskInfo {
        id: TaskId::OrderObt,
        cli: "bondprop order obt",
        name: "轨道键级",
        desc: "指定分子轨道对键级的贡献",
    },
    TaskInfo {
        id: TaskId::PiOrder,
        cli: "bondprop order pi-order",
        name: "π 键级",
        desc: "POCV / MOCV π 键级",
    },
    TaskInfo {
        id: TaskId::Decom,
        cli: "bondprop order decom",
        name: "键级分解",
        desc: "σ/π/δ/φ 键级分量",
    },
    TaskInfo {
        id: TaskId::BondEnergy,
        cli: "bondprop energy",
        name: "键能量分解",
        desc: "bond / rest / cross 能量分解",
    },
    TaskInfo {
        id: TaskId::Aromacity,
        cli: "moleprop aromacity",
        name: "芳香性",
        desc: "基于 π 键级的 PISD 芳香性",
    },
    TaskInfo {
        id: TaskId::PiEnergy,
        cli: "moleprop pi-energy",
        name: "π 电子能量",
        desc: "π 电子能量分解",
    },
    TaskInfo {
        id: TaskId::PiElesPocv,
        cli: "orbtprop pi-eles-pocv",
        name: "轨道 π 电子数 (POCV)",
        desc: "每个分子轨道的 π 电子数",
    },
    TaskInfo {
        id: TaskId::PiElesMocv,
        cli: "orbtprop pi-eles-mocv",
        name: "轨道 π 电子数 (MOCV)",
        desc: "每个分子轨道的 π 电子数",
    },
    TaskInfo {
        id: TaskId::Composition,
        cli: "orbtprop composition",
        name: "轨道成分",
        desc: "指定轨道的原子电子布居",
    },
    TaskInfo {
        id: TaskId::GridDensity,
        cli: "gridprop density",
        name: "实空间函数 cub",
        desc: "密度、RDG、ELF、福井函数等导出 cub",
    },
    TaskInfo {
        id: TaskId::GridPotential,
        cli: "gridprop potential",
        name: "静电势 cub",
        desc: "静电势导出 cub",
    },
    TaskInfo {
        id: TaskId::GridTopology,
        cli: "gridprop topology",
        name: "拓扑临界点",
        desc: "寻找核/键/环/笼临界点",
    },
    TaskInfo {
        id: TaskId::QmHf,
        cli: "qmcal hf",
        name: "Hartree-Fock",
        desc: "闭壳层 RHF 计算",
    },
    TaskInfo {
        id: TaskId::QmDft,
        cli: "qmcal dft",
        name: "DFT",
        desc: "RKS / UKS 密度泛函计算",
    },
    TaskInfo {
        id: TaskId::Coord,
        cli: "coord",
        name: "提取坐标",
        desc: "打印/导出分子坐标",
    },
    TaskInfo {
        id: TaskId::Trans,
        cli: "trans",
        name: "格式转换",
        desc: "log/out→gjf、gjf→dmi",
    },
    TaskInfo {
        id: TaskId::Dhmo,
        cli: "qmcal dhmo",
        name: "DHMO 计算",
        desc: "简单 Hückel 计算",
    },
    TaskInfo {
        id: TaskId::Si,
        cli: "si",
        name: "支持信息汇总",
        desc: "汇总能量、频率与结构",
    },
];

fn task_info(id: TaskId) -> &'static TaskInfo {
    TASKS.iter().find(|task| task.id == id).expect("任务表缺少条目")
}

// ===== 分析菜单树 =====

/// 菜单元数据；层级严格对应 `src/cli.rs` 中的子命令结构。
#[derive(Debug, Clone, Copy)]
enum NodeKind {
    Group(&'static str),
    Task(TaskId),
    /// CLI 的 `view` 子命令：打开结构预览窗口
    Viewer,
}

#[derive(Debug)]
struct Node {
    kind: NodeKind,
    children: &'static [Node],
}

const fn leaf(task: TaskId) -> Node {
    Node {
        kind: NodeKind::Task(task),
        children: &[],
    }
}

/// CLI 的 `view` 子命令：在 TUI 中打开结构预览窗口
const fn viewer() -> Node {
    Node { kind: NodeKind::Viewer, children: &[] }
}

const fn branch(label: &'static str, children: &'static [Node]) -> Node {
    Node { kind: NodeKind::Group(label), children }
}

static TREE: &[Node] = &[
    branch(
        "atomprop · 原子属性",
        &[
            branch(
                "charge · 电荷",
                &[leaf(TaskId::Mulliken), leaf(TaskId::Lowdin), leaf(TaskId::Hirshfeld), leaf(TaskId::Spin), leaf(TaskId::Resp), leaf(TaskId::PiPop)],
            ),
            branch("activity · 活性", &[leaf(TaskId::Fukui), leaf(TaskId::FukuiPi), leaf(TaskId::Freev), leaf(TaskId::FreevPi)]),
        ],
    ),
    branch(
        "bondprop · 键属性",
        &[
            branch("order · 键级", &[leaf(TaskId::Mayer), leaf(TaskId::Wiberg), leaf(TaskId::OrderObt), leaf(TaskId::PiOrder), leaf(TaskId::Decom)]),
            leaf(TaskId::BondEnergy),
        ],
    ),
    branch("moleprop · 分子属性", &[leaf(TaskId::Aromacity), leaf(TaskId::PiEnergy)]),
    branch("orbtprop · 轨道属性", &[leaf(TaskId::PiElesPocv), leaf(TaskId::PiElesMocv), leaf(TaskId::Composition)]),
    branch("gridprop · 格点属性", &[leaf(TaskId::GridDensity), leaf(TaskId::GridPotential), leaf(TaskId::GridTopology)]),
    leaf(TaskId::Coord),
    leaf(TaskId::Trans),
    branch("qmcal · 量子化学", &[leaf(TaskId::QmHf), leaf(TaskId::QmDft), leaf(TaskId::Dhmo)]),
    leaf(TaskId::Si),
    viewer(),
];

/// 可见树的一行
#[derive(Debug, Clone, Copy)]
struct TreeRow {
    node: &'static Node,
    depth: usize,
    parent: Option<usize>,
    is_group: bool,
    expanded: bool,
}

fn node_key(node: &'static Node) -> usize {
    node as *const Node as usize
}

fn node_label(node: &'static Node) -> &'static str {
    match node.kind {
        NodeKind::Group(label) => label,
        NodeKind::Task(task) => task_info(task).name,
        NodeKind::Viewer => "结构预览（view）",
    }
}

fn node_task(node: &'static Node) -> Option<TaskId> {
    match node.kind {
        NodeKind::Task(task) => Some(task),
        NodeKind::Group(_) | NodeKind::Viewer => None,
    }
}

fn flatten(nodes: &'static [Node], expanded: &std::collections::HashSet<usize>, depth: usize, parent: Option<usize>, out: &mut Vec<TreeRow>) {
    for node in nodes {
        let index = out.len();
        let is_group = !node.children.is_empty();
        let is_expanded = expanded.contains(&node_key(node));
        out.push(TreeRow {
            node,
            depth,
            parent,
            is_group,
            expanded: is_expanded,
        });
        if is_group && is_expanded {
            flatten(node.children, expanded, depth + 1, Some(index), out);
        }
    }
}

// ===== 表单字段 =====

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FileNeed {
    Required,
    Optional,
}

#[derive(Debug, Clone, Copy)]
enum FieldKind {
    Choice(&'static [&'static str], usize),
    Text(&'static str),
    Int(&'static str),
    Float(&'static str),
    Bool(bool),
    File { need: FileNeed, prefill: bool },
}

#[derive(Debug, Clone, Copy)]
struct FieldSpec {
    label: &'static str,
    kind: FieldKind,
}

fn spec(label: &'static str, kind: FieldKind) -> FieldSpec {
    FieldSpec { label, kind }
}

fn choice(options: &'static [&'static str], default: usize) -> FieldKind {
    FieldKind::Choice(options, default)
}

fn text(default: &'static str) -> FieldKind {
    FieldKind::Text(default)
}

fn int(default: &'static str) -> FieldKind {
    FieldKind::Int(default)
}

fn float(default: &'static str) -> FieldKind {
    FieldKind::Float(default)
}

fn boolean(default: bool) -> FieldKind {
    FieldKind::Bool(default)
}

/// 必填文件，默认填入当前打开的文件
fn file_req() -> FieldKind {
    FieldKind::File { need: FileNeed::Required, prefill: true }
}

/// 必填文件，不预填（例如福井函数需要的 N+1 / N-1 文件）
fn mole_req() -> FieldKind {
    FieldKind::File { need: FileNeed::Required, prefill: false }
}

/// 可选文件
fn file_opt() -> FieldKind {
    FieldKind::File { need: FileNeed::Optional, prefill: false }
}

fn field_specs(id: TaskId) -> Vec<FieldSpec> {
    match id {
        TaskId::Mulliken | TaskId::Lowdin | TaskId::Hirshfeld | TaskId::Mayer | TaskId::Wiberg | TaskId::PiElesPocv | TaskId::PiElesMocv => vec![],
        TaskId::Spin => vec![spec("布居方法", choice(CTYPES, 0))],
        TaskId::Resp => vec![spec("MK 点密度", float("6.0")), spec("仅第一阶段", boolean(false))],
        TaskId::PiPop => vec![spec("布居方法", choice(CTYPES, 0)), spec("轨道改造", choice(METHODS, 0))],
        TaskId::Fukui => vec![spec("阴离子 N+1", mole_req()), spec("阳离子 N-1", mole_req()), spec("布居方法", choice(CTYPES, 0))],
        TaskId::FukuiPi => vec![
            spec("阴离子 N+1", mole_req()),
            spec("阳离子 N-1", mole_req()),
            spec("布居方法", choice(CTYPES, 0)),
            spec("轨道改造", choice(METHODS, 0)),
        ],
        TaskId::Freev | TaskId::FreevPi => vec![spec("键级类型", choice(OTYPES, 0))],
        TaskId::OrderObt => vec![spec("轨道序号", int("1")), spec("键级类型", choice(OTYPES, 0))],
        TaskId::PiOrder => vec![spec("轨道改造", choice(METHODS, 0)), spec("键级类型", choice(OTYPES, 0))],
        TaskId::Decom | TaskId::BondEnergy => vec![spec("原子1", int("1")), spec("原子2", int("2"))],
        TaskId::Aromacity => vec![spec("环原子（留空自动检测）", text(""))],
        TaskId::PiEnergy => vec![spec("轨道改造", choice(METHODS, 0))],
        TaskId::Composition => vec![spec("轨道序号", int("1")), spec("布居方法", choice(CTYPES, 0))],
        TaskId::GridDensity => vec![
            spec("函数", choice(DENSITY_FUNCS, 0)),
            spec("步长 (Bohr)", float("0.2")),
            spec("输出 cub（留空自动命名）", text("")),
            spec("阴离子 N+1（fukui_* 用）", file_opt()),
            spec("阳离子 N-1（fukui_* 用）", file_opt()),
        ],
        TaskId::GridPotential => vec![spec("步长 (Bohr)", float("0.2")), spec("输出 cub（留空自动命名）", text(""))],
        TaskId::GridTopology => vec![spec("函数", choice(TOPO_FUNCS, 0))],
        TaskId::QmHf => vec![spec("输入文件", file_req()), spec("基组（留空用文件内基组）", text("")), spec("输出 fch", text(""))],
        TaskId::QmDft => vec![
            spec("输入文件", file_req()),
            spec("基组（留空用文件内基组）", text("")),
            spec("泛函", choice(DFT_FUNCS, 0)),
            spec("输出 fch", text("")),
        ],
        TaskId::Coord => vec![spec("坐标单位", choice(UNITS, 0)), spec("输出 xyz", text(""))],
        TaskId::Trans => vec![spec("输入文件", file_req()), spec("输出文件", text(""))],
        TaskId::Dhmo => vec![spec("输入 .dmi", file_req()), spec("输出 .dmo", text(""))],
        TaskId::Si => vec![spec("输入文件（用 ; 分隔）", text("")), spec("输出文本", text(""))],
    }
}

// ===== 文本输入框 =====

#[derive(Debug, Clone, Default)]
struct TextInput {
    value: String,
    cursor: usize,
}

impl TextInput {
    fn new(value: &str) -> Self {
        Self {
            value: value.to_string(),
            cursor: value.len(),
        }
    }

    fn insert(&mut self, ch: char) {
        self.value.insert(self.cursor, ch);
        self.cursor += ch.len_utf8();
    }

    fn backspace(&mut self) {
        if self.cursor > 0 {
            let prev = self.prev_boundary();
            self.value.replace_range(prev..self.cursor, "");
            self.cursor = prev;
        }
    }

    fn delete(&mut self) {
        if self.cursor < self.value.len() {
            let next = self.next_boundary();
            self.value.replace_range(self.cursor..next, "");
        }
    }

    fn left(&mut self) {
        if self.cursor > 0 {
            self.cursor = self.prev_boundary();
        }
    }

    fn right(&mut self) {
        if self.cursor < self.value.len() {
            self.cursor = self.next_boundary();
        }
    }

    fn home(&mut self) {
        self.cursor = 0;
    }

    fn end(&mut self) {
        self.cursor = self.value.len();
    }

    fn prev_boundary(&self) -> usize {
        self.value[..self.cursor].char_indices().last().map(|(i, _)| i).unwrap_or(0)
    }

    fn next_boundary(&self) -> usize {
        self.value[self.cursor..].chars().next().map(|c| self.cursor + c.len_utf8()).unwrap_or(self.value.len())
    }
}

// ===== 表单 =====

#[derive(Debug, Clone)]
enum Field {
    Choice { label: &'static str, options: &'static [&'static str], index: usize },
    Text { label: &'static str, input: TextInput },
    Bool { label: &'static str, value: bool },
    File { label: &'static str, input: TextInput, need: FileNeed },
}

impl Field {
    fn from_spec(spec: &FieldSpec, current: Option<&Path>) -> Self {
        match spec.kind {
            FieldKind::Choice(options, index) => Field::Choice { label: spec.label, options, index },
            FieldKind::Bool(value) => Field::Bool { label: spec.label, value },
            FieldKind::Text(default) | FieldKind::Int(default) | FieldKind::Float(default) => Field::Text {
                label: spec.label,
                input: TextInput::new(default),
            },
            FieldKind::File { need, prefill } => {
                let initial = if prefill { current.map(|p| p.to_string_lossy().to_string()).unwrap_or_default() } else { String::new() };
                Field::File {
                    label: spec.label,
                    input: TextInput::new(&initial),
                    need,
                }
            }
        }
    }

    fn label(&self) -> &'static str {
        match self {
            Field::Choice { label, .. } | Field::Text { label, .. } | Field::Bool { label, .. } | Field::File { label, .. } => label,
        }
    }

    fn text(&self) -> &str {
        match self {
            Field::Choice { options, index, .. } => options.get(*index).copied().unwrap_or(""),
            Field::Text { input, .. } | Field::File { input, .. } => &input.value,
            Field::Bool { value, .. } => {
                if *value {
                    "true"
                } else {
                    "false"
                }
            }
        }
    }

    fn value_line(&self, active: bool) -> Line<'static> {
        match self {
            Field::Choice { options, index, .. } => {
                let value = options.get(*index).copied().unwrap_or("");
                Line::from(Span::styled(format!("‹ {value} ›"), Style::default().fg(Color::Yellow)))
            }
            Field::Bool { value, .. } => Line::from(if *value { "[x]" } else { "[ ]" }),
            Field::Text { input, .. } | Field::File { input, .. } => {
                if !active {
                    return Line::from(input.value.clone());
                }
                let (before, after) = input.value.split_at(input.cursor.min(input.value.len()));
                Line::from(vec![Span::raw(before.to_string()), Span::styled("│", Style::default().fg(Color::Yellow)), Span::raw(after.to_string())])
            }
        }
    }

    fn editable(&self) -> bool {
        matches!(self, Field::Text { .. } | Field::File { .. })
    }
}

#[derive(Debug, Clone)]
struct Form {
    task: TaskId,
    title: &'static str,
    desc: &'static str,
    fields: Vec<Field>,
    selected: usize,
}

impl Form {
    fn new(task: TaskId, current: Option<&Path>) -> Self {
        let info = task_info(task);
        let mut fields: Vec<Field> = field_specs(task).into_iter().map(|spec| Field::from_spec(&spec, current)).collect();
        // “支持信息汇总”的第一个字段是文件列表，默认填入当前文件
        if task == TaskId::Si {
            if let (Some(path), Some(Field::Text { input, .. })) = (current, fields.get_mut(0)) {
                if input.value.is_empty() {
                    *input = TextInput::new(&path.to_string_lossy());
                }
            }
        }
        Self {
            task,
            title: info.name,
            desc: info.desc,
            fields,
            selected: 0,
        }
    }

    fn text(&self, i: usize) -> &str {
        self.fields.get(i).map(Field::text).unwrap_or("")
    }

    fn choice(&self, i: usize) -> &'static str {
        match self.fields.get(i) {
            Some(Field::Choice { options, index, .. }) => options.get(*index).copied().unwrap_or(""),
            _ => "",
        }
    }

    fn flag(&self, i: usize) -> bool {
        matches!(self.fields.get(i), Some(Field::Bool { value: true, .. }))
    }

    fn selected_field(&self) -> Option<&Field> {
        self.fields.get(self.selected)
    }

    fn next(&mut self) {
        if !self.fields.is_empty() {
            self.selected = (self.selected + 1) % self.fields.len();
        }
    }

    fn prev(&mut self) {
        if !self.fields.is_empty() {
            self.selected = (self.selected + self.fields.len() - 1) % self.fields.len();
        }
    }

    fn edit(&mut self, key: KeyEvent) {
        let Some(field) = self.fields.get_mut(self.selected) else { return };
        match field {
            Field::Choice { options, index, .. } => match key.code {
                KeyCode::Left => *index = (*index + options.len() - 1) % options.len(),
                KeyCode::Right | KeyCode::Char(' ') => *index = (*index + 1) % options.len(),
                _ => {}
            },
            Field::Bool { value, .. } => {
                if matches!(key.code, KeyCode::Char(' ')) {
                    *value = !*value;
                }
            }
            Field::Text { input, .. } | Field::File { input, .. } => match key.code {
                KeyCode::Char(ch) => input.insert(ch),
                KeyCode::Backspace => input.backspace(),
                KeyCode::Delete => input.delete(),
                KeyCode::Left => input.left(),
                KeyCode::Right => input.right(),
                KeyCode::Home => input.home(),
                KeyCode::End => input.end(),
                _ => {}
            },
        }
    }

    fn hint(&self) -> &'static str {
        match self.selected_field() {
            Some(Field::Choice { .. }) => "←/→ 切换选项 · Enter 运行 · ↑/↓ 或 Tab 切换字段 · Esc 返回",
            Some(Field::Bool { .. }) => "空格 切换 · Enter 运行 · ↑/↓ 或 Tab 切换字段 · Esc 返回",
            _ => "输入文字 · ←/→ 移动光标 · Enter 运行 · Ctrl+O 选择文件 · Esc 返回",
        }
    }
}

// ===== 结果与任务 =====

#[derive(Debug)]
enum JobOutcome {
    Ok(Report),
    Err(String),
}

struct Job {
    label: String,
    started: Instant,
    rx: Receiver<JobOutcome>,
    cancelled: bool,
}

#[derive(Debug, Clone)]
struct Header {
    formula: String,
    charge: String,
    multip: String,
    natoms: String,
    /// 供 F4 结构预览使用的分子快照
    mol: viewer::Molecule,
}

impl Header {
    /// 读取文件并提取界面头部信息（顺带准备结构预览用的分子快照）
    fn from_path(path: &Path) -> Option<Header> {
        let mole = service::load_mole(path).or_else(|_| service::load_any_mole(path)).ok()?;
        Some(Header {
            formula: mole.formula(),
            charge: mole.charge().to_string(),
            multip: mole.multip().to_string(),
            natoms: mole.atoms().len().to_string(),
            mol: viewer::Molecule::from_mole(&mole),
        })
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Mode {
    Welcome,
    Form,
    Running,
    Result,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Focus {
    Files,
    Tree,
}

// ===== 文件浏览器 =====

#[derive(Debug, Clone)]
struct Entry {
    name: String,
    path: PathBuf,
    is_dir: bool,
}

#[derive(Debug, Clone, Copy)]
enum BrowserTarget {
    CurrentFile,
    Field(usize),
    None,
}

#[derive(Debug)]
enum BrowserEvent {
    Stay,
    Close,
    Choose(PathBuf),
}

/// 文件列表（左侧面板与弹出浏览器共用）
struct FilePanel {
    dir: PathBuf,
    entries: Vec<Entry>,
    selected: usize,
    state: ListState,
    error: Option<String>,
}

impl FilePanel {
    fn new(dir: PathBuf) -> Self {
        let mut panel = Self {
            dir,
            entries: Vec::new(),
            selected: 0,
            state: ListState::default(),
            error: None,
        };
        if !panel.dir.is_dir() {
            panel.dir = std::env::current_dir().unwrap_or_else(|_| PathBuf::from("."));
        }
        panel.refresh();
        panel
    }

    fn refresh(&mut self) {
        match read_entries(&self.dir) {
            Ok(entries) => {
                self.entries = entries;
                self.error = None;
            }
            Err(err) => {
                self.entries.clear();
                self.error = Some(format!("无法读取目录: {err}"));
            }
        }
        self.selected = 0;
        self.state.select(Some(0));
    }

    fn enter_dir(&mut self, dir: PathBuf) {
        self.dir = dir;
        self.refresh();
    }

    fn go_parent(&mut self) {
        if let Some(parent) = self.dir.parent() {
            let parent = parent.to_path_buf();
            self.enter_dir(parent);
        }
    }

    fn move_by(&mut self, delta: isize) {
        if self.entries.is_empty() {
            self.selected = 0;
            return;
        }
        let len = self.entries.len() as isize;
        self.selected = ((self.selected as isize + delta).rem_euclid(len)) as usize;
    }

    fn current(&self) -> Option<&Entry> {
        self.entries.get(self.selected)
    }

    fn select_path(&mut self, path: &Path) {
        if let Some(index) = self.entries.iter().position(|entry| entry.path == path) {
            self.selected = index;
        }
    }

    /// 选中目录则进入，选中文件则返回该文件路径
    fn enter_selected(&mut self) -> Option<PathBuf> {
        let entry = self.current().cloned()?;
        if entry.is_dir {
            self.enter_dir(entry.path);
            None
        } else {
            Some(entry.path)
        }
    }

    fn items(&self, current: Option<&Path>) -> Vec<ListItem<'static>> {
        if self.entries.is_empty() {
            let msg = self.error.clone().unwrap_or_else(|| "（空目录）".to_string());
            return vec![ListItem::new(Line::from(Span::styled(msg, Style::default().fg(Color::DarkGray))))];
        }
        self.entries
            .iter()
            .map(|entry| {
                let is_current = current.map(|path| path == entry.path.as_path()).unwrap_or(false);
                let (text, style) = if entry.is_dir {
                    (format!("  {}/", entry.name), Style::default().fg(Color::Cyan))
                } else if is_current {
                    (format!("* {}", entry.name), Style::default().fg(Color::Green).add_modifier(Modifier::BOLD))
                } else {
                    (format!("  {}", entry.name), Style::default())
                };
                ListItem::new(Line::from(Span::styled(text, style)))
            })
            .collect()
    }

    fn dir_title(&self) -> String {
        self.dir.file_name().map(|name| name.to_string_lossy().to_string()).unwrap_or_else(|| self.dir.display().to_string())
    }
}

struct Browser {
    panel: FilePanel,
    target: BrowserTarget,
}

impl Browser {
    fn new(dir: PathBuf, target: BrowserTarget) -> Self {
        Self { panel: FilePanel::new(dir), target }
    }

    fn on_key(&mut self, key: KeyEvent) -> BrowserEvent {
        match key.code {
            KeyCode::Esc | KeyCode::Char('q') => BrowserEvent::Close,
            KeyCode::Up | KeyCode::Char('k') => {
                self.panel.move_by(-1);
                BrowserEvent::Stay
            }
            KeyCode::Down | KeyCode::Char('j') => {
                self.panel.move_by(1);
                BrowserEvent::Stay
            }
            KeyCode::PageUp => {
                self.panel.move_by(-10);
                BrowserEvent::Stay
            }
            KeyCode::PageDown => {
                self.panel.move_by(10);
                BrowserEvent::Stay
            }
            KeyCode::Home => {
                self.panel.selected = 0;
                BrowserEvent::Stay
            }
            KeyCode::End => {
                self.panel.selected = self.panel.entries.len().saturating_sub(1);
                BrowserEvent::Stay
            }
            KeyCode::Backspace | KeyCode::Char('h') => {
                self.panel.go_parent();
                BrowserEvent::Stay
            }
            KeyCode::Enter => match self.panel.enter_selected() {
                Some(path) => BrowserEvent::Choose(path),
                None => BrowserEvent::Stay,
            },
            _ => BrowserEvent::Stay,
        }
    }
}

fn read_entries(dir: &Path) -> io::Result<Vec<Entry>> {
    let mut dirs = Vec::new();
    let mut files = Vec::new();
    for entry in std::fs::read_dir(dir)? {
        let entry = entry?;
        let name = entry.file_name().to_string_lossy().to_string();
        if name.starts_with('.') {
            continue;
        }
        let path = entry.path();
        let is_dir = entry.file_type().map(|kind| kind.is_dir()).unwrap_or(false);
        if is_dir {
            dirs.push(Entry { name, path, is_dir: true });
        } else if is_supported_file(&path) {
            files.push(Entry { name, path, is_dir: false });
        }
    }
    dirs.sort_by_key(|entry| entry.name.to_lowercase());
    files.sort_by_key(|entry| entry.name.to_lowercase());
    dirs.extend(files);
    Ok(dirs)
}

fn is_supported_file(path: &Path) -> bool {
    match path.extension() {
        Some(ext) => {
            let ext = ext.to_string_lossy().to_ascii_lowercase();
            EXTENSIONS.contains(&ext.as_str())
        }
        None => false,
    }
}

// ===== 覆盖层 =====

/// OSC 52 需要标准 base64（不引入依赖，自己编）
fn base64(data: &[u8]) -> String {
    const TABLE: &[u8; 64] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    let mut out = String::with_capacity(data.len().div_ceil(3) * 4);
    for chunk in data.chunks(3) {
        let b = [chunk[0], *chunk.get(1).unwrap_or(&0), *chunk.get(2).unwrap_or(&0)];
        let n = ((b[0] as u32) << 16) | ((b[1] as u32) << 8) | b[2] as u32;
        out.push(TABLE[(n >> 18) as usize & 63] as char);
        out.push(TABLE[(n >> 12) as usize & 63] as char);
        out.push(if chunk.len() > 1 { TABLE[(n >> 6) as usize & 63] as char } else { '=' });
        out.push(if chunk.len() > 2 { TABLE[n as usize & 63] as char } else { '=' });
    }
    out
}

/// 预览配色循环：真彩色 → 256 色 → 无色 ASCII
fn next_depth(depth: viewer::Depth) -> viewer::Depth {
    match depth {
        viewer::Depth::True => viewer::Depth::Ansi256,
        viewer::Depth::Ansi256 => viewer::Depth::Mono,
        viewer::Depth::Mono => viewer::Depth::True,
    }
}

/// 分子结构预览（F4）的界面状态
struct ViewerState {
    title: String,
    mol: viewer::Molecule,
    view: View,
    /// 色彩深度（`c` 键切换，影响预览的着色方式）
    depth: viewer::Depth,
}

/// 四个面板，同一时刻只有一个处于激活状态
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Panel {
    Files,
    Menu,
    Params,
    Result,
}

/// 激活面板的边框颜色：用黄色，避开文件夹/分组用的青色
fn panel_border(active: bool) -> Style {
    if active {
        Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD)
    } else {
        Style::default().fg(Color::DarkGray)
    }
}

/// 一次计算结果的快照：F3 弹窗按时间顺序保存全部结果
#[derive(Debug, Clone)]
struct Snapshot {
    label: String,
    /// 等价 CLI（子命令 + 文件）
    cli: String,
    time: String,
    text: String,
}

enum Overlay {
    Help,
    Log,
    Error(String),
    Browser(Browser),
    Viewer(ViewerState),
    /// F3：历次结果弹窗，index 指向 history（最新的在末尾）
    Results {
        index: usize,
        scroll: u16,
    },
}

// ===== 应用状态 =====

thread_local! {
    /// worker 线程里的 panic 已被 catch_unwind 捕获，不应破坏全屏界面
    static QUIET_PANIC: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

struct App {
    log_lines: Arc<Mutex<VecDeque<String>>>,
    current_file: Option<PathBuf>,
    header: Option<Header>,
    header_rx: Option<Receiver<Option<Header>>>,
    mode: Mode,
    focus: Focus,
    files: FilePanel,
    expanded: HashSet<usize>,
    menu_index: usize,
    tree_state: ListState,
    form: Option<Form>,
    result: Option<Report>,
    /// 历次结果（F3 弹窗用），最新的在末尾
    history: Vec<Snapshot>,
    result_scroll: u16,
    overlay: Option<Overlay>,
    help_scroll: u16,
    log_scroll: u16,
    job: Option<Job>,
    status: String,
    quit: bool,
    tick: usize,
}

impl App {
    fn new(log_lines: Arc<Mutex<VecDeque<String>>>, initial: Option<PathBuf>) -> Self {
        let cwd = std::env::current_dir().unwrap_or_else(|_| PathBuf::from("."));
        let mut current_file = None;
        let mut files = FilePanel::new(cwd.clone());

        if let Some(path) = initial {
            let path = std::path::absolute(&path).unwrap_or(path);
            if path.is_file() {
                current_file = Some(path.clone());
                let parent = path.parent().unwrap_or(&cwd).to_path_buf();
                files = FilePanel::new(parent);
                files.select_path(&path);
            } else if path.is_dir() {
                files = FilePanel::new(path);
            }
        }

        // 默认折叠，先展示 CLI 的顶层结构；用 ←/→ 展开需要的分组
        let expanded = HashSet::new();

        let mut app = Self {
            log_lines,
            current_file,
            header: None,
            header_rx: None,
            mode: Mode::Welcome,
            focus: Focus::Files,
            files,
            expanded,
            menu_index: 0,
            tree_state: ListState::default(),
            form: None,
            result: None,
            history: Vec::new(),
            result_scroll: 0,
            overlay: None,
            help_scroll: 0,
            log_scroll: 0,
            job: None,
            status: "Tab 切换 文件/分析菜单，Enter 选择；o 打开文件浏览器".to_string(),
            quit: false,
            tick: 0,
        };

        if app.current_file.is_some() {
            app.focus = Focus::Tree;
            if let Some(path) = app.current_file.clone() {
                app.spawn_header(path);
            }
        }
        app
    }

    // ---------- 事件循环 ----------

    fn event_loop(&mut self, terminal: &mut Terminal<CrosstermBackend<io::Stdout>>) -> io::Result<()> {
        while !self.quit {
            self.poll();
            terminal.draw(|frame| self.draw(frame))?;
            if event::poll(Duration::from_millis(100))? {
                match event::read()? {
                    Event::Key(key) if matches!(key.kind, KeyEventKind::Press | KeyEventKind::Repeat) => self.on_key(key),
                    Event::Paste(text) => self.on_paste(text),
                    _ => {}
                }
            }
            self.tick = self.tick.wrapping_add(1);
        }
        Ok(())
    }

    fn poll(&mut self) {
        self.poll_header();
        self.poll_job();
    }

    fn poll_header(&mut self) {
        if let Some(rx) = self.header_rx.as_ref() {
            match rx.try_recv() {
                Ok(header) => {
                    self.header = header;
                    self.header_rx = None;
                }
                Err(TryRecvError::Empty) => {}
                Err(TryRecvError::Disconnected) => self.header_rx = None,
            }
        }
    }

    fn poll_job(&mut self) {
        let mut outcome = None;
        if let Some(job) = self.job.as_ref() {
            match job.rx.try_recv() {
                Ok(value) => outcome = Some(value),
                Err(TryRecvError::Empty) => {}
                Err(TryRecvError::Disconnected) => outcome = Some(JobOutcome::Err("任务线程异常退出".to_string())),
            }
        }
        let Some(outcome) = outcome else { return };
        let job = self.job.take().expect("job 存在");
        match outcome {
            JobOutcome::Ok(report) => {
                self.result_scroll = 0;
                let cli = self.form.as_ref().map(|form| task_info(form.task).cli).unwrap_or("");
                let file = self.current_file.as_ref().map(|path| path.display().to_string()).unwrap_or_default();
                self.history.push(Snapshot {
                    label: job.label.clone(),
                    cli: format!("{cli} {file}").trim().to_string(),
                    time: chrono::Local::now().format("%H:%M:%S").to_string(),
                    text: report.to_string(),
                });
                self.result = Some(report);
                if job.cancelled {
                    self.status = format!("任务完成（后台）: {}，按 F3 查看", job.label);
                } else {
                    self.mode = Mode::Result;
                    self.status = format!("完成: {}", job.label);
                }
            }
            JobOutcome::Err(err) => {
                if job.cancelled {
                    self.status = format!("任务失败（后台）: {err}");
                } else {
                    self.mode = Mode::Form;
                    self.overlay = Some(Overlay::Error(err));
                }
            }
        }
    }

    // ---------- 输入 ----------

    fn on_key(&mut self, key: KeyEvent) {
        if key.modifiers.contains(KeyModifiers::CONTROL) && key.code == KeyCode::Char('c') {
            self.quit = true;
            return;
        }
        // F1/F2/F3 全局可用：弹层里也能直接切过去（F3 会先关掉弹层）
        match key.code {
            KeyCode::F(1) => {
                self.help_scroll = 0;
                self.overlay = Some(Overlay::Help);
                return;
            }
            KeyCode::F(2) => {
                self.log_scroll = 0;
                self.overlay = Some(Overlay::Log);
                return;
            }
            KeyCode::F(3) => {
                // F3 是独立弹窗：任意时刻都能打开，里面是历次全部结果
                if self.history.is_empty() {
                    self.status = "还没有计算结果：在功能层选一个任务并按 Enter".to_string();
                } else {
                    let last = self.history.len() - 1;
                    self.overlay = Some(Overlay::Results { index: last, scroll: 0 });
                }
                return;
            }
            _ => {}
        }
        if self.overlay.is_some() {
            self.on_key_overlay(key);
            return;
        }
        if key.code == KeyCode::F(4) {
            self.open_viewer();
            return;
        }
        match self.mode {
            Mode::Welcome => self.on_key_menu(key),
            Mode::Form => self.on_key_form(key),
            Mode::Running => self.on_key_running(key),
            Mode::Result => self.on_key_result(key),
        }
    }

    fn on_key_overlay(&mut self, key: KeyEvent) {
        if matches!(self.overlay, Some(Overlay::Viewer(_))) {
            self.on_key_viewer(key);
            return;
        }
        if matches!(self.overlay, Some(Overlay::Results { .. })) {
            self.on_key_results(key);
            return;
        }
        let is_browser = matches!(self.overlay.as_ref(), Some(Overlay::Browser(_)));
        if is_browser {
            let event = {
                let Some(Overlay::Browser(browser)) = self.overlay.as_mut() else { return };
                browser.on_key(key)
            };
            match event {
                BrowserEvent::Stay => {}
                BrowserEvent::Close => self.overlay = None,
                BrowserEvent::Choose(path) => self.choose_path(path),
            }
            return;
        }

        let help = matches!(self.overlay.as_ref(), Some(Overlay::Help));
        let log = matches!(self.overlay.as_ref(), Some(Overlay::Log));
        match key.code {
            KeyCode::Esc | KeyCode::Enter | KeyCode::Char('q') => self.overlay = None,
            KeyCode::Up | KeyCode::Char('k') => {
                if help {
                    self.help_scroll = self.help_scroll.saturating_sub(1);
                }
                if log {
                    self.log_scroll = self.log_scroll.saturating_sub(1);
                }
            }
            KeyCode::Down | KeyCode::Char('j') => {
                if help {
                    self.help_scroll = self.help_scroll.saturating_add(1);
                }
                if log {
                    self.log_scroll = self.log_scroll.saturating_add(1);
                }
            }
            KeyCode::PageUp => {
                if help {
                    self.help_scroll = self.help_scroll.saturating_sub(10);
                }
                if log {
                    self.log_scroll = self.log_scroll.saturating_sub(10);
                }
            }
            KeyCode::PageDown => {
                if help {
                    self.help_scroll = self.help_scroll.saturating_add(10);
                }
                if log {
                    self.log_scroll = self.log_scroll.saturating_add(10);
                }
            }
            _ => {}
        }
    }

    /// 结构预览窗口的按键：旋转、缩放、切换模型与显示开关
    fn on_key_viewer(&mut self, key: KeyEvent) {
        if matches!(key.code, KeyCode::Esc | KeyCode::Enter | KeyCode::Char('q') | KeyCode::F(4)) {
            self.overlay = None;
            return;
        }
        let pi = std::f64::consts::FRAC_PI_2;
        let Some(Overlay::Viewer(state)) = self.overlay.as_mut() else { return };
        match key.code {
            KeyCode::Left | KeyCode::Char('h') => state.view.rotate(-0.12, 0.0),
            KeyCode::Right | KeyCode::Char('l') => state.view.rotate(0.12, 0.0),
            KeyCode::Up | KeyCode::Char('k') => state.view.rotate(0.0, 0.12),
            KeyCode::Down | KeyCode::Char('j') => state.view.rotate(0.0, -0.12),
            KeyCode::Char('+') | KeyCode::Char('=') => state.view.zoom_by(1.12),
            KeyCode::Char('-') | KeyCode::Char('_') => state.view.zoom_by(1.0 / 1.12),
            KeyCode::Char('r') => state.view = View::default(),
            KeyCode::Char('s') => state.view.next_style(),
            KeyCode::Char('t') => state.view.labels = !state.view.labels,
            KeyCode::Char('d') => state.view.depth_cue = !state.view.depth_cue,
            KeyCode::Char('c') => state.depth = next_depth(state.depth),
            KeyCode::Char('1') => state.view.angles(0.0, 0.0),
            KeyCode::Char('2') => state.view.angles(pi, 0.0),
            KeyCode::Char('3') => state.view.angles(0.0, pi),
            _ => {}
        }
    }

    fn active_panel(&self) -> Panel {
        match self.mode {
            Mode::Welcome => {
                if self.focus == Focus::Files {
                    Panel::Files
                } else {
                    Panel::Menu
                }
            }
            // 参数区提交后焦点直接交给结果区
            Mode::Form => Panel::Params,
            Mode::Running | Mode::Result => Panel::Result,
        }
    }

    /// `+` 展开全部分组，`-` 收起全部分组
    fn expand_all(&mut self, expand: bool) {
        self.expanded.clear();
        if expand {
            fn walk(nodes: &'static [Node], out: &mut HashSet<usize>) {
                for node in nodes {
                    if !node.children.is_empty() {
                        out.insert(node_key(node));
                        walk(node.children, out);
                    }
                }
            }
            walk(TREE, &mut self.expanded);
        }
        self.menu_index = 0;
    }

    fn on_key_menu(&mut self, key: KeyEvent) {
        match key.code {
            KeyCode::Char('q') => self.quit = true,
            KeyCode::Esc => match self.focus {
                // 功能层 → 回上层（文件选择）
                Focus::Tree => {
                    self.focus = Focus::Files;
                    self.status = "已返回上一层：文件选择".to_string();
                }
                // 最外层：不退出（防误触）
                Focus::Files => self.status = "已在最外层（文件选择）；q 或 Ctrl+C 退出".to_string(),
            },
            KeyCode::Char('+') | KeyCode::Char('=') => {
                self.expand_all(true);
                self.status = "已展开全部分组".to_string();
            }
            KeyCode::Char('-') | KeyCode::Char('_') => {
                self.expand_all(false);
                self.status = "已收起全部分组".to_string();
            }
            KeyCode::Tab | KeyCode::BackTab => {
                self.focus = match self.focus {
                    Focus::Files => Focus::Tree,
                    Focus::Tree => Focus::Files,
                };
            }
            KeyCode::Up | KeyCode::Char('k') => self.move_selection(-1),
            KeyCode::Down | KeyCode::Char('j') => self.move_selection(1),
            KeyCode::Left | KeyCode::Char('h') => self.collapse_or_parent(),
            KeyCode::Right | KeyCode::Char('l') => self.expand_or_child(),
            KeyCode::Enter | KeyCode::Char(' ') => self.activate_selection(),
            KeyCode::Char('o') => self.open_browser(BrowserTarget::CurrentFile),
            KeyCode::Char('p') => self.open_viewer(),
            KeyCode::Char('v') => {
                if self.result.is_some() {
                    self.mode = Mode::Result;
                } else {
                    self.status = "还没有计算结果".to_string();
                }
            }
            _ => {}
        }
    }

    fn on_key_form(&mut self, key: KeyEvent) {
        if self.form.is_none() {
            self.mode = Mode::Welcome;
            return;
        }
        if key.modifiers.contains(KeyModifiers::CONTROL) {
            if key.code == KeyCode::Char('o') {
                let selected = self.form.as_ref().map(|form| form.selected).unwrap_or(0);
                let editable = self.form.as_ref().and_then(Form::selected_field).map(Field::editable).unwrap_or(false);
                if editable {
                    self.open_browser(BrowserTarget::Field(selected));
                } else {
                    self.status = "当前字段不能选择文件".to_string();
                }
            }
            return;
        }
        if key.modifiers.contains(KeyModifiers::ALT) {
            return;
        }
        match key.code {
            KeyCode::Esc => {
                self.mode = Mode::Welcome;
                self.focus = Focus::Tree;
            }
            KeyCode::Up | KeyCode::BackTab => {
                if let Some(form) = self.form.as_mut() {
                    form.prev();
                }
            }
            KeyCode::Down | KeyCode::Tab => {
                if let Some(form) = self.form.as_mut() {
                    form.next();
                }
            }
            // 新计算开始前先清空，避免把上一次的结果误认为本次已完成
            KeyCode::Enter => {
                self.result = None;
                self.result_scroll = 0;
                self.start_job();
            }
            _ => {
                if let Some(form) = self.form.as_mut() {
                    form.edit(key);
                }
            }
        }
    }

    fn on_key_running(&mut self, key: KeyEvent) {
        match key.code {
            KeyCode::Esc => {
                if let Some(job) = self.job.as_mut() {
                    job.cancelled = true;
                }
                self.mode = Mode::Form;
                self.status = "计算已转入后台，完成后按 F3 查看".to_string();
            }
            KeyCode::Char('q') => self.quit = true,
            _ => {}
        }
    }

    fn on_key_result(&mut self, key: KeyEvent) {
        match key.code {
            // 层级回退：结果 → 参数 → 分析菜单
            KeyCode::Esc => {
                if self.form.is_some() {
                    self.mode = Mode::Form;
                    self.status = "已返回上一层：参数输入".to_string();
                } else {
                    self.mode = Mode::Welcome;
                    self.focus = Focus::Tree;
                    self.status = "已返回上一层：分析菜单".to_string();
                }
            }
            KeyCode::Tab | KeyCode::BackTab => self.mode = Mode::Form,
            KeyCode::Char('c') => self.copy_result(),
            KeyCode::Char('s') => self.export(false),
            KeyCode::Up => self.result_scroll = self.result_scroll.saturating_sub(1),
            KeyCode::Down => self.result_scroll = self.result_scroll.saturating_add(1),
            KeyCode::PageUp => self.result_scroll = self.result_scroll.saturating_sub(10),
            KeyCode::PageDown => self.result_scroll = self.result_scroll.saturating_add(10),
            KeyCode::Home => self.result_scroll = 0,
            KeyCode::End => self.result_scroll = u16::MAX,
            KeyCode::Char('e') => self.export(false),
            KeyCode::Char('j') => self.export(true),
            KeyCode::Char('o') => self.open_browser(BrowserTarget::CurrentFile),
            KeyCode::Char('p') => self.open_viewer(),
            KeyCode::Char('q') => self.quit = true,
            _ => {}
        }
    }

    /// F3 弹窗的按键：切结果、滚动、复制/保存
    fn on_key_results(&mut self, key: KeyEvent) {
        let Some(Overlay::Results { index, scroll }) = self.overlay.as_ref() else { return };
        let mut index = *index;
        let mut scroll = *scroll;
        let len = self.history.len();
        match key.code {
            KeyCode::Esc | KeyCode::Char('q') => {
                self.overlay = None;
                return;
            }
            KeyCode::Up | KeyCode::Char('k') => scroll = scroll.saturating_sub(1),
            KeyCode::Down | KeyCode::Char('j') => scroll = scroll.saturating_add(1),
            KeyCode::PageUp => scroll = scroll.saturating_sub(10),
            KeyCode::PageDown => scroll = scroll.saturating_add(10),
            KeyCode::Home => scroll = 0,
            KeyCode::Left | KeyCode::Char('[') => {
                if index > 0 {
                    index -= 1;
                    scroll = 0;
                }
            }
            KeyCode::Right | KeyCode::Char(']') => {
                if index + 1 < len {
                    index += 1;
                    scroll = 0;
                }
            }
            KeyCode::Char('c') => {
                if let Some(snapshot) = self.history.get(index) {
                    let text = snapshot.text.clone();
                    self.copy_text(&text);
                }
            }
            KeyCode::Char('s') => {
                if let Some(snapshot) = self.history.get(index) {
                    let text = snapshot.text.clone();
                    self.save_text(&text);
                }
            }
            _ => {}
        }
        if let Some(Overlay::Results { index: i, scroll: s }) = self.overlay.as_mut() {
            *i = index;
            *s = scroll;
        }
    }

    /// `c`：用 OSC 52 让本地终端把结果写进剪贴板，同时落盘一份做兜底
    fn copy_result(&mut self) {
        let Some(report) = self.result.as_ref() else { return };
        let text = report.to_string();
        self.copy_text(&text);
    }

    fn copy_text(&mut self, text: &str) {
        let stamp = chrono::Local::now().format("%Y%m%d_%H%M%S");
        let path = format!("rswfn_copy_{stamp}.txt");
        let saved = std::fs::write(&path, text).is_ok();

        use std::io::Write as _;
        let mut out = std::io::stdout();
        let _ = write!(out, "\x1b]52;c;{}\x07", base64(text.as_bytes()));
        let _ = out.flush();

        self.status = if saved {
            format!("已发送 OSC 52 复制（终端需支持）；同时已存 {path}")
        } else {
            "已发送 OSC 52 复制请求（终端需支持）".to_string()
        };
    }

    fn save_text(&mut self, text: &str) {
        let stamp = chrono::Local::now().format("%Y%m%d_%H%M%S");
        let path = format!("rswfn_result_{stamp}.txt");
        self.status = match std::fs::write(&path, text) {
            Ok(()) => format!("已保存: {path}"),
            Err(err) => format!("保存失败: {err}"),
        };
    }

    fn on_paste(&mut self, text: String) {
        let text = text.trim().trim_matches('"').to_string();
        if text.is_empty() {
            return;
        }

        let mut choose = None;
        if let Some(Overlay::Browser(browser)) = self.overlay.as_mut() {
            let path = PathBuf::from(&text);
            if path.is_dir() {
                browser.panel.enter_dir(path);
            } else if path.is_file() {
                choose = Some(path);
            }
        }
        if let Some(path) = choose {
            self.choose_path(path);
            return;
        }

        if self.mode == Mode::Form {
            if let Some(form) = self.form.as_mut() {
                let selected = form.selected;
                if let Some(Field::Text { input, .. } | Field::File { input, .. }) = form.fields.get_mut(selected) {
                    for ch in text.chars() {
                        input.insert(ch);
                    }
                }
            }
        }
    }

    fn move_selection(&mut self, delta: isize) {
        match self.focus {
            Focus::Files => self.files.move_by(delta),
            Focus::Tree => {
                let len = self.tree_rows().len();
                if len == 0 {
                    self.menu_index = 0;
                    return;
                }
                self.menu_index = (self.menu_index as isize + delta).clamp(0, len as isize - 1) as usize;
            }
        }
    }

    fn activate_selection(&mut self) {
        match self.focus {
            Focus::Files => {
                if let Some(path) = self.files.enter_selected() {
                    self.set_current_file(path);
                    self.focus = Focus::Tree;
                }
            }
            Focus::Tree => {
                let rows = self.tree_rows();
                let Some(row) = rows.get(self.menu_index) else { return };
                if row.is_group {
                    let key = node_key(row.node);
                    if !self.expanded.remove(&key) {
                        self.expanded.insert(key);
                    }
                    self.clamp_menu_index();
                } else if matches!(row.node.kind, NodeKind::Viewer) {
                    self.open_viewer();
                } else if let Some(task) = node_task(row.node) {
                    let info = task_info(task);
                    self.form = Some(Form::new(task, self.current_file.as_deref()));
                    self.mode = Mode::Form;
                    // 换了任务就把上一次的结果清掉，免得误以为是这次的结果
                    self.result = None;
                    self.result_scroll = 0;
                    self.status = format!("rswfn {} …", info.cli);
                }
            }
        }
    }

    fn expand_or_child(&mut self) {
        if self.focus == Focus::Files {
            if let Some(entry) = self.files.current().cloned() {
                if entry.is_dir {
                    self.files.enter_dir(entry.path);
                }
            }
            return;
        }
        let rows = self.tree_rows();
        let Some(row) = rows.get(self.menu_index) else { return };
        if row.is_group {
            let key = node_key(row.node);
            if !self.expanded.contains(&key) {
                self.expanded.insert(key);
            } else if self.menu_index + 1 < rows.len() {
                self.menu_index += 1;
            }
        }
    }

    fn collapse_or_parent(&mut self) {
        if self.focus == Focus::Files {
            self.files.go_parent();
            return;
        }
        let rows = self.tree_rows();
        let Some(row) = rows.get(self.menu_index) else { return };
        if row.is_group && row.expanded {
            self.expanded.remove(&node_key(row.node));
            self.clamp_menu_index();
        } else if let Some(parent) = row.parent {
            self.menu_index = parent;
        }
    }

    fn tree_rows(&self) -> Vec<TreeRow> {
        let mut rows = Vec::new();
        flatten(TREE, &self.expanded, 0, None, &mut rows);
        rows
    }

    fn clamp_menu_index(&mut self) {
        let len = self.tree_rows().len();
        self.menu_index = self.menu_index.min(len.saturating_sub(1));
    }

    // ---------- 动作 ----------

    fn start_job(&mut self) {
        if self.job.is_some() {
            self.status = "已有任务在运行，请稍候".to_string();
            return;
        }
        let Some(form) = self.form.clone() else { return };
        let task = match build_task(&form, self.current_file.as_deref()) {
            Ok(task) => task,
            Err(err) => {
                self.overlay = Some(Overlay::Error(err));
                return;
            }
        };
        let label = form.title.to_string();
        let (tx, rx) = mpsc::channel();
        thread::spawn(move || {
            QUIET_PANIC.with(|flag| flag.set(true));
            let result = catch_unwind(AssertUnwindSafe(task));
            QUIET_PANIC.with(|flag| flag.set(false));
            let outcome = match result {
                Ok(Ok(report)) => JobOutcome::Ok(report),
                Ok(Err(err)) => JobOutcome::Err(err),
                Err(payload) => JobOutcome::Err(format!("内部错误: {}", panic_message(payload))),
            };
            let _ = tx.send(outcome);
        });
        self.job = Some(Job {
            label,
            started: Instant::now(),
            rx,
            cancelled: false,
        });
        self.mode = Mode::Running;
        self.status.clear();
    }

    fn open_browser(&mut self, target: BrowserTarget) {
        let mut start = self.files.dir.clone();

        if let BrowserTarget::Field(i) = target {
            if let Some(form) = self.form.as_ref() {
                if let Some(field) = form.fields.get(i) {
                    let text = field.text().trim();
                    if !text.is_empty() {
                        if let Some(parent) = Path::new(text).parent() {
                            if parent.is_dir() {
                                start = parent.to_path_buf();
                            }
                        }
                    }
                }
            }
        }
        self.overlay = Some(Overlay::Browser(Browser::new(start, target)));
    }

    fn choose_path(&mut self, path: PathBuf) {
        let target = match self.overlay.as_ref() {
            Some(Overlay::Browser(browser)) => browser.target,
            _ => BrowserTarget::None,
        };
        self.overlay = None;
        match target {
            BrowserTarget::CurrentFile => {
                self.set_current_file(path);
                self.mode = Mode::Welcome;
            }
            BrowserTarget::Field(i) => {
                if let Some(form) = self.form.as_mut() {
                    if let Some(Field::Text { input, .. } | Field::File { input, .. }) = form.fields.get_mut(i) {
                        *input = TextInput::new(&path.to_string_lossy());
                    }
                }
                self.mode = Mode::Form;
            }
            BrowserTarget::None => {}
        }
    }

    /// 打开分子结构预览（F4 / p）
    fn open_viewer(&mut self) {
        // 文件面板里 F4 预览光标所在的分子（不一定是已打开的那个）
        if self.active_panel() == Panel::Files {
            if let Some(entry) = self.files.current().filter(|entry| !entry.is_dir).cloned() {
                let name = entry.name.clone();
                match viewer::load(&entry.path) {
                    Ok(mol) if !mol.is_empty() => {
                        self.overlay = Some(Overlay::Viewer(ViewerState {
                            title: name,
                            mol,
                            view: View::default(),
                            depth: viewer::Depth::detect(),
                        }));
                    }
                    Ok(_) => self.status = format!("{name} 里没有可预览的原子坐标"),
                    Err(err) => self.status = format!("无法预览 {name}: {err}"),
                }
                return;
            }
        }
        let Some(header) = self.header.as_ref() else {
            self.status = if self.header_rx.is_some() {
                "正在读取分子信息，请稍候…".to_string()
            } else {
                "请先打开分子文件（o 选择文件）".to_string()
            };
            return;
        };
        if header.mol.is_empty() {
            self.status = "当前文件没有可预览的原子坐标".to_string();
            return;
        }
        let name = self
            .current_file
            .as_ref()
            .and_then(|path| path.file_name())
            .map(|name| name.to_string_lossy().to_string())
            .unwrap_or_else(|| header.formula.clone());
        let title = format!("{} · {} 个原子", name, header.natoms);
        let mol = header.mol.clone();
        self.overlay = Some(Overlay::Viewer(ViewerState {
            title,
            mol,
            view: View::default(),
            depth: viewer::Depth::detect(),
        }));
    }

    fn set_current_file(&mut self, path: PathBuf) {
        self.current_file = Some(path.clone());
        self.header = None;
        if let Some(parent) = path.parent() {
            if parent != self.files.dir.as_path() {
                self.files.enter_dir(parent.to_path_buf());
            }
            self.files.select_path(&path);
        }
        self.status = format!("已选择: {}", path.display());
        self.spawn_header(path);
    }

    fn spawn_header(&mut self, path: PathBuf) {
        let (tx, rx) = mpsc::channel();
        thread::spawn(move || {
            QUIET_PANIC.with(|flag| flag.set(true));
            let header = catch_unwind(AssertUnwindSafe(|| Header::from_path(&path))).unwrap_or(None);
            QUIET_PANIC.with(|flag| flag.set(false));
            let _ = tx.send(header);
        });
        self.header_rx = Some(rx);
    }

    fn export(&mut self, json: bool) {
        let Some(report) = self.result.as_ref() else { return };
        let name = sanitize(&report.title);
        let stamp = chrono::Local::now().format("%Y%m%d_%H%M%S");
        let filename = if json { format!("rswfn_{name}_{stamp}.json") } else { format!("rswfn_{name}_{stamp}.txt") };
        let content = if json { report.to_json() } else { report.to_string() };
        match std::fs::write(&filename, content) {
            Ok(()) => self.status = format!("已导出: {filename}"),
            Err(err) => self.status = format!("导出失败: {err}"),
        }
    }

    fn spinner(&self) -> &'static str {
        ["|", "/", "-", "\\"][self.tick % 4]
    }

    // ---------- 绘制 ----------

    fn draw(&mut self, frame: &mut Frame) {
        let area = frame.area();
        let chunks = Layout::vertical([Constraint::Length(1), Constraint::Min(3), Constraint::Length(2)]).split(area);
        self.draw_header(frame, chunks[0]);

        let body = Layout::horizontal([Constraint::Length(26), Constraint::Length(34), Constraint::Min(20)]).split(chunks[1]);
        self.draw_files(frame, body[0]);
        self.draw_tree(frame, body[1]);
        match self.mode {
            Mode::Welcome => self.draw_welcome(frame, body[2]),
            // 运行层：参数在上、结果在下，中间一条横线
            Mode::Form | Mode::Running | Mode::Result => self.draw_run(frame, body[2]),
        }
        self.draw_status(frame, chunks[2]);
        // 弹层盖住内容区即可，不要压住顶部信息条和底部状态栏
        self.draw_overlay(frame, chunks[1]);
    }

    fn draw_header(&self, frame: &mut Frame, area: Rect) {
        let file = self.current_file.as_ref().map(|path| path.display().to_string()).unwrap_or_else(|| "未打开文件（按 o 选择）".to_string());
        let mut spans = vec![
            Span::styled(" rswfn ", Style::default().fg(Color::Black).bg(Color::Cyan).add_modifier(Modifier::BOLD)),
            Span::raw(" "),
            Span::styled(file, Style::default().fg(Color::White)),
        ];
        if let Some(header) = &self.header {
            spans.push(Span::styled(
                format!("   分子式 {}   电荷 {}   多重度 {}   原子数 {}", header.formula, header.charge, header.multip, header.natoms),
                Style::default().fg(Color::Gray),
            ));
        } else if self.header_rx.is_some() {
            spans.push(Span::styled("   读取分子信息…", Style::default().fg(Color::DarkGray)));
        }
        frame.render_widget(Paragraph::new(Line::from(spans)), area);
    }

    fn draw_files(&mut self, frame: &mut Frame, area: Rect) {
        let items = self.files.items(self.current_file.as_deref());
        let border = panel_border(self.active_panel() == Panel::Files);
        let title = format!(" 文件 · {} ", self.files.dir_title());
        let list = List::new(items)
            .block(Block::new().borders(Borders::ALL).border_style(border).title(title))
            .highlight_style(Style::default().bg(Color::DarkGray).fg(Color::White))
            .highlight_symbol("▸");
        self.files.state.select(Some(self.files.selected));
        frame.render_stateful_widget(list, area, &mut self.files.state);
    }

    fn draw_tree(&mut self, frame: &mut Frame, area: Rect) {
        let rows = self.tree_rows();
        let mut items: Vec<ListItem> = Vec::new();
        for row in &rows {
            let indent = "  ".repeat(row.depth);
            let marker = if row.is_group { if row.expanded { "▼ " } else { "▶ " } } else { "· " };
            let style = if row.is_group { Style::default().fg(Color::Cyan) } else { Style::default() };
            items.push(ListItem::new(Line::from(Span::styled(format!("{indent}{marker}{}", node_label(row.node)), style))));
        }
        let border = panel_border(self.active_panel() == Panel::Menu);
        let list = List::new(items)
            .block(Block::new().borders(Borders::ALL).border_style(border).title(" 分析菜单 "))
            .highlight_style(Style::default().bg(Color::DarkGray).fg(Color::White).add_modifier(Modifier::BOLD))
            .highlight_symbol("▸");
        let selected = self.menu_index.min(rows.len().saturating_sub(1));
        self.tree_state.select(Some(selected));
        frame.render_stateful_widget(list, area, &mut self.tree_state);
    }

    fn draw_welcome(&self, frame: &mut Frame, area: Rect) {
        let lines = vec![
            Line::from(Span::styled(" 使用说明", Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD))),
            Line::from(""),
            Line::from(" 1. 左侧文件列表选择文件；也可用 rswfn <文件或目录> 直接启动"),
            Line::from(" 2. Tab 在 文件列表 / 分析菜单 之间切换"),
            Line::from(" 3. 分析菜单中 ←/→ 折叠/展开，Enter 打开参数表单"),
            Line::from(" 4. 表单中 Enter 运行，Ctrl+O 选择文件；结果页 e/j 导出"),
            Line::from(""),
            Line::from(Span::styled(" F1 帮助    F2 运行日志    F3 结果    q / Ctrl+C 退出（Esc 不会退出）", Style::default().fg(Color::Gray))),
        ];
        let paragraph = Paragraph::new(lines).block(Block::new().borders(Borders::ALL).title(" 欢迎 "));
        frame.render_widget(paragraph, area);
    }

    /// 运行层：参数在上、结果在下，中间只隔一条横线
    fn draw_run(&mut self, frame: &mut Frame, area: Rect) {
        // 参数区按字段数量定高（标题 + 字段，最多 8 行），剩下的都给结果区
        let field_rows = self.form.as_ref().map(|form| form.fields.len()).unwrap_or(0) as u16;
        let param_height = (field_rows + 1).clamp(2, 8);
        let chunks = Layout::vertical([Constraint::Length(param_height), Constraint::Length(1), Constraint::Min(3)]).split(area);
        let active = self.active_panel();
        self.draw_params(frame, chunks[0], active == Panel::Params);
        let rule = "─".repeat(chunks[1].width as usize);
        frame.render_widget(Paragraph::new(Line::from(Span::styled(rule, Style::default().fg(Color::DarkGray)))), chunks[1]);
        if self.mode == Mode::Running {
            self.draw_running_plain(frame, chunks[2]);
        } else {
            self.draw_result_plain(frame, chunks[2], active == Panel::Result);
        }
    }

    /// 参数区：无边框；焦点用左侧 ▸ 与颜色区分
    fn draw_params(&self, frame: &mut Frame, area: Rect, focused: bool) {
        let Some(form) = &self.form else { return };
        let parts = Layout::vertical([Constraint::Length(1), Constraint::Min(1)]).split(area);
        let accent = if focused { Color::Yellow } else { Color::DarkGray };
        let mut head = vec![
            Span::styled(format!("{} 参数 · {}", if focused { "▸" } else { " " }, form.title), Style::default().fg(accent).add_modifier(Modifier::BOLD)),
            Span::styled(format!("   {}", form.desc), Style::default().fg(Color::Gray)),
        ];
        // 只有文件参数（或完全没有参数）时，提示直接回车即可
        if form.fields.iter().all(|field| matches!(field, Field::File { .. })) {
            head.push(Span::styled("    不需要输入参数，按 Enter 执行计算", Style::default().fg(Color::Yellow)));
        }
        frame.render_widget(Paragraph::new(Line::from(head)), parts[0]);

        let rows: Vec<Row> = form
            .fields
            .iter()
            .enumerate()
            .map(|(i, field)| {
                let active = focused && i == form.selected;
                let label = if active { format!("▸ {}", field.label()) } else { format!("  {}", field.label()) };
                Row::new(vec![Cell::from(label), Cell::from(field.value_line(active))])
            })
            .collect();
        let table = Table::new(rows, [Constraint::Length(22), Constraint::Min(20)])
            .row_highlight_style(Style::default().bg(Color::DarkGray))
            .highlight_symbol("");
        let mut state = TableState::default();
        if focused {
            state.select(Some(form.selected));
        }
        frame.render_stateful_widget(table, parts[1], &mut state);
    }

    /// 结果区：无边框，首行给等价 CLI 命令
    fn draw_result_plain(&mut self, frame: &mut Frame, area: Rect, focused: bool) {
        let parts = Layout::vertical([Constraint::Length(1), Constraint::Min(1)]).split(area);
        let accent = if focused { Color::Yellow } else { Color::DarkGray };
        let cli = self.form.as_ref().map(|form| task_info(form.task).cli).unwrap_or("");
        let file = self.current_file.as_ref().map(|path| path.display().to_string()).unwrap_or_default();
        let head = Line::from(vec![
            Span::styled(format!("{} 结果", if focused { "▸" } else { " " }), Style::default().fg(accent).add_modifier(Modifier::BOLD)),
            Span::styled(format!("   等价命令: rswfn {cli} {file}"), Style::default().fg(Color::Gray)),
            Span::styled("    c 复制 · s 保存 · Tab 回到参数", Style::default().fg(Color::DarkGray)),
        ]);
        frame.render_widget(Paragraph::new(head), parts[0]);

        let Some(report) = self.result.as_ref() else {
            frame.render_widget(Paragraph::new(Span::styled(" 尚无结果：在参数区按 Enter 开始计算", Style::default().fg(Color::DarkGray))), parts[1]);
            return;
        };
        let lines: Vec<Line> = report
            .to_string()
            .lines()
            .map(|line| {
                let style = if line.starts_with("===") || line.starts_with("---") {
                    Style::default().fg(Color::DarkGray)
                } else {
                    Style::default()
                };
                Line::from(Span::styled(line.to_string(), style))
            })
            .collect();
        let inner = parts[1].height as usize;
        let max_scroll = lines.len().saturating_sub(inner) as u16;
        self.result_scroll = self.result_scroll.min(max_scroll);
        frame.render_widget(Paragraph::new(lines).scroll((self.result_scroll, 0)), parts[1]);
    }

    /// 运行中：同样无边框，保留上面参数区的可见性
    fn draw_running_plain(&self, frame: &mut Frame, area: Rect) {
        let job = self.job.as_ref();
        let last_log = self.log_lines.lock().ok().and_then(|lines| lines.back().cloned()).unwrap_or_default();
        let lines = vec![
            Line::from(Span::styled(
                format!("{} 正在计算（Esc 转入后台，F2 看日志）", self.spinner()),
                Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD),
            )),
            Line::from(job.map(|job| format!("任务: {}   已用时: {:.1} s", job.label, job.started.elapsed().as_secs_f64())).unwrap_or_default()),
            Line::from(Span::styled(last_log, Style::default().fg(Color::DarkGray))),
        ];
        frame.render_widget(Paragraph::new(lines), area);
    }

    fn draw_form(&self, frame: &mut Frame, area: Rect) {
        let Some(form) = &self.form else { return };
        let chunks = Layout::vertical([Constraint::Min(3), Constraint::Length(1)]).split(area);

        let rows: Vec<Row> = form
            .fields
            .iter()
            .enumerate()
            .map(|(i, field)| {
                let active = i == form.selected;
                let label = if active { format!("▸ {}", field.label()) } else { format!("  {}", field.label()) };
                Row::new(vec![Cell::from(label), Cell::from(field.value_line(active))])
            })
            .collect();

        let table = Table::new(rows, [Constraint::Length(22), Constraint::Min(20)])
            .block(Block::new().borders(Borders::ALL).title(format!(" {} · {} ", form.title, form.desc)))
            .row_highlight_style(Style::default().bg(Color::DarkGray))
            .highlight_symbol("");
        let mut state = TableState::default();
        state.select(Some(form.selected));
        frame.render_stateful_widget(table, chunks[0], &mut state);

        let hint = format!(" {}   |   rswfn {} <file>   |   F1 帮助  F2 日志", form.hint(), task_info(form.task).cli);
        frame.render_widget(Paragraph::new(Line::from(Span::styled(hint, Style::default().fg(Color::Gray)))), chunks[1]);
    }

    fn draw_running(&self, frame: &mut Frame, area: Rect) {
        let job = self.job.as_ref();
        let last_log = self.log_lines.lock().ok().and_then(|lines| lines.back().cloned()).unwrap_or_default();
        let lines = vec![
            Line::from(""),
            Line::from(Span::styled(format!("{} 正在计算", self.spinner()), Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD))),
            Line::from(""),
            Line::from(job.map(|job| format!("任务: {}", job.label)).unwrap_or_default()),
            Line::from(job.map(|job| format!("已用时: {:.1} s", job.started.elapsed().as_secs_f64())).unwrap_or_default()),
            Line::from(""),
            Line::from(Span::styled(last_log, Style::default().fg(Color::DarkGray))),
            Line::from(""),
            Line::from(Span::styled("Esc 转入后台运行 · F2 查看日志 · Ctrl+C 退出", Style::default().fg(Color::Gray))),
        ];
        let paragraph = Paragraph::new(lines).alignment(Alignment::Center).block(Block::new().borders(Borders::ALL).title(" 计算中 "));
        frame.render_widget(paragraph, area);
    }

    fn draw_result(&mut self, frame: &mut Frame, area: Rect) {
        // 标题里直接给出等价 CLI 命令，方便复制到脚本里重跑
        let cli = self.form.as_ref().map(|form| task_info(form.task).cli).unwrap_or("");
        let file = self.current_file.as_ref().map(|path| path.display().to_string()).unwrap_or_default();
        let heading = if cli.is_empty() { " 结果 ".to_string() } else { format!(" 结果 · rswfn {cli} {file} ") };

        let Some(report) = self.result.as_ref() else {
            let hint = vec![
                Line::from(Span::styled(" 尚无结果", Style::default().fg(Color::DarkGray))),
                Line::from(""),
                Line::from(" 在参数区填好后按 Enter 开始计算；结果会一直保留在这里（F3 也随时可看）"),
                Line::from(Span::styled(" c 复制结果 · s 保存文本", Style::default().fg(Color::Gray))),
            ];
            frame.render_widget(Paragraph::new(hint).block(Block::new().borders(Borders::ALL).title(heading)), area);
            return;
        };
        let title = report.title.clone();
        let _ = title;
        let lines: Vec<Line> = report
            .to_string()
            .lines()
            .map(|line| {
                let style = if line.starts_with("===") || line.starts_with("---") {
                    Style::default().fg(Color::DarkGray)
                } else {
                    Style::default()
                };
                Line::from(Span::styled(line.to_string(), style))
            })
            .collect();

        let inner_height = area.height.saturating_sub(2) as usize;
        let max_scroll = lines.len().saturating_sub(inner_height) as u16;
        let scroll = self.result_scroll.min(max_scroll);
        self.result_scroll = scroll;

        let paragraph = Paragraph::new(lines)
            .scroll((scroll, 0))
            .block(Block::new().borders(Borders::ALL).title(format!("{heading} c 复制 · s 保存 · ↑↓ 滚动 ")));
        frame.render_widget(paragraph, area);
    }

    fn draw_status(&self, frame: &mut Frame, area: Rect) {
        let mode = match self.mode {
            Mode::Welcome => "菜单",
            Mode::Form => "表单",
            Mode::Running => "计算",
            Mode::Result => "结果",
        };
        let mut line1 = format!(" [{mode}] ");
        if let Some(job) = &self.job {
            line1.push_str(&format!("{} {} {:.1}s ", self.spinner(), job.label, job.started.elapsed().as_secs_f64()));
        }
        if !self.status.is_empty() {
            line1.push_str(&format!("{} ", self.status));
        }
        if self.mode == Mode::Welcome && self.focus == Focus::Tree {
            let rows = self.tree_rows();
            if let Some(row) = rows.get(self.menu_index) {
                if let Some(task) = node_task(row.node) {
                    line1.push_str(&format!("  rswfn {} <file>", task_info(task).cli));
                } else if matches!(row.node.kind, NodeKind::Viewer) {
                    line1.push_str("  rswfn view <file>");
                }
            }
        }

        let hints = match self.mode {
            Mode::Welcome => match self.focus {
                Focus::Files => "↑/↓ 选择文件  Enter 打开/进入  Backspace 上级  Tab 分析菜单  o 浏览器  F4 结构预览  F1 帮助  q 退出（Esc 不退出）",
                Focus::Tree => "↑/↓ 移动  ←/→ 折叠/展开  Enter 打开  Tab 文件列表  o 浏览器  F4 结构预览  F1 帮助  q 退出（Esc 不退出）",
            },
            Mode::Form => "↑/↓ 切换字段  Enter 运行  Ctrl+O 选择文件  Esc 返回  F4 结构预览  F1 帮助  F3 结果",
            Mode::Running => "Esc 转后台  F2 日志  F4 结构预览  F3 结果  Ctrl+C 退出",
            Mode::Result => "↑/↓/PgUp/PgDn 滚动  Home/End  e 文本  j JSON  p 结构预览  Esc 返回  o 打开文件",
        };

        let paragraph = Paragraph::new(vec![Line::from(line1), Line::from(Span::styled(hints.to_string(), Style::default().fg(Color::DarkGray)))]);
        frame.render_widget(paragraph, area);
    }

    fn draw_overlay(&mut self, frame: &mut Frame, area: Rect) {
        let current_file = self.current_file.clone();
        let Some(overlay) = self.overlay.as_mut() else { return };
        match overlay {
            Overlay::Results { index, scroll } => {
                let popup = centered(area, 86, 84);
                frame.render_widget(Clear, popup);
                let total = self.history.len();
                let snapshot = self.history.get(*index);
                let title = match snapshot {
                    Some(snapshot) => format!(" F3 计算结果 {}/{} · {} · {} ", *index + 1, total, snapshot.time, snapshot.label),
                    None => " F3 计算结果 ".to_string(),
                };
                let block = Block::new().borders(Borders::ALL).border_style(Style::default().fg(Color::Cyan)).title(title);
                let inner = block.inner(popup);
                frame.render_widget(block, popup);
                let parts = Layout::vertical([Constraint::Length(1), Constraint::Min(1), Constraint::Length(1)]).split(inner);
                let cli = snapshot.map(|snapshot| snapshot.cli.clone()).unwrap_or_default();
                frame.render_widget(Paragraph::new(Line::from(Span::styled(format!(" 等价命令: rswfn {cli}"), Style::default().fg(Color::Gray)))), parts[0]);
                let body: Vec<Line> = snapshot.map(|snapshot| snapshot.text.clone()).unwrap_or_default().lines().map(|line| Line::from(line.to_string())).collect();
                frame.render_widget(Paragraph::new(body).scroll((*scroll, 0)), parts[1]);
                let hint = " ↑↓ 滚动   ←/→ 切换结果   c 复制   s 保存   Esc 关闭 ";
                frame.render_widget(Paragraph::new(Line::from(Span::styled(hint, Style::default().fg(Color::Gray)))), parts[2]);
            }
            Overlay::Browser(browser) => {
                let popup = centered(area, 80, 80);
                frame.render_widget(Clear, popup);
                frame.render_widget(Block::new().borders(Borders::ALL).title(" 打开文件 "), popup);
                let inner = Rect {
                    x: popup.x + 1,
                    y: popup.y + 1,
                    width: popup.width.saturating_sub(2),
                    height: popup.height.saturating_sub(2),
                };
                let chunks = Layout::vertical([Constraint::Length(1), Constraint::Min(1), Constraint::Length(1)]).split(inner);

                let path_line = Line::from(vec![Span::styled(" 路径: ", Style::default().fg(Color::Gray)), Span::raw(browser.panel.dir.display().to_string())]);
                frame.render_widget(Paragraph::new(path_line), chunks[0]);

                let items = browser.panel.items(current_file.as_deref());
                let list = List::new(items).highlight_style(Style::default().bg(Color::DarkGray).fg(Color::White)).highlight_symbol("▸");
                browser.panel.state.select(Some(browser.panel.selected));
                frame.render_stateful_widget(list, chunks[1], &mut browser.panel.state);

                let hint = browser.panel.error.clone().unwrap_or_else(|| "Enter 选择/进入   Backspace 上级   Esc 取消   也可直接粘贴路径".to_string());
                frame.render_widget(Paragraph::new(Line::from(Span::styled(hint, Style::default().fg(Color::Gray)))), chunks[2]);
            }
            Overlay::Viewer(state) => {
                let popup = area;
                frame.render_widget(Clear, popup);
                let block = Block::new()
                    .borders(Borders::ALL)
                    .border_style(Style::default().fg(Color::Cyan))
                    .style(Style::default().bg(viewer::BG.color_in(state.depth)))
                    .title(format!(" 分子结构预览 · {} · {} ", state.title, state.view.style.name()));
                let inner = block.inner(popup);
                frame.render_widget(block, popup);
                let art = Rect {
                    x: inner.x,
                    y: inner.y,
                    width: inner.width,
                    height: inner.height.saturating_sub(1),
                };
                let lines = state.view.render(&state.mol, art.width, art.height).lines_in(state.depth);
                frame.render_widget(Paragraph::new(lines), art);
                if inner.height >= 1 && inner.width > 8 {
                    let footer = Rect {
                        x: inner.x,
                        y: inner.y + inner.height - 1,
                        width: inner.width,
                        height: 1,
                    };
                    let hint = " ↑↓←→ 旋转  +/- 缩放  s 模型  t 标签  d 景深  c 配色  1/2/3 视角  r 复位  Esc 关闭";
                    let paragraph = Paragraph::new(Line::from(Span::styled(hint, Style::default().fg(Color::Gray)))).style(Style::default().bg(viewer::BG.color_in(state.depth)));
                    frame.render_widget(paragraph, footer);
                }
            }
            Overlay::Help => {
                let popup = centered(area, 70, 80);
                frame.render_widget(Clear, popup);
                let lines = help_lines();
                let inner = popup.height.saturating_sub(2) as usize;
                let max_scroll = lines.len().saturating_sub(inner) as u16;
                self.help_scroll = self.help_scroll.min(max_scroll);
                let paragraph = Paragraph::new(lines)
                    .scroll((self.help_scroll, 0))
                    .wrap(Wrap { trim: false })
                    .block(Block::new().borders(Borders::ALL).title(" 帮助（Esc 关闭） "));
                frame.render_widget(paragraph, popup);
            }
            Overlay::Log => {
                let popup = centered(area, 80, 80);
                frame.render_widget(Clear, popup);
                let lines: Vec<Line> = {
                    let guard = self.log_lines.lock().ok();
                    let mut lines: Vec<Line> = guard.map(|guard| guard.iter().map(|line| Line::from(line.clone())).collect()).unwrap_or_default();
                    if lines.is_empty() {
                        lines.push(Line::from(Span::styled("（暂无日志）", Style::default().fg(Color::DarkGray))));
                    }
                    lines
                };
                let inner = popup.height.saturating_sub(2) as usize;
                let max_scroll = lines.len().saturating_sub(inner) as u16;
                self.log_scroll = self.log_scroll.min(max_scroll);
                let paragraph = Paragraph::new(lines).scroll((self.log_scroll, 0)).block(Block::new().borders(Borders::ALL).title(" 运行日志（Esc 关闭） "));
                frame.render_widget(paragraph, popup);
            }
            Overlay::Error(message) => {
                let popup = centered(area, 70, 40);
                frame.render_widget(Clear, popup);
                let paragraph = Paragraph::new(message.clone()).wrap(Wrap { trim: true }).block(Block::new().borders(Borders::ALL).title(" 错误（Esc 关闭） "));
                frame.render_widget(paragraph, popup);
            }
        }
    }
}

fn centered(area: Rect, percent_x: u16, percent_y: u16) -> Rect {
    let vertical = Layout::vertical([Constraint::Percentage((100 - percent_y) / 2), Constraint::Percentage(percent_y), Constraint::Percentage((100 - percent_y) / 2)]).split(area);
    let horizontal = Layout::horizontal([Constraint::Percentage((100 - percent_x) / 2), Constraint::Percentage(percent_x), Constraint::Percentage((100 - percent_x) / 2)]).split(vertical[1]);
    horizontal[1]
}

fn help_lines() -> Vec<Line<'static>> {
    fn section(lines: &mut Vec<Line<'static>>, title: &str) {
        lines.push(Line::from(""));
        lines.push(Line::from(Span::styled(format!(" {title}"), Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD))));
    }

    let mut lines = Vec::new();
    section(&mut lines, "全局");
    lines.push(Line::from(" F1 / ?            帮助"));
    lines.push(Line::from(" F2                运行日志"));
    lines.push(Line::from(" F3                查看上次结果"));
    lines.push(Line::from(" F4 / p            分子结构预览"));
    lines.push(Line::from(" q / Ctrl+C        退出（Esc 不会退出）"));
    lines.push(Line::from(" Esc               关闭弹窗 / 返回上一步"));
    section(&mut lines, "文件列表 / 分析菜单");
    lines.push(Line::from(" Tab               切换 文件列表 / 分析菜单"));
    lines.push(Line::from(" ↑/↓ 或 j/k        移动"));
    lines.push(Line::from(" ←/→ 或 h/l        折叠 / 展开（分析菜单）"));
    lines.push(Line::from(" Enter / 空格      打开文件 / 进入目录 / 打开表单"));
    lines.push(Line::from(" Backspace        上级目录（文件列表）"));
    lines.push(Line::from(" o                 打开文件浏览器"));
    lines.push(Line::from(" v                 查看上次计算结果"));
    section(&mut lines, "表单");
    lines.push(Line::from(" ↑/↓ 或 Tab        切换字段"));
    lines.push(Line::from(" ←/→               切换选项 / 移动光标"));
    lines.push(Line::from(" Enter             运行任务"));
    lines.push(Line::from(" Ctrl+O            为当前字段选择文件"));
    lines.push(Line::from(" Esc               返回菜单"));
    section(&mut lines, "结构预览（F4 / p）");
    lines.push(Line::from(" ↑/↓/←/→ 或 hjkl   旋转视角"));
    lines.push(Line::from(" + / -             缩放"));
    lines.push(Line::from(" s                 切换 球棍 / 空间填充 / 线框"));
    lines.push(Line::from(" t                 元素符号开关"));
    lines.push(Line::from(" d                 景深雾化开关"));
    lines.push(Line::from(" 1 / 2 / 3         正视 / 侧视 / 俯视"));
    lines.push(Line::from(" r                 复位视角"));
    lines.push(Line::from(" Esc               关闭预览"));
    section(&mut lines, "结果");
    lines.push(Line::from(" ↑/↓ / PgUp/PgDn   滚动"));
    lines.push(Line::from(" Home / End        跳到顶部 / 底部"));
    lines.push(Line::from(" e                 导出文本报告"));
    lines.push(Line::from(" j                 导出 JSON 报告"));
    lines.push(Line::from(" Esc               返回表单"));
    lines
}

// ===== 任务构建 =====

type TaskFn = Box<dyn FnOnce() -> Result<Report, String> + Send>;

fn build_task(form: &Form, current: Option<&Path>) -> Result<TaskFn, String> {
    let current = current.map(|path| path.to_path_buf());
    let need_current = || {
        let path = current.clone().ok_or_else(|| "请先打开波函数文件（按 o 选择）".to_string())?;
        if !path.is_file() {
            return Err(format!("当前文件不存在: {}", path.display()));
        }
        Ok(path)
    };

    match form.task {
        TaskId::Mulliken => {
            let path = need_current()?;
            Ok(Box::new(move || {
                let mole = service::load_mole(&path)?;
                atomprop::mulliken(&mole)
            }))
        }
        TaskId::Lowdin => {
            let path = need_current()?;
            Ok(Box::new(move || {
                let mole = service::load_mole(&path)?;
                atomprop::lowdin(&mole)
            }))
        }
        TaskId::Hirshfeld => {
            let path = need_current()?;
            Ok(Box::new(move || {
                let mole = service::load_mole(&path)?;
                atomprop::hirshfeld(&mole)
            }))
        }
        TaskId::Spin => {
            let path = need_current()?;
            let ctype = form.choice(0).to_string();
            Ok(Box::new(move || {
                let mole = service::load_mole(&path)?;
                atomprop::spin(&mole, &ctype)
            }))
        }
        TaskId::Resp => {
            let path = need_current()?;
            let density = parse_f64(form.text(0), "MK 点密度")?;
            let one_stage = form.flag(1);
            Ok(Box::new(move || {
                let mole = service::load_mole(&path)?;
                atomprop::resp(&mole, density, one_stage)
            }))
        }
        TaskId::PiPop => {
            let path = need_current()?;
            let ctype = form.choice(0).to_string();
            let method = form.choice(1).to_string();
            Ok(Box::new(move || {
                let mole = service::load_mole(&path)?;
                atomprop::pi_pop(&mole, &ctype, &method)
            }))
        }
        TaskId::Fukui => {
            let path = need_current()?;
            let anion = req_path(form.text(0), "阴离子文件")?;
            let cation = req_path(form.text(1), "阳离子文件")?;
            let ctype = form.choice(2).to_string();
            Ok(Box::new(move || {
                let mole = service::load_mole(&path)?;
                let mole_n = service::load_mole(&anion)?;
                let mole_p = service::load_mole(&cation)?;
                atomprop::fukui(&mole, &mole_n, &mole_p, &ctype)
            }))
        }
        TaskId::FukuiPi => {
            let path = need_current()?;
            let anion = req_path(form.text(0), "阴离子文件")?;
            let cation = req_path(form.text(1), "阳离子文件")?;
            let ctype = form.choice(2).to_string();
            let method = form.choice(3).to_string();
            Ok(Box::new(move || {
                let mole = service::load_mole(&path)?;
                let mole_n = service::load_mole(&anion)?;
                let mole_p = service::load_mole(&cation)?;
                atomprop::fukui_pi(&mole, &mole_n, &mole_p, &ctype, &method)
            }))
        }
        TaskId::Freev => {
            let path = need_current()?;
            let otype = form.choice(0).to_string();
            Ok(Box::new(move || {
                let mole = service::load_mole(&path)?;
                atomprop::freev(&mole, &otype)
            }))
        }
        TaskId::FreevPi => {
            let path = need_current()?;
            let otype = form.choice(0).to_string();
            Ok(Box::new(move || {
                let mole = service::load_mole(&path)?;
                atomprop::freev_pi(&mole, &otype)
            }))
        }
        TaskId::Mayer => {
            let path = need_current()?;
            Ok(Box::new(move || {
                let mole = service::load_mole(&path)?;
                bondprop::mayer(&mole)
            }))
        }
        TaskId::Wiberg => {
            let path = need_current()?;
            Ok(Box::new(move || {
                let mole = service::load_mole(&path)?;
                bondprop::wiberg(&mole)
            }))
        }
        TaskId::OrderObt => {
            let path = need_current()?;
            let obt = parse_usize(form.text(0), "轨道序号")?;
            let otype = form.choice(1).to_string();
            Ok(Box::new(move || {
                let mole = service::load_mole(&path)?;
                bondprop::order_obt(&mole, obt, &otype)
            }))
        }
        TaskId::PiOrder => {
            let path = need_current()?;
            let method = form.choice(0).to_string();
            let otype = form.choice(1).to_string();
            Ok(Box::new(move || {
                let mole = service::load_mole(&path)?;
                bondprop::pi_order(&mole, &method, &otype)
            }))
        }
        TaskId::Decom => {
            let path = need_current()?;
            let atm1 = parse_usize(form.text(0), "原子1")?;
            let atm2 = parse_usize(form.text(1), "原子2")?;
            Ok(Box::new(move || {
                let mole = service::load_mole(&path)?;
                bondprop::decom(&mole, atm1, atm2)
            }))
        }
        TaskId::BondEnergy => {
            let path = need_current()?;
            let atm1 = parse_usize(form.text(0), "原子1")?;
            let atm2 = parse_usize(form.text(1), "原子2")?;
            Ok(Box::new(move || {
                let mole = service::load_mole(&path)?;
                bondprop::bond_energy(&mole, atm1, atm2)
            }))
        }
        TaskId::Aromacity => {
            let path = need_current()?;
            let rings = split_list(form.text(0));
            Ok(Box::new(move || {
                let mole = service::load_mole(&path)?;
                moleprop::aromacity(&mole, &rings)
            }))
        }
        TaskId::PiEnergy => {
            let path = need_current()?;
            let method = form.choice(0).to_string();
            Ok(Box::new(move || {
                let mole = service::load_mole(&path)?;
                moleprop::pi_energy(&mole, &method)
            }))
        }
        TaskId::PiElesPocv => {
            let path = need_current()?;
            Ok(Box::new(move || {
                let mole = service::load_mole(&path)?;
                orbtprop::pi_eles_pocv(&mole)
            }))
        }
        TaskId::PiElesMocv => {
            let path = need_current()?;
            Ok(Box::new(move || {
                let mole = service::load_mole(&path)?;
                orbtprop::pi_eles_mocv(&mole)
            }))
        }
        TaskId::Composition => {
            let path = need_current()?;
            let obt = parse_usize(form.text(0), "轨道序号")?;
            let ctype = form.choice(1).to_string();
            Ok(Box::new(move || {
                let mole = service::load_mole(&path)?;
                orbtprop::composition(&mole, obt, &ctype)
            }))
        }
        TaskId::GridDensity => {
            let path = need_current()?;
            let func = form.choice(0).to_string();
            let step = parse_f64(form.text(1), "步长")?;
            let output = opt_path(form.text(2));
            let mole_n = opt_file_path(form.text(3), "阴离子文件")?;
            let mole_p = opt_file_path(form.text(4), "阳离子文件")?;
            Ok(Box::new(move || {
                let mole = service::load_mole(&path)?;
                let mole_n = match &mole_n {
                    Some(path) => Some(service::load_mole(path)?),
                    None => None,
                };
                let mole_p = match &mole_p {
                    Some(path) => Some(service::load_mole(path)?),
                    None => None,
                };
                gridprop::density(&mole, output.as_deref(), step, &func, mole_n.as_ref(), mole_p.as_ref())
            }))
        }
        TaskId::GridPotential => {
            let path = need_current()?;
            let step = parse_f64(form.text(0), "步长")?;
            let output = opt_path(form.text(1));
            Ok(Box::new(move || {
                let mole = service::load_mole(&path)?;
                gridprop::potential(&mole, output.as_deref(), step)
            }))
        }
        TaskId::GridTopology => {
            let path = need_current()?;
            let func = form.choice(0).to_string();
            Ok(Box::new(move || {
                let mole = service::load_mole(&path)?;
                gridprop::topology(&mole, &func)
            }))
        }
        TaskId::QmHf => {
            let input = req_path(form.text(0), "输入文件")?;
            let basis = opt_string(form.text(1));
            let output = opt_path(form.text(2));
            Ok(Box::new(move || qmcal::hf(&input, basis.as_deref(), output.as_deref())))
        }
        TaskId::QmDft => {
            let input = req_path(form.text(0), "输入文件")?;
            let basis = opt_string(form.text(1));
            let func = form.choice(2).to_string();
            let output = opt_path(form.text(3));
            Ok(Box::new(move || qmcal::dft(&input, basis.as_deref(), &func, output.as_deref())))
        }
        TaskId::Coord => {
            let path = need_current()?;
            let unit = form.choice(0).to_string();
            let output = opt_path(form.text(1));
            Ok(Box::new(move || file::coord(&path, output.as_deref(), &unit)))
        }
        TaskId::Trans => {
            let input = req_path(form.text(0), "输入文件")?;
            let output = out_path(form.text(1), "输出文件")?;
            Ok(Box::new(move || file::trans(&input, &output)))
        }
        TaskId::Dhmo => {
            let input = req_path(form.text(0), "输入文件")?;
            let output = opt_path(form.text(1));
            Ok(Box::new(move || qmcal::dhmo(&input, output.as_deref())))
        }
        TaskId::Si => {
            let inputs: Vec<PathBuf> = split_list(form.text(0)).iter().map(|path| req_path(path, "输入文件")).collect::<Result<_, _>>()?;
            if inputs.is_empty() {
                return Err("至少需要一个输入文件".to_string());
            }
            let output = opt_path(form.text(1));
            Ok(Box::new(move || file::si(&inputs, output.as_deref())))
        }
    }
}

// ===== 小工具 =====

fn parse_usize(value: &str, label: &str) -> Result<usize, String> {
    value.trim().parse().map_err(|_| format!("{label} 必须是正整数，收到: {}", value.trim()))
}

fn parse_f64(value: &str, label: &str) -> Result<f64, String> {
    let parsed: f64 = value.trim().parse().map_err(|_| format!("{label} 必须是数字，收到: {}", value.trim()))?;
    if !parsed.is_finite() {
        return Err(format!("{label} 必须是有限数字，收到: {}", value.trim()));
    }
    Ok(parsed)
}

fn req_path(value: &str, label: &str) -> Result<PathBuf, String> {
    let value = value.trim();
    if value.is_empty() {
        return Err(format!("{label} 不能为空"));
    }
    let path = PathBuf::from(value);
    if !path.is_file() {
        return Err(format!("{label} 不存在: {}", path.display()));
    }
    Ok(path)
}

fn opt_path(value: &str) -> Option<PathBuf> {
    let value = value.trim();
    if value.is_empty() { None } else { Some(PathBuf::from(value)) }
}

fn opt_file_path(value: &str, label: &str) -> Result<Option<PathBuf>, String> {
    let value = value.trim();
    if value.is_empty() {
        return Ok(None);
    }
    let path = PathBuf::from(value);
    if !path.is_file() {
        return Err(format!("{label} 不存在: {}", path.display()));
    }
    Ok(Some(path))
}

fn out_path(value: &str, label: &str) -> Result<PathBuf, String> {
    let value = value.trim();
    if value.is_empty() { Err(format!("{label} 不能为空")) } else { Ok(PathBuf::from(value)) }
}

fn opt_string(value: &str) -> Option<String> {
    let value = value.trim();
    if value.is_empty() { None } else { Some(value.to_string()) }
}

fn split_list(value: &str) -> Vec<String> {
    value
        .split(|ch: char| ch == ',' || ch == ';' || ch == '|')
        .map(str::trim)
        .filter(|item| !item.is_empty())
        .map(str::to_string)
        .collect()
}

fn sanitize(name: &str) -> String {
    let mut out: String = name.chars().map(|ch| if ch.is_alphanumeric() || ch == '-' || ch == '_' { ch } else { '_' }).collect();
    let keep = out.char_indices().nth(40).map(|(index, _)| index).unwrap_or(out.len());
    out.truncate(keep);
    if out.is_empty() { "report".to_string() } else { out }
}

fn panic_message(payload: Box<dyn std::any::Any + Send>) -> String {
    if let Some(message) = payload.downcast_ref::<&str>() {
        (*message).to_string()
    } else if let Some(message) = payload.downcast_ref::<String>() {
        message.clone()
    } else {
        "未知 panic".to_string()
    }
}

// ===== 日志缓冲与入口 =====

struct BufferLogger {
    lines: Arc<Mutex<VecDeque<String>>>,
}

impl log::Log for BufferLogger {
    fn enabled(&self, metadata: &log::Metadata) -> bool {
        metadata.level() <= log::Level::Info
    }

    fn log(&self, record: &log::Record) {
        if !self.enabled(record.metadata()) {
            return;
        }
        if let Ok(mut lines) = self.lines.lock() {
            lines.push_back(format!("[{}] {}", record.level(), record.args()));
            while lines.len() > 500 {
                lines.pop_front();
            }
        }
    }

    fn flush(&self) {}
}

fn init_buffer_log() -> Arc<Mutex<VecDeque<String>>> {
    let lines = Arc::new(Mutex::new(VecDeque::new()));
    let logger = BufferLogger { lines: Arc::clone(&lines) };
    let _ = log::set_boxed_logger(Box::new(logger));
    log::set_max_level(log::LevelFilter::Info);
    lines
}

/// TUI 入口：初始化终端并进入事件循环。
///
/// `initial` 可以是待打开的波函数文件或目录；目录会作为左侧文件列表的根目录。
pub fn run(initial: Option<PathBuf>) {
    if let Err(err) = run_inner(initial) {
        eprintln!("TUI 启动失败: {err}");
    }
}

fn run_inner(initial: Option<PathBuf>) -> io::Result<()> {
    let log_lines = init_buffer_log();

    enable_raw_mode()?;
    let mut stdout = io::stdout();
    execute!(stdout, EnterAlternateScreen, EnableBracketedPaste)?;
    let backend = CrosstermBackend::new(stdout);
    let mut terminal = Terminal::new(backend)?;
    terminal.clear()?;

    // panic 时恢复终端，避免留下 raw mode/备用屏幕；
    // worker 线程的 panic 已被 catch_unwind 捕获，不再破坏界面。
    std::panic::set_hook(Box::new(|info| {
        if QUIET_PANIC.with(|flag| flag.get()) {
            let _ = std::fs::write("rswfn-panic.log", format!("{info}\n"));
            return;
        }
        let _ = disable_raw_mode();
        let mut out = io::stdout();
        let _ = execute!(out, LeaveAlternateScreen, DisableBracketedPaste);
        let _ = writeln!(out, "rswfn TUI 发生内部错误: {info}");
        let _ = std::fs::write("rswfn-panic.log", format!("{info}\n"));
    }));

    let mut app = App::new(log_lines, initial);
    let result = app.event_loop(&mut terminal);

    let _ = disable_raw_mode();
    let _ = execute!(terminal.backend_mut(), LeaveAlternateScreen, DisableBracketedPaste);
    let _ = terminal.show_cursor();
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn text_input_edits() {
        let mut input = TextInput::new("abc");
        input.left();
        input.insert('X');
        assert_eq!(input.value, "abXc");
        input.backspace();
        assert_eq!(input.value, "abc");
        input.home();
        input.delete();
        assert_eq!(input.value, "bc");
    }

    #[test]
    fn splits_lists() {
        assert_eq!(split_list("a, b; c"), vec!["a", "b", "c"]);
        assert!(split_list("  ").is_empty());
    }

    #[test]
    fn sanitizes_names() {
        assert_eq!(sanitize("a b/c"), "a_b_c");
        assert_eq!(sanitize(""), "report");
    }

    #[test]
    fn tasks_are_unique() {
        for (i, task) in TASKS.iter().enumerate() {
            assert!(!task.name.is_empty());
            for other in TASKS.iter().skip(i + 1) {
                assert_ne!(task.id, other.id);
            }
        }
    }

    fn render(app: &mut App) -> String {
        let backend = ratatui::backend::TestBackend::new(100, 32);
        let mut terminal = Terminal::new(backend).expect("terminal");
        terminal.draw(|frame| app.draw(frame)).expect("draw");
        terminal.backend().buffer().content().iter().map(|cell| cell.symbol()).collect()
    }

    /// 宽字符会占用两列，第二列的符号是空格；断言前先把空白去掉
    fn compact(text: &str) -> String {
        text.chars().filter(|ch| !ch.is_whitespace()).collect()
    }

    #[test]
    fn draws_welcome() {
        let mut app = App::new(Arc::new(Mutex::new(VecDeque::new())), None);
        let text = compact(&render(&mut app));
        assert!(text.contains("rswfn"));
        assert!(text.contains("atomprop"));
        assert!(text.contains("使用说明"));
    }

    #[test]
    fn opens_and_draws_viewer() {
        let b = crate::datas::consts::BOHR;
        let water = viewer::Molecule {
            syms: vec!["O".to_string(), "H".to_string(), "H".to_string()],
            coords: vec![[0.0, 0.0, 0.0], [0.96 * b, 0.0, 0.0], [-0.24 * b, 0.93 * b, 0.0]],
            rads: vec![0.63 * b, 0.32 * b, 0.32 * b],
            colors: vec![viewer::Rgb(229, 0, 0), viewer::Rgb(204, 204, 204), viewer::Rgb(204, 204, 204)],
            bonds: vec![[0, 1], [0, 2]],
            axes: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        };
        let mut app = App::new(Arc::new(Mutex::new(VecDeque::new())), None);
        app.header = Some(Header {
            formula: "H2O".to_string(),
            charge: "0".to_string(),
            multip: "1".to_string(),
            natoms: "3".to_string(),
            mol: water,
        });
        app.open_viewer();
        assert!(matches!(app.overlay, Some(Overlay::Viewer(_))));
        let text = compact(&render(&mut app));
        assert!(text.contains("分子结构预览"), "预览标题应出现在界面上");
        assert!(text.contains("球棍"), "标题应显示当前绘制模型");
    }

    /// 把 ratatui 的终端缓冲转成 ANSI，供文档截图使用（平时不写文件）
    fn buffer_to_ansi(buffer: &ratatui::buffer::Buffer) -> String {
        fn rgb(color: ratatui::style::Color) -> (u8, u8, u8) {
            use ratatui::style::Color;
            match color {
                Color::Rgb(r, g, b) => (r, g, b),
                Color::Black => (0, 0, 0),
                Color::Red => (205, 49, 49),
                Color::Green => (13, 188, 121),
                Color::Yellow => (229, 229, 16),
                Color::Blue => (36, 114, 200),
                Color::Magenta => (188, 63, 188),
                Color::Cyan => (17, 168, 205),
                Color::Gray => (204, 204, 204),
                Color::DarkGray => (118, 118, 118),
                Color::White => (255, 255, 255),
                _ => (205, 205, 212),
            }
        }
        let area = *buffer.area();
        let mut out = String::new();
        for y in 0..area.height {
            for x in 0..area.width {
                let cell = &buffer[(x, y)];
                let (fr, fg, fb) = rgb(cell.fg);
                let (br, bg, bb) = match cell.bg {
                    ratatui::style::Color::Reset => (18, 18, 22),
                    other => rgb(other),
                };
                out.push_str(&format!("\x1b[38;2;{fr};{fg};{fb};48;2;{br};{bg};{bb}m"));
                out.push_str(cell.symbol());
            }
            out.push_str("\x1b[0m\n");
        }
        out
    }

    /// 生成文档用的预览窗口截图：`RSWFN_SHOT=<文件>` 时把真实绘制结果写成 ANSI。
    /// 未设置环境变量时直接返回，不参与常规测试。
    #[test]
    fn dumps_viewer_screenshot() {
        let Ok(path) = std::env::var("RSWFN_SHOT") else {
            return;
        };
        let file = std::env::var("RSWFN_SHOT_FILE").unwrap_or_else(|_| "moles/C6H6.fch".to_string());
        let Some(header) = Header::from_path(Path::new(&file)) else {
            panic!("无法读取截图用分子: {file}");
        };
        let mut app = App::new(Arc::new(Mutex::new(VecDeque::new())), Some(PathBuf::from(&file)));
        app.header = Some(header);
        app.open_viewer();
        let backend = ratatui::backend::TestBackend::new(112, 38);
        let mut terminal = Terminal::new(backend).expect("terminal");
        terminal.draw(|frame| app.draw(frame)).expect("draw");
        std::fs::write(&path, buffer_to_ansi(terminal.backend().buffer())).expect("写截图失败");
    }

    #[test]
    fn draws_every_form() {
        for task in TASKS {
            let mut app = App::new(Arc::new(Mutex::new(VecDeque::new())), None);
            app.form = Some(Form::new(task.id, None));
            app.mode = Mode::Form;
            let text = compact(&render(&mut app));
            let name = compact(task.name);
            assert!(text.contains(&name), "表单缺少标题: {:?}", task.id);
        }
    }

    #[test]
    fn grid_form_defaults_and_validation() {
        let form = Form::new(TaskId::GridDensity, None);
        assert_eq!(form.text(0), "density");
        assert_eq!(form.text(1), "0.2");
        assert!(build_task(&form, None).is_err());
    }

    #[test]
    fn qm_form_prefills_current_file() {
        let form = Form::new(TaskId::QmHf, Some(Path::new("moles/C6H6.fch")));
        assert_eq!(form.text(0), "moles/C6H6.fch");
        assert_eq!(form.text(1), "");
    }

    #[test]
    fn missing_current_file_is_rejected() {
        let form = Form::new(TaskId::GridDensity, None);
        assert!(build_task(&form, Some(Path::new("no_such_file.fch"))).is_err());
    }

    #[test]
    fn runs_mulliken_job() {
        let mut app = App::new(Arc::new(Mutex::new(VecDeque::new())), None);
        let path = PathBuf::from("moles/C6H6.fch");
        app.current_file = Some(path.clone());
        app.form = Some(Form::new(TaskId::Mulliken, Some(&path)));
        app.start_job();
        assert!(matches!(app.mode, Mode::Running), "任务应立即进入运行状态");

        let deadline = Instant::now() + Duration::from_secs(30);
        while app.result.is_none() && Instant::now() < deadline {
            std::thread::sleep(Duration::from_millis(20));
            app.poll_job();
        }

        let report = app.result.as_ref().expect("任务未在超时前完成");
        assert!(report.title.contains("Mulliken"));
        assert!(matches!(app.mode, Mode::Result));
    }

    #[test]
    fn rejects_missing_file_before_running() {
        let mut app = App::new(Arc::new(Mutex::new(VecDeque::new())), None);
        let path = PathBuf::from("no_such_file.fch");
        app.current_file = Some(path.clone());
        app.form = Some(Form::new(TaskId::Mulliken, Some(&path)));
        app.start_job();
        assert!(matches!(app.overlay, Some(Overlay::Error(_))));
        assert!(app.job.is_none());
    }

    #[test]
    fn tree_covers_every_task() {
        fn collect(nodes: &'static [Node], out: &mut Vec<TaskId>) {
            for node in nodes {
                if let Some(task) = node_task(node) {
                    out.push(task);
                }
                collect(node.children, out);
            }
        }
        let mut tasks = Vec::new();
        collect(TREE, &mut tasks);
        assert_eq!(tasks.len(), TASKS.len());
        for task in TASKS {
            assert!(tasks.contains(&task.id), "菜单缺少 {:?}", task.id);
        }
    }

    #[test]
    fn tree_collapses_and_expands() {
        let mut app = App::new(Arc::new(Mutex::new(VecDeque::new())), None);
        app.focus = Focus::Tree;
        app.menu_index = 0;
        let collapsed = app.tree_rows().len();
        app.activate_selection();
        let expanded = app.tree_rows().len();
        assert!(expanded > collapsed, "展开后行数应增加");
        app.activate_selection();
        assert_eq!(app.tree_rows().len(), collapsed);
    }

    #[test]
    fn file_panel_lists_moles() {
        if !Path::new("moles").is_dir() {
            return;
        }
        let panel = FilePanel::new(PathBuf::from("moles"));
        assert!(panel.entries.iter().any(|entry| entry.name == "C6H6.fch"));
    }

    #[test]
    fn initial_directory_sets_file_panel() {
        if !Path::new("moles").is_dir() {
            return;
        }
        let app = App::new(Arc::new(Mutex::new(VecDeque::new())), Some(PathBuf::from("moles")));
        assert!(app.files.dir.ends_with("moles"));
        assert!(matches!(app.focus, Focus::Files));
    }

    #[test]
    fn initial_file_is_selected() {
        if !Path::new("moles/C6H6.fch").is_file() {
            return;
        }
        let app = App::new(Arc::new(Mutex::new(VecDeque::new())), Some(PathBuf::from("moles/C6H6.fch")));
        assert!(app.current_file.is_some());
        assert!(matches!(app.focus, Focus::Tree));
        assert!(app.files.entries.iter().any(|entry| entry.name == "C6H6.fch"));
    }

    #[test]
    fn draws_file_panel_and_tree() {
        if !Path::new("moles").is_dir() {
            return;
        }
        let mut app = App::new(Arc::new(Mutex::new(VecDeque::new())), Some(PathBuf::from("moles")));
        let text = compact(&render(&mut app));
        assert!(text.contains("C6H6.fch"), "文件列表应显示 moles 下的文件");
        assert!(text.contains("分析菜单"));
        assert!(text.contains("atomprop"), "树应包含 CLI 的 atomprop 分组");
        assert!(text.contains("结构预览"), "树应包含 CLI 的 view 入口");
    }

    #[test]
    fn viewer_tree_entry_opens_viewer() {
        let mut app = App::new(Arc::new(Mutex::new(VecDeque::new())), None);
        app.focus = Focus::Tree;
        let rows = app.tree_rows();
        let index = rows.iter().position(|row| matches!(row.node.kind, NodeKind::Viewer)).expect("树中应有 view 入口");
        app.menu_index = index;
        app.activate_selection();
        assert!(app.status.contains("请先打开"), "无分子时应提示先打开文件，实际: {}", app.status);
        assert!(app.overlay.is_none(), "没有分子时不应打开空白预览窗口");
    }

    #[test]
    fn browser_popup_has_border() {
        let mut app = App::new(Arc::new(Mutex::new(VecDeque::new())), None);
        app.open_browser(BrowserTarget::CurrentFile);
        let text = render(&mut app);
        assert!(text.contains("┌") && text.contains("┘"), "弹出层应有边框，实际渲染:\n{text}");
        assert!(compact(&text).contains("打开文件"), "弹出层标题缺失，实际渲染:\n{text}");
    }

    #[test]
    fn f3_opens_results_popup() {
        let mut app = App::new(Arc::new(Mutex::new(VecDeque::new())), None);
        app.history.push(Snapshot {
            label: "Mulliken 布居".to_string(),
            cli: "atomprop charge mulliken a.fch".to_string(),
            time: "12:00:00".to_string(),
            text: "=== 结果 ===\nC1  0.123".to_string(),
        });
        app.on_key(KeyEvent::new(KeyCode::F(3), KeyModifiers::NONE));
        assert!(matches!(app.overlay, Some(Overlay::Results { .. })), "F3 应打开结果弹窗");
        let text = render(&mut app);
        let flat: String = text.chars().filter(|ch| !ch.is_whitespace()).collect();
        assert!(flat.contains("F3计算结果"), "弹窗标题缺失，实际:\n{text}");
        assert!(flat.contains("等价命令"), "弹窗里应显示等价 CLI，实际:\n{text}");
        app.on_key(KeyEvent::new(KeyCode::Esc, KeyModifiers::NONE));
        assert!(app.overlay.is_none(), "Esc 应关闭弹窗");
    }

    #[test]
    fn esc_in_menu_does_not_quit() {
        let mut app = App::new(Arc::new(Mutex::new(VecDeque::new())), None);
        app.on_key(KeyEvent::new(KeyCode::Esc, KeyModifiers::NONE));
        assert!(!app.quit, "最外层 Esc 不应退出");
        assert!(app.status.contains('q'), "应提示退出方式，实际: {}", app.status);
        // 功能层按 Esc 先回文件层（层级回退），不退出
        app.focus = Focus::Tree;
        app.on_key(KeyEvent::new(KeyCode::Esc, KeyModifiers::NONE));
        assert!(!app.quit, "功能层 Esc 不应退出");
        assert_eq!(app.focus, Focus::Files, "功能层 Esc 应回到文件层");
        app.on_key(KeyEvent::new(KeyCode::Char('q'), KeyModifiers::NONE));
        assert!(app.quit, "q 仍应退出");
    }
}
