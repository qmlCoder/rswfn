use std::process::exit;
// use std::char;

use regex::Regex;

use crate::base::{Atoms, Basis, Coefs, Stm};
use crate::datas::consts::BOHR;
use crate::datas::elems::sym2atm;
use crate::reader::charge_number;
use crate::reader::{Reader, read_lines};
use crate::rsbse;
use crate::rsxcs::XcFunc;

/// 高斯输入文件(.gjf)读取器
pub struct GjfReader {
    path: String,
    lines: Vec<String>,
}

impl GjfReader {
    /// 创建GjfReader实例
    pub fn new(path: &str) -> Self {
        let lines = read_lines(path);
        Self {
            path: path.to_string(),
            lines,
        }
    }

    /// 读取几何结构
    pub fn read_geome(&self) -> (Vec<usize>, Vec<[f64; 3]>) {
        let mut atms = vec![];
        let mut xyzs = vec![];
        let s1 = Regex::new(r" *([A-Za-z\d]+) +(-?\d+.\d+) +(-?\d+.\d+) +(-?\d+.\d+)").expect("正则表达式生成失败");
        for line in self.lines.iter() {
            if s1.is_match(line) {
                let caps = s1.captures(line).expect("捕获失败");
                let sym = caps[1].to_string();
                let atm = sym2atm(&sym);
                let x = caps[2].parse::<f64>().expect("解析浮点数失败") * BOHR;
                let y = caps[3].parse::<f64>().expect("解析浮点数失败") * BOHR;
                let z = caps[4].parse::<f64>().expect("解析浮点数失败") * BOHR;
                atms.push(atm);
                xyzs.push([x, y, z]);
            }
        }
        (atms, xyzs)
    }

    /// 读取电荷和自旋多重度
    pub fn read_charge(&self) -> (isize, usize) {
        let s0 = Regex::new(r"^ *(-?\d+) +(-?\d+)").expect("正则表达式生成失败");
        let mut C = 0; // 电荷
        let mut S = 0; // 自旋
        let mut matched = false;
        for line in self.lines.iter() {
            if s0.is_match(line) {
                let caps = s0.captures(line).expect("捕获失败");
                C = caps[1].parse::<isize>().expect("转换失败");
                S = caps[2].parse::<usize>().expect("转换失败");
                matched = true;
                break;
            }
        }
        if !matched {
            log::error!("无法读取电荷和自旋");
            exit(0)
        }
        let (atms, _) = self.read_geome();
        (C, S)
    }

    /// 读取作业关键词
    pub fn read_job(&self) -> String {
        let mut job = String::new();
        for line in self.lines.iter() {
            if line.starts_with("#") {
                job.push_str(&line.replace("#", "").trim().to_string());
            }
        }
        job
    }

    /// 解析作业关键词，得到任务类型、泛函与基组
    ///
    /// 见 [`parse_job`]；泛函/基组的识别由 `rsxcs` / `rsbse` 支持名单完成。
    pub fn read_gau_job(&self) -> Option<GauJob> {
        parse_job(&self.read_job())
    }
}

impl Reader for GjfReader {
    fn get_atoms(&self) -> Atoms {
        let (atms, xyzs) = self.read_geome();
        let crgs: Vec<f64> = atms.iter().map(|atm| *atm as f64).collect();
        let mut geome = Atoms::new(&atms, &xyzs, &crgs);
        let N = atms.iter().map(|atm| *atm as usize).sum::<usize>();
        let (C, S) = self.read_charge();

        geome.neles = charge_number(C, S, N);
        geome
    }
    fn get_basis(&self) -> Option<Basis> {
        None
    }

    fn get_coefs(&self) -> Option<Coefs> {
        None
    }
}

/// 从高斯作业关键词(route section)解析出的作业信息
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GauJob {
    /// 任务类型: `energy` / `opt` / `opt+freq`
    pub job_type: String,
    /// 交换相关泛函名（如 `b3lyp`）；未识别时为 `None`（通常表示 HF 或本库不支持）
    pub functional: Option<String>,
    /// 基组名（如 `6-31g(d)`）；未识别时为 `None`
    pub basis_set: Option<String>,
}

/// 解析高斯作业关键词行(形如 `#p b3lyp/6-31g(d) opt freq`)
///
/// - 任务类型由 `opt` / `freq` 关键词判断；
/// - 泛函与基组分别借助 [`crate::rsxcs`]（泛函）与 [`crate::rsbse`]（基组）
///   支持的名单来识别：优先按 `泛函/基组` 形式拆分，否则在词元中逐个识别。
///   泛函用 [`XcFunc::from_name`] 判断（自动忽略受限/非受限前缀 `R`/`U`），
///   基组则与 [`crate::rsbse::supported_names`] 做归一化比较。
pub fn parse_job(line: &str) -> Option<GauJob> {
    let line = line.trim();
    if line.is_empty() {
        return None;
    }
    let lower = line.to_lowercase();

    // 任务类型
    let job_type = if lower.contains("opt") {
        if lower.contains("freq") { "opt+freq" } else { "opt" }
    } else {
        "energy"
    };

    let supported = rsbse::supported_names();
    let (functional, basis_set) = parse_func_basis(&lower, &supported);

    Some(GauJob {
        job_type: job_type.to_string(),
        functional,
        basis_set,
    })
}

/// 从作业行中提取泛函与基组名
fn parse_func_basis(lower: &str, supported: &[String]) -> (Option<String>, Option<String>) {
    // 1) 优先按 "泛函/基组" 形式解析, 如 b3lyp/6-31g(d)、RB3LYP/6-31G(d)
    for token in lower.split_whitespace() {
        // 排除 iop(3/33=1) 之类的带 '=' 关键词
        if token.contains('=') {
            continue;
        }
        let Some((left, right)) = token.split_once('/') else {
            continue;
        };
        let left = trim_name(left);
        let right = trim_name(right);
        if !is_name_like(left) || !is_basis_like(right) {
            continue;
        }
        let functional = recognize_functional(left).unwrap_or_else(|| left.to_string());
        let basis_set = canonical_basis(right, supported);
        return (Some(functional), Some(basis_set));
    }

    // 2) 否则在词元中逐个识别（如 "m06l def2svp"）
    let mut functional = None;
    let mut basis_set = None;
    for token in lower.split_whitespace() {
        if token.contains('=') {
            continue;
        }
        let token = trim_name(token);
        if token.is_empty() {
            continue;
        }
        if functional.is_none() {
            if let Some(name) = recognize_functional(token) {
                functional = Some(name);
                continue;
            }
        }
        if basis_set.is_none() {
            if let Some(name) = match_basis(token, supported) {
                basis_set = Some(name);
            }
        }
    }
    (functional, basis_set)
}

/// 去掉名字两端空白与逗号
fn trim_name(s: &str) -> &str {
    s.trim().trim_matches(',')
}

/// 是否包含字母
fn has_alpha(s: &str) -> bool {
    s.chars().any(|c| c.is_ascii_alphabetic())
}

/// 是否像一个泛函名（字母数字与连字符，且含字母）
fn is_name_like(s: &str) -> bool {
    has_alpha(s) && s.chars().all(|c| c.is_ascii_alphanumeric() || c == '-')
}

/// 是否像一个基组名（允许括号、逗号、加减号、星号等）
fn is_basis_like(s: &str) -> bool {
    has_alpha(s)
        && s.chars()
            .all(|c| c.is_ascii_alphanumeric() || matches!(c, '-' | '(' | ')' | ',' | '+' | '*'))
}

/// 借助 [`XcFunc::from_name`] 识别泛函；自动去掉受限/非受限前缀 `R`/`U`
fn recognize_functional(token: &str) -> Option<String> {
    let token = token.to_lowercase();
    if XcFunc::from_name(&token).is_some() {
        return Some(token);
    }
    if let Some(rest) = token.strip_prefix('r').or_else(|| token.strip_prefix('u')) {
        if XcFunc::from_name(rest).is_some() {
            return Some(rest.to_string());
        }
    }
    None
}

/// 基组名归一化：小写并去掉非字母数字字符（`def2-SVP` 与 `def2svp` 视为同名）
fn basis_key(s: &str) -> String {
    s.chars()
        .filter(|c| c.is_ascii_alphanumeric())
        .map(|c| c.to_ascii_lowercase())
        .collect()
}

/// 在受支持基组中查找匹配项，返回其规范名
fn match_basis(token: &str, supported: &[String]) -> Option<String> {
    let key = basis_key(token);
    if key.is_empty() {
        return None;
    }
    supported.iter().find(|name| basis_key(name) == key).cloned()
}

/// 优先返回规范名，找不到则原样返回（小写）
fn canonical_basis(token: &str, supported: &[String]) -> String {
    match_basis(token, supported).unwrap_or_else(|| token.to_lowercase())
}

#[cfg(test)]
mod tests {
    use super::parse_job;

    #[test]
    fn parse_route_line() {
        let job = parse_job("# opt freq b3lyp/6-31g(d) pop=full gfinput iop(3/33=1)").unwrap();
        assert_eq!(job.job_type, "opt+freq");
        assert_eq!(job.functional.as_deref(), Some("b3lyp"));
        assert_eq!(super::basis_key(&job.basis_set.unwrap()), "631gd");
    }

    #[test]
    fn parse_route_with_restricted_prefix() {
        let job = parse_job("#N Geom=AllCheck Guess=TCheck RB3LYP/6-31G(d)").unwrap();
        assert_eq!(job.job_type, "energy");
        assert_eq!(job.functional.as_deref(), Some("b3lyp"));
        assert_eq!(super::basis_key(&job.basis_set.unwrap()), "631gd");
    }

    #[test]
    fn parse_basis_only_line() {
        let job = parse_job("# m06l def2svp pop=full gfinput nosymm").unwrap();
        assert_eq!(job.job_type, "energy");
        // def2svp 归一化后匹配内嵌基组 def2-svp
        assert_eq!(super::basis_key(&job.basis_set.unwrap()), "def2svp");
    }

    #[test]
    fn empty_line_returns_none() {
        assert!(parse_job("   ").is_none());
    }
}
