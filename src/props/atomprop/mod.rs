//! 原子性质计算器

/// 分子活性（福井函数、自由价）计算。
pub mod activity;
/// 原子电荷与电子布居分析。
pub mod charge;
/// 原子局部方向与法向量计算。
pub mod direction;

/// RESP（限制性静电势）电荷，实现位于 `_resp.rs`
#[path = "_resp.rs"]
pub mod resp;
