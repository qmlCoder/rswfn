//! 波函数分析库
pub mod base;
pub mod datas;

pub mod editor;
pub mod loger;
pub mod maths;
pub mod props;
pub mod qmcal;
pub mod reader;
pub mod rsbse;
pub mod rsint;
pub mod rsxcs; // 泛函
pub mod space;
pub mod tools;
pub mod utils;
pub mod writer;

// 保持旧路径可用（`rswfn::march` 等），子模块已归入 `maths`
pub use maths::{march, matrix, optim};
