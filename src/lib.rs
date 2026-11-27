//! 波函数分析库
//!
//! 提供分子结构、基组、分子轨道系数等核心数据结构，以及量子化学计算和
//! 各类波函数分析（电荷、键级、轨道、格点等）的实现。
//!
//! 各子模块大致分工如下：
//! - [`base`]：核心数据结构（分子、原子、基组、轨道系数、三维向量）
//! - [`datas`]：元素、物理常数、转换表等基础数据
//! - [`editor`]：分子结构编辑（键的旋转、伸缩等几何操作）
//! - [`loger`]：日志器
//! - [`maths`]：矩阵运算、优化、等值面提取等数学工具
//! - [`props`]：电荷、键级、能量、格点、轨道等性质计算
//! - [`qmcal`]：SCF/DFT 等量子化学计算
//! - [`reader`]：各种格式文件的读取器
//! - [`rsbse`]：基组库（内置基组与 Basis Set Exchange 数据库）
//! - [`rsint`]：单/双电子积分计算
//! - [`rsxcs`]：交换相关泛函
//! - [`space`]：空间格点生成与原子权重分配
//! - [`tools`]：辅助工具
//! - [`utils`]：通用工具函数
//! - [`writer`]：各种格式文件的写入器

/// 核心数据结构：分子、原子、基组、轨道系数与三维向量
pub mod base;
/// 元素、物理常数、转换表等基础数据
pub mod datas;

/// 分子结构编辑器，提供键的旋转和伸缩等几何操作
pub mod editor;
/// 日志器（仅写 stderr 的轻量实现）
pub mod loger;
/// 数学工具：矩阵运算、优化、等值面提取等
pub mod maths;
/// 各类性质（电荷、键级、能量、格点、轨道等）计算
pub mod props;
/// 量子化学计算（SCF/DFT 等）
pub mod qmcal;
pub mod reader;
/// 基组库（内置基组与 Basis Set Exchange 数据库）
pub mod rsbse;
/// 单/双电子积分计算
pub mod rsint;
pub mod rsxcs; // 泛函
/// 空间格点生成与原子权重分配
pub mod space;
/// 辅助工具
pub mod tools;
/// 通用工具函数
pub mod utils;
/// 各种格式文件的写入器
pub mod writer;

// 保持旧路径可用（`rswfn::march` 等），子模块已归入 `maths`
pub use maths::{march, matrix, optim};
