//! 键性质计算：键级、键能分解与键方向。

/// 键轴周围方向向量的生成。
pub mod direction;
/// 键能量的可加分解。
pub mod energy;
mod energy_v1;
/// 键级（Mayer、Wiberg 等）计算。
pub mod order;
