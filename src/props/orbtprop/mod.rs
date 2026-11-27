//! 轨道级别的性质：MOCV 投影（`mocv`）、分子轨道系数矩阵操作（`obtmat`）
//! 以及分子轨道电子数等轨道取值（`obtval`）。
/// 局部坐标系下的 MOCV 投影与坐标变换
pub mod mocv;
/// 分子轨道系数矩阵操作（POCV、MOCV、轨道定域化等）
pub mod obtmat;
/// 分子轨道属性与轨道电子数计算
pub mod obtval;
