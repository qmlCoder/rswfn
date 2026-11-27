//! 分子级别的性质：芳香性、能量（π/σ 能量分解）、多极矩、布居分析，以及
//! 自旋分离、π 电子投影等分子变体构造。
/// 芳香性计算
pub mod aromacity;
/// 能量计算（π/σ 能量分解等）
pub mod energy;
/// 多极矩计算
pub mod multipole;
/// 布居分析
pub mod population;
/// 分子变体构造（自旋分离、π 电子投影）
pub mod variants;
