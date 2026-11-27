//! 各类性质（电荷、键级、能量、格点、轨道等）计算模块。

/// 原子性质：原子电荷/电子布居、分子活性与局部方向。
pub mod atomprop;
/// 键性质：键级、键能分解与键方向。
pub mod bondprop;
/// 片段性质：片段（原子集合）间的相互作用能。
pub mod fragprop;
/// 格点性质：密度格点、DFT 格点与拓扑分析。
pub mod gridprop;
/// 分子整体性质：能量、多极矩、芳香性与电子布居。
pub mod moleprop;
/// 轨道性质：轨道变换、轨道占据与成键分析。
pub mod orbtprop;
