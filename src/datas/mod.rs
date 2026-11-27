/// 基函数变换矩阵（笛卡尔/球谐互转）与轨道符号、角动量分量映射
pub mod btrans;
/// 物理常数
pub mod consts;
/// LDA 交换相关泛函（Dirac）能量与导数计算
pub mod dftfc;
/// 双原子分子轨道（DHMO）参数
pub mod dhmo2;
/// 元素数据（符号、原子序数、半径、周期表位置、颜色、电子排布与氧化数）
pub mod elems;
/// 高斯积分所需的 Hermite 多项式节点与权重
pub mod integ;
/// Lebedev 球面积分网格生成
pub mod lebedev;
/// 原子的径向电子密度分布数据
pub mod rdens;
/// Legendre/Hermite 求积节点权重与 Rys 求根（插值法）
pub mod roots;
/// Rys 多项式求根数据与算法（数据来自 libcint 源代码）
pub mod rys_root;
/// 输出文件（GJF/XYZ 等）文本模板
pub mod temps;

// use std::collections::HashMap;

// let mut cache: HashMap<String, CacheTypes> = HashMap::new();

// enum CacheTypes {
//     Int(usize),
//     Float(f64),
// }

// pub static RAND_VECT: LazyLock<Vector3<f64>> = LazyLock::new(|| {
//     //随机向量
//     let x = rand::random();
//     let y = rand::random();
//     let z = rand::random();
//     Vector3::new(x, y, z).normalize()
// });

/// 双阶乘 (2i-1)!! 的前几项
pub const FACS: [f64; 6] = [1., 1., 3., 15., 105., 945.];
