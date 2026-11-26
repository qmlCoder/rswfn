// 所有量子化学计算器的目的都是根据分子结构和基组信息来计算波函数信息
// 波函数=f(分子结构，基组信息)

use std::f64::consts::PI;

pub mod dhmo;
pub mod dhmo2;
pub mod diis;
pub mod guess;
pub mod hf;
pub mod ks;
pub mod lks;
pub mod se;

/// 计算Slater交换泛函的数值
///
/// 基于均匀电子气模型的交换能密度公式:
/// \(E_x[\rho] = -\frac{3}{4}\left(\frac{3}{\pi}\right)^{1/3} \rho^{4/3}\)
pub fn SlaterX(rho: f64) -> f64 {
    -3.0 / 4.0 * (3.0 / PI).powf(1.0 / 3.0) * rho.powf(4.0 / 3.0)
}
