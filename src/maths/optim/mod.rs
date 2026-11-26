//! 无约束优化算法
//!
//! 目前提供最速下降法（[`sd`]）、非线性共轭梯度法（[`cg`]）、牛顿法（[`nt`]）、
//! 阻尼牛顿法（[`dn`]）和拟牛顿法（[`qn`]）。
//!
//! 目标函数以闭包形式传入，闭包接受自变量 `paras: &Vec<f64>`，返回
//! [`FuncEval`]。函数值必须提供，梯度和 Hessian 为可选项：
//! - 未提供梯度时，使用中心差分近似梯度；
//! - Hessian 未提供时，若已有解析梯度则对梯度做中心差分，否则对函数值做二阶中心差分；
//! - 最速下降/共轭梯度/拟牛顿法只把 Hessian 用于估计初始步长，牛顿法直接使用 Hessian；
//!   拟牛顿法（BFGS）通过梯度差商逐步构造逆 Hessian 的近似。
//!
//! 收敛判据统一为梯度的无穷范数：`max_i |∂F/∂x_i| <= conv_thresh`。
//!
//! # 示例
//!
//! ```ignore
//! use rswfn::maths::optim::{FuncEval, steepest_descent};
//!
//! // F(x, y) = x^2 + 2y^2，梯度为 (2x, 4y)
//! let func = |paras: &Vec<f64>| {
//!     let x = paras[0];
//!     let y = paras[1];
//!     FuncEval::from_grad(x * x + 2.0 * y * y, vec![2.0 * x, 4.0 * y])
//! };
//!
//! let paras = vec![1.0, 1.0];
//! let result = steepest_descent(func, &paras, 1e-6).unwrap();
//! ```
//!
//! 牛顿法/阻尼牛顿法的调用方式相同，只需将 `steepest_descent` 替换为
//! `newton`、`damped_newton` 或 `quasi_newton`，并提供（或让算法数值近似）
//! Hessian；拟牛顿法本身只需要梯度。

use ndarray::Array2;

pub mod cg;
pub mod dn;
pub mod nt;
pub mod qn;
pub mod sd;
mod utils;

pub use cg::minimize as conjugate_gradient;
pub use dn::minimize as damped_newton;
pub use nt::minimize as newton;
pub use qn::minimize as quasi_newton;
pub use sd::minimize as steepest_descent;

/// 目标函数在某个自变量处的求值结果。
///
/// `grad` 和 `hessian` 均可为 `None`，算法会按需计算数值梯度和数值 Hessian。
#[derive(Clone, Debug)]
pub struct FuncEval {
    /// 函数值
    pub value: f64,
    /// 梯度，长度必须与自变量维数相同
    pub grad: Option<Vec<f64>>,
    /// Hessian 矩阵，形状必须为 `(n, n)`
    pub hessian: Option<Array2<f64>>,
}

impl FuncEval {
    /// 创建完整的求值结果。
    pub fn new(value: f64, grad: Option<Vec<f64>>, hessian: Option<Array2<f64>>) -> Self {
        Self { value, grad, hessian }
    }

    /// 仅提供函数值。
    pub fn from_value(value: f64) -> Self {
        Self { value, grad: None, hessian: None }
    }

    /// 提供函数值和梯度。
    pub fn from_grad(value: f64, grad: Vec<f64>) -> Self {
        Self { value, grad: Some(grad), hessian: None }
    }

    /// 提供函数值和 Hessian。
    pub fn from_hessian(value: f64, hessian: Array2<f64>) -> Self {
        Self {
            value,
            grad: None,
            hessian: Some(hessian),
        }
    }
}

impl From<(f64, Option<Vec<f64>>, Option<Array2<f64>>)> for FuncEval {
    fn from((value, grad, hessian): (f64, Option<Vec<f64>>, Option<Array2<f64>>)) -> Self {
        Self { value, grad, hessian }
    }
}

impl From<f64> for FuncEval {
    fn from(value: f64) -> Self {
        Self::from_value(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// F(x, y) = x^2 + 2y^2，梯度 (2x, 4y)
    fn quadratic_grad(paras: &Vec<f64>) -> FuncEval {
        let x = paras[0];
        let y = paras[1];
        FuncEval::from_grad(x * x + 2.0 * y * y, vec![2.0 * x, 4.0 * y])
    }

    /// F(x, y) = x^2 + 2y^2，仅提供函数值
    fn quadratic_value(paras: &Vec<f64>) -> FuncEval {
        let x = paras[0];
        let y = paras[1];
        FuncEval::from_value(x * x + 2.0 * y * y)
    }

    /// F(x, y) = x^2 + 2y^2，仅提供函数值和 Hessian
    fn quadratic_hess(paras: &Vec<f64>) -> FuncEval {
        let x = paras[0];
        let y = paras[1];
        let hessian = Array2::from_shape_vec((2, 2), vec![2.0, 0.0, 0.0, 4.0]).unwrap();
        FuncEval::from_hessian(x * x + 2.0 * y * y, hessian)
    }

    /// Rosenbrock 函数及其解析梯度
    fn rosenbrock(paras: &Vec<f64>) -> FuncEval {
        let x = paras[0];
        let y = paras[1];
        let value = (1.0 - x).powi(2) + 100.0 * (y - x * x).powi(2);
        let gx = -2.0 * (1.0 - x) - 400.0 * x * (y - x * x);
        let gy = 200.0 * (y - x * x);
        FuncEval::from_grad(value, vec![gx, gy])
    }

    /// Rosenbrock 函数及其解析梯度和 Hessian
    fn rosenbrock_hess(paras: &Vec<f64>) -> FuncEval {
        let x = paras[0];
        let y = paras[1];
        let value = (1.0 - x).powi(2) + 100.0 * (y - x * x).powi(2);
        let gx = -2.0 * (1.0 - x) - 400.0 * x * (y - x * x);
        let gy = 200.0 * (y - x * x);
        let hxx = 2.0 - 400.0 * y + 1200.0 * x * x;
        let hxy = -400.0 * x;
        let hyy = 200.0;
        let hess = Array2::from_shape_vec((2, 2), vec![hxx, hxy, hxy, hyy]).unwrap();
        FuncEval::new(value, Some(vec![gx, gy]), Some(hess))
    }

    /// F(x, y) = x^4 + y^4 + x^2 + y^2，提供解析梯度
    fn quartic_grad(paras: &Vec<f64>) -> FuncEval {
        let x = paras[0];
        let y = paras[1];
        let value = x.powi(4) + y.powi(4) + x * x + y * y;
        let gx = 4.0 * x.powi(3) + 2.0 * x;
        let gy = 4.0 * y.powi(3) + 2.0 * y;
        FuncEval::from_grad(value, vec![gx, gy])
    }

    /// F(x, y) = x^4 + y^4 + x^2 + y^2，仅提供函数值
    fn quartic_value(paras: &Vec<f64>) -> FuncEval {
        let x = paras[0];
        let y = paras[1];
        FuncEval::from_value(x.powi(4) + y.powi(4) + x * x + y * y)
    }

    #[test]
    fn test_sd_quadratic_with_grad() {
        let start = vec![1.0, 1.0];
        let result = sd::minimize(quadratic_grad, &start, 1e-8).unwrap();
        assert!(result[0].abs() < 1e-4, "{result:?}");
        assert!(result[1].abs() < 1e-4, "{result:?}");
    }

    #[test]
    fn test_sd_quadratic_numeric_grad() {
        let start = vec![1.0, 1.0];
        let result = sd::minimize(quadratic_value, &start, 1e-6).unwrap();
        assert!(result[0].abs() < 1e-3, "{result:?}");
        assert!(result[1].abs() < 1e-3, "{result:?}");
    }

    #[test]
    fn test_cg_quadratic_with_grad() {
        let start = vec![1.0, 1.0];
        let result = cg::minimize(quadratic_grad, &start, 1e-8).unwrap();
        assert!(result[0].abs() < 1e-5, "{result:?}");
        assert!(result[1].abs() < 1e-5, "{result:?}");
    }

    #[test]
    fn test_cg_quadratic_numeric_grad() {
        let start = vec![1.0, 1.0];
        let result = cg::minimize(quadratic_value, &start, 1e-6).unwrap();
        assert!(result[0].abs() < 1e-3, "{result:?}");
        assert!(result[1].abs() < 1e-3, "{result:?}");
    }

    #[test]
    fn test_sd_quadratic_hess_only() {
        let start = vec![1.0, 1.0];
        let result = sd::minimize(quadratic_hess, &start, 1e-6).unwrap();
        assert!(result[0].abs() < 1e-3, "{result:?}");
        assert!(result[1].abs() < 1e-3, "{result:?}");
    }

    #[test]
    fn test_cg_quadratic_hess_only() {
        let start = vec![1.0, 1.0];
        let result = cg::minimize(quadratic_hess, &start, 1e-6).unwrap();
        assert!(result[0].abs() < 1e-3, "{result:?}");
        assert!(result[1].abs() < 1e-3, "{result:?}");
    }

    #[test]
    fn test_cg_rosenbrock() {
        let start = vec![-1.2, 1.0];
        let result = cg::minimize(rosenbrock, &start, 1e-7).unwrap();
        assert!((result[0] - 1.0).abs() < 1e-2, "{result:?}");
        assert!((result[1] - 1.0).abs() < 1e-2, "{result:?}");
    }

    #[test]
    fn test_nt_quadratic_hess_only() {
        let start = vec![1.0, 1.0];
        let result = nt::minimize(quadratic_hess, &start, 1e-8).unwrap();
        assert!(result[0].abs() < 1e-8, "{result:?}");
        assert!(result[1].abs() < 1e-8, "{result:?}");
    }

    #[test]
    fn test_nt_quartic_numeric_hess() {
        let start = vec![1.0, 1.0];
        let result = nt::minimize(quartic_grad, &start, 1e-8).unwrap();
        assert!(result[0].abs() < 1e-5, "{result:?}");
        assert!(result[1].abs() < 1e-5, "{result:?}");
    }

    #[test]
    fn test_dn_rosenbrock() {
        let start = vec![-1.2, 1.0];
        let result = dn::minimize(rosenbrock_hess, &start, 1e-7).unwrap();
        assert!((result[0] - 1.0).abs() < 1e-3, "{result:?}");
        assert!((result[1] - 1.0).abs() < 1e-3, "{result:?}");
    }

    #[test]
    fn test_dn_quartic_value_only() {
        let start = vec![1.0, 1.0];
        let result = dn::minimize(quartic_value, &start, 1e-6).unwrap();
        assert!(result[0].abs() < 1e-3, "{result:?}");
        assert!(result[1].abs() < 1e-3, "{result:?}");
    }

    #[test]
    fn test_nt_rejects_indefinite_hessian() {
        // F(x) = -x^2，Hessian 为负，牛顿方向不是下降方向
        let func = |paras: &Vec<f64>| {
            let x = paras[0];
            let hess = Array2::from_shape_vec((1, 1), vec![-2.0]).unwrap();
            FuncEval::new(-x * x, Some(vec![-2.0 * x]), Some(hess))
        };
        assert!(nt::minimize(func, &vec![1.0], 1e-8).is_err());
    }

    #[test]
    fn test_qn_quadratic_with_grad() {
        let start = vec![1.0, 1.0];
        let result = qn::minimize(quadratic_grad, &start, 1e-8).unwrap();
        assert!(result[0].abs() < 1e-6, "{result:?}");
        assert!(result[1].abs() < 1e-6, "{result:?}");
    }

    #[test]
    fn test_qn_quadratic_numeric_grad() {
        let start = vec![1.0, 1.0];
        let result = qn::minimize(quadratic_value, &start, 1e-6).unwrap();
        assert!(result[0].abs() < 1e-3, "{result:?}");
        assert!(result[1].abs() < 1e-3, "{result:?}");
    }

    #[test]
    fn test_qn_quadratic_hess_only() {
        let start = vec![1.0, 1.0];
        let result = qn::minimize(quadratic_hess, &start, 1e-8).unwrap();
        assert!(result[0].abs() < 1e-6, "{result:?}");
        assert!(result[1].abs() < 1e-6, "{result:?}");
    }

    #[test]
    fn test_qn_rosenbrock() {
        let start = vec![-1.2, 1.0];
        let result = qn::minimize(rosenbrock, &start, 1e-7).unwrap();
        assert!((result[0] - 1.0).abs() < 1e-3, "{result:?}");
        assert!((result[1] - 1.0).abs() < 1e-3, "{result:?}");
    }

    #[test]
    fn test_invalid_threshold() {
        let start = vec![1.0, 1.0];
        assert!(sd::minimize(quadratic_grad, &start, 0.0).is_err());
        assert!(cg::minimize(quadratic_grad, &start, -1.0).is_err());
        assert!(nt::minimize(quadratic_grad, &start, 0.0).is_err());
        assert!(dn::minimize(quadratic_grad, &start, -1.0).is_err());
        assert!(qn::minimize(quadratic_grad, &start, 0.0).is_err());
    }
}
