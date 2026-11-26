//! 牛顿法（Newton Method）
//!
//! 每步采用完整牛顿步 `x_{k+1} = x_k - H_k^{-1} g_k`，不做线搜索。
//! 该方法在极小点附近具有二次收敛速度，但要求初始点足够接近极小点，
//! 且 Hessian 正定、可逆；否则可能返回错误。
//!
//! 收敛判据为梯度的无穷范数：`max_i |∂F/∂x_i| <= conv_thresh`。

use super::FuncEval;
use super::utils::{MAX_ITER, axpy, dot, gradient, hessian, newton_step, norm_inf, validate_inputs};

/// 使用牛顿法最小化目标函数。
///
/// # 参数
/// - `func`: 目标函数，输入为自变量 `&Vec<f64>`，输出为 [`FuncEval`]。
///   若未提供梯度或 Hessian，则自动用有限差分近似。
/// - `paras`: 初始自变量。
/// - `conv_thresh`: 收敛阈值，当梯度的无穷范数小于该值时认为收敛。
///
/// # 返回
/// 收敛后的自变量。若 Hessian 奇异、牛顿方向不是下降方向，或达到最大迭代
/// 次数仍未收敛，返回 `Err`。
pub fn minimize<F>(func: F, paras: &Vec<f64>, conv_thresh: f64) -> Result<Vec<f64>, String>
where
    F: Fn(&Vec<f64>) -> FuncEval,
{
    validate_inputs(paras, conv_thresh)?;

    let mut x = paras.clone();
    let mut eval = func(&x);
    if !eval.value.is_finite() {
        return Err(format!("初始点的函数值非有限: {}", eval.value));
    }

    let mut grad = gradient(&func, &x, &eval)?;
    let mut grad_norm = norm_inf(&grad);
    if grad_norm <= conv_thresh {
        return Ok(x);
    }

    for iter in 0..MAX_ITER {
        let hess = hessian(&func, &x, &eval)?;
        let step = newton_step(&hess, &grad)?;
        let slope = dot(&grad, &step);
        if !(slope < 0.0) {
            return Err(format!(
                "牛顿方向不是下降方向: g·p = {slope:.3e}，Hessian 可能不是正定的"
            ));
        }

        x = axpy(&x, 1.0, &step);
        eval = func(&x);
        if !eval.value.is_finite() {
            return Err(format!("牛顿步后函数值非有限: {}", eval.value));
        }
        grad = gradient(&func, &x, &eval)?;
        grad_norm = norm_inf(&grad);

        log::debug!(
            "牛顿法 第 {} 次迭代: F = {:.12}, |g|_inf = {:.3e}, |p|_inf = {:.3e}",
            iter + 1,
            eval.value,
            grad_norm,
            norm_inf(&step)
        );

        if grad_norm <= conv_thresh {
            return Ok(x);
        }
    }

    Err(format!(
        "牛顿法在 {MAX_ITER} 次迭代后未收敛，当前梯度无穷范数为 {grad_norm:.3e}（阈值 {conv_thresh:.3e}）"
    ))
}
