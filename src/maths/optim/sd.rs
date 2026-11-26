//! 最速下降法（Steepest Descent）
//!
//! 每步沿负梯度方向 `p = -g` 搜索，并使用 Armijo 回溯线搜索确定步长。
//! 收敛判据为梯度的无穷范数：`max_i |∂F/∂x_i| <= conv_thresh`。

use super::FuncEval;
use super::utils::{MAX_ITER, backtracking, gradient, initial_step, neg, norm_inf, validate_inputs};

/// 使用最速下降法最小化目标函数。
///
/// # 参数
/// - `func`: 目标函数，输入为自变量 `&Vec<f64>`，输出为 [`FuncEval`]。
///   若未提供梯度，则自动使用中心差分近似。
/// - `paras`: 初始自变量。
/// - `conv_thresh`: 收敛阈值，当梯度的无穷范数小于该值时认为收敛。
///
/// # 返回
/// 收敛后的自变量。若达到最大迭代次数仍未收敛，返回 `Err`。
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
        // 最速下降方向
        let dir = neg(&grad);
        let alpha0 = initial_step(&grad, &dir, eval.hessian.as_ref());
        let line = backtracking(&func, &x, eval.value, &grad, &dir, alpha0)?;

        let new_grad = gradient(&func, &line.x, &line.eval)?;
        x = line.x;
        eval = line.eval;
        grad = new_grad;
        grad_norm = norm_inf(&grad);

        log::debug!(
            "最速下降法 第 {} 次迭代: F = {:.12}, |g|_inf = {:.3e}, alpha = {:.3e}",
            iter + 1,
            eval.value,
            grad_norm,
            line.alpha
        );

        if grad_norm <= conv_thresh {
            return Ok(x);
        }
    }

    Err(format!(
        "最速下降法在 {MAX_ITER} 次迭代后未收敛，当前梯度无穷范数为 {grad_norm:.3e}（阈值 {conv_thresh:.3e}）"
    ))
}
