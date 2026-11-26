//! 阻尼牛顿法（Damped Newton Method）
//!
//! 在牛顿方向上进行 Armijo 回溯线搜索。若 Hessian 非正定，则对其进行特征值平移
//! `H + shift I`，使修正后的矩阵正定，从而保证牛顿方向是下降方向。
//! 相比纯牛顿法，阻尼牛顿法具有更好的全局收敛性。
//!
//! 收敛判据为梯度的无穷范数：`max_i |∂F/∂x_i| <= conv_thresh`。

use super::FuncEval;
use super::utils::{
    MAX_ITER, backtracking, dot, eig_hessian, gradient, hessian, neg, norm_inf, solve_shifted, validate_inputs,
};

/// 使用阻尼牛顿法最小化目标函数。
///
/// # 参数
/// - `func`: 目标函数，输入为自变量 `&Vec<f64>`，输出为 [`FuncEval`]。
///   若未提供梯度或 Hessian，则自动用有限差分近似。
/// - `paras`: 初始自变量。
/// - `conv_thresh`: 收敛阈值，当梯度的无穷范数小于该值时认为收敛。
///
/// # 返回
/// 收敛后的自变量。若线搜索失败或达到最大迭代次数仍未收敛，返回 `Err`。
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
        let eig = eig_hessian(&hess, x.len())?;
        let min_eig = eig.values.iter().cloned().fold(f64::INFINITY, f64::min);
        let scale = eig.values.iter().fold(0.0_f64, |max, v| max.max(v.abs())).max(1.0);

        // Hessian 非正定时做特征值平移，保证 H + shift I 正定
        let shift = if min_eig <= 0.0 { -min_eig + 1e-8 * scale } else { 0.0 };
        let mut dir = solve_shifted(&eig, &grad, shift, f64::MIN_POSITIVE)?;
        if dot(&grad, &dir) >= 0.0 {
            log::debug!("修正后的牛顿方向不是下降方向，回退到负梯度方向");
            dir = neg(&grad);
        }

        let line = backtracking(&func, &x, eval.value, &grad, &dir, 1.0)?;
        let new_grad = gradient(&func, &line.x, &line.eval)?;

        x = line.x;
        eval = line.eval;
        grad = new_grad;
        grad_norm = norm_inf(&grad);

        log::debug!(
            "阻尼牛顿法 第 {} 次迭代: F = {:.12}, |g|_inf = {:.3e}, alpha = {:.3e}, shift = {:.3e}",
            iter + 1,
            eval.value,
            grad_norm,
            line.alpha,
            shift
        );

        if grad_norm <= conv_thresh {
            return Ok(x);
        }
    }

    Err(format!(
        "阻尼牛顿法在 {MAX_ITER} 次迭代后未收敛，当前梯度无穷范数为 {grad_norm:.3e}（阈值 {conv_thresh:.3e}）"
    ))
}
