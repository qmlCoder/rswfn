//! 非线性共轭梯度法（Conjugate Gradient）
//!
//! 使用 Polak–Ribière+ 公式计算共轭方向，并用强 Wolfe 线搜索确定步长：
//!
//! ```text
//! beta = max(0, g_new·(g_new - g_old) / (g_old·g_old))
//! p_new = -g_new + beta * p_old
//! ```
//!
//! 当线搜索失败或新方向不是下降方向时，自动退化为负梯度方向并重新开始。
//! 收敛判据为梯度的无穷范数：`max_i |∂F/∂x_i| <= conv_thresh`。

use super::FuncEval;
use super::utils::{
    MAX_ITER, add_scaled, backtracking, dot, gradient, initial_step, neg, norm_inf, strong_wolfe, validate_inputs,
};

/// 使用非线性共轭梯度法最小化目标函数。
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

    // 初始方向取负梯度；restart 为 true 时下一步重新使用最速下降方向
    let mut dir = neg(&grad);
    let mut restart = true;

    for iter in 0..MAX_ITER {
        if restart {
            dir = neg(&grad);
        }
        // 保证搜索方向为下降方向
        if dot(&grad, &dir) >= 0.0 {
            dir = neg(&grad);
            restart = true;
        }

        let alpha0 = initial_step(&grad, &dir, eval.hessian.as_ref());
        let (line, used_fallback) = match strong_wolfe(&func, &x, eval.value, &grad, &dir, alpha0)? {
            Some(line) => (line, false),
            None => {
                log::debug!("强 Wolfe 线搜索失败，回退到 Armijo 回溯线搜索");
                (backtracking(&func, &x, eval.value, &grad, &dir, alpha0)?, true)
            }
        };

        let new_grad = match &line.grad {
            Some(grad) => grad.clone(),
            None => gradient(&func, &line.x, &line.eval)?,
        };
        let new_grad_norm = norm_inf(&new_grad);

        if new_grad_norm <= conv_thresh {
            log::debug!(
                "共轭梯度法 第 {} 次迭代: F = {:.12}, |g|_inf = {:.3e}, alpha = {:.3e}",
                iter + 1,
                line.eval.value,
                new_grad_norm,
                line.alpha
            );
            return Ok(line.x);
        }

        // Polak–Ribière+：负的 beta 截断为 0，等价于重启
        let gg_old = dot(&grad, &grad);
        let beta = if gg_old > 0.0 {
            let num = dot(&new_grad, &new_grad) - dot(&new_grad, &grad);
            (num / gg_old).max(0.0)
        } else {
            0.0
        };

        let mut new_dir = add_scaled(&neg(&new_grad), beta, &dir);
        let is_descent = dot(&new_grad, &new_dir) < 0.0;
        restart = used_fallback || beta <= 0.0 || !is_descent;
        if restart {
            new_dir = neg(&new_grad);
        }

        x = line.x;
        eval = line.eval;
        grad = new_grad;
        grad_norm = new_grad_norm;
        dir = new_dir;

        log::debug!(
            "共轭梯度法 第 {} 次迭代: F = {:.12}, |g|_inf = {:.3e}, alpha = {:.3e}, beta = {:.3e}",
            iter + 1,
            eval.value,
            grad_norm,
            line.alpha,
            beta
        );
    }

    Err(format!(
        "共轭梯度法在 {MAX_ITER} 次迭代后未收敛，当前梯度无穷范数为 {grad_norm:.3e}（阈值 {conv_thresh:.3e}）"
    ))
}
