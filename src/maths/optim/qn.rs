//! 拟牛顿法（Quasi-Newton Method，BFGS）
//!
//! 使用 BFGS 公式逐步构造逆 Hessian 的近似矩阵 `B`：
//!
//! ```text
//! s = x_new - x,  y = g_new - g,  rho = 1 / (y·s)
//! B_new = (I - rho s y^T) B (I - rho y s^T) + rho s s^T
//! ```
//!
//! 搜索方向为 `p = -B g`，使用强 Wolfe 线搜索确定步长；若强 Wolfe 条件
//! 未满足，则回退到 Armijo 回溯线搜索。该方法只需要梯度，不需要 Hessian。
//!
//! 收敛判据为梯度的无穷范数：`max_i |∂F/∂x_i| <= conv_thresh`。

use ndarray::Array2;

use super::FuncEval;
use super::utils::{
    MAX_ITER, backtracking, dot, gradient, initial_step, mat_vec, neg, norm_inf, strong_wolfe, sub, validate_inputs,
};

/// 使用拟牛顿法（BFGS）最小化目标函数。
///
/// # 参数
/// - `func`: 目标函数，输入为自变量 `&Vec<f64>`，输出为 [`FuncEval`]。
///   若未提供梯度，则自动用中心差分近似；Hessian 可选，仅用于估计线搜索
///   的初始步长。
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

    let dim = x.len();
    // 逆 Hessian 近似矩阵，初始取单位矩阵
    let mut inv_hessian = Array2::<f64>::eye(dim);

    for iter in 0..MAX_ITER {
        // 拟牛顿方向 p = -B g
        let mut dir = mat_vec(&inv_hessian, &grad);
        for value in dir.iter_mut() {
            *value = -*value;
        }

        // 数值误差可能破坏正定性，此时重置为单位矩阵
        if dot(&grad, &dir) >= 0.0 {
            log::debug!("拟牛顿方向不是下降方向，重置逆 Hessian 近似为单位矩阵");
            inv_hessian = Array2::<f64>::eye(dim);
            dir = neg(&grad);
        }

        let alpha0 = initial_step(&grad, &dir, eval.hessian.as_ref());
        let line = match strong_wolfe(&func, &x, eval.value, &grad, &dir, alpha0)? {
            Some(line) => line,
            None => {
                log::debug!("强 Wolfe 线搜索失败，回退到 Armijo 回溯线搜索");
                backtracking(&func, &x, eval.value, &grad, &dir, alpha0)?
            }
        };

        let new_grad = match &line.grad {
            Some(grad) => grad.clone(),
            None => gradient(&func, &line.x, &line.eval)?,
        };
        let new_grad_norm = norm_inf(&new_grad);

        if new_grad_norm <= conv_thresh {
            log::debug!(
                "拟牛顿法(BFGS) 第 {} 次迭代: F = {:.12}, |g|_inf = {:.3e}, alpha = {:.3e}",
                iter + 1,
                line.eval.value,
                new_grad_norm,
                line.alpha
            );
            return Ok(line.x);
        }

        let s = sub(&line.x, &x);
        let y = sub(&new_grad, &grad);
        let updated = bfgs_update(&mut inv_hessian, &s, &y);
        if !updated {
            log::debug!("BFGS 曲率条件不满足，跳过本次逆 Hessian 更新");
        }

        x = line.x;
        eval = line.eval;
        grad = new_grad;
        grad_norm = new_grad_norm;

        log::debug!(
            "拟牛顿法(BFGS) 第 {} 次迭代: F = {:.12}, |g|_inf = {:.3e}, alpha = {:.3e}, 更新 = {}",
            iter + 1,
            eval.value,
            grad_norm,
            line.alpha,
            updated
        );
    }

    Err(format!(
        "拟牛顿法(BFGS) 在 {MAX_ITER} 次迭代后未收敛，当前梯度无穷范数为 {grad_norm:.3e}（阈值 {conv_thresh:.3e}）"
    ))
}

/// BFGS 逆 Hessian 更新。
///
/// 当 `y·s` 不是显著正数时跳过更新（保持当前近似矩阵不变）。
/// 返回是否成功更新。
fn bfgs_update(inv_hessian: &mut Array2<f64>, s: &[f64], y: &[f64]) -> bool {
    let dim = s.len();
    if inv_hessian.nrows() != dim || inv_hessian.ncols() != dim {
        return false;
    }

    let s_norm = norm_inf(s);
    let y_norm = norm_inf(y);
    if s_norm == 0.0 || y_norm == 0.0 {
        return false;
    }

    let ys = dot(y, s);
    if !ys.is_finite() || ys <= 1e-12 * s_norm * y_norm {
        return false;
    }

    let rho = 1.0 / ys;
    let hy = mat_vec(inv_hessian, y);
    let yhy = dot(y, &hy);
    let coeff = rho * rho * yhy + rho;
    if !coeff.is_finite() {
        return false;
    }

    for i in 0..dim {
        for j in 0..dim {
            inv_hessian[(i, j)] += coeff * s[i] * s[j] - rho * (s[i] * hy[j] + hy[i] * s[j]);
        }
    }

    // 消除数值噪声造成的非对称性
    for i in 0..dim {
        for j in 0..i {
            let avg = 0.5 * (inv_hessian[(i, j)] + inv_hessian[(j, i)]);
            inv_hessian[(i, j)] = avg;
            inv_hessian[(j, i)] = avg;
        }
    }

    if inv_hessian.iter().any(|value| !value.is_finite()) {
        *inv_hessian = Array2::<f64>::eye(dim);
        return false;
    }

    true
}
