//! 优化算法的公共内部工具
//!
//! 包含数值梯度与 Hessian、线搜索、向量运算等。

use crate::maths;
use ndarray::Array2;

use super::FuncEval;

/// 默认最大迭代次数
pub(crate) const MAX_ITER: usize = 1000;

/// Armijo 线搜索的充分下降系数
const ARMIJO_C1: f64 = 1e-4;
/// 回溯线搜索的步长收缩因子
const BACKTRACK_RHO: f64 = 0.5;
/// 强 Wolfe 线搜索的曲率条件系数
const WOLFE_C2: f64 = 0.4;
/// 单次线搜索允许的最大迭代次数
const MAX_LINE_SEARCH_ITER: usize = 60;

/// 线性代数向量点积
#[inline]
pub(crate) fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// 向量的无穷范数
#[inline]
pub(crate) fn norm_inf(a: &[f64]) -> f64 {
    a.iter().fold(0.0_f64, |max, v| max.max(v.abs()))
}

/// 计算 `x + alpha * p`
#[inline]
pub(crate) fn axpy(x: &[f64], alpha: f64, p: &[f64]) -> Vec<f64> {
    x.iter().zip(p).map(|(xi, pi)| xi + alpha * pi).collect()
}

/// 计算 `a + alpha * b`
#[inline]
pub(crate) fn add_scaled(a: &[f64], alpha: f64, b: &[f64]) -> Vec<f64> {
    a.iter().zip(b).map(|(ai, bi)| ai + alpha * bi).collect()
}

/// 向量取负
#[inline]
pub(crate) fn neg(a: &[f64]) -> Vec<f64> {
    a.iter().map(|v| -v).collect()
}

/// 计算 `a - b`
#[inline]
pub(crate) fn sub(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b).map(|(ai, bi)| ai - bi).collect()
}

/// 矩阵与向量相乘
pub(crate) fn mat_vec(mat: &Array2<f64>, vec: &[f64]) -> Vec<f64> {
    let nrow = mat.nrows();
    let ncol = mat.ncols();
    let mut out = vec![0.0; nrow];
    for i in 0..nrow {
        for j in 0..ncol {
            out[i] += mat[(i, j)] * vec[j];
        }
    }
    out
}

/// 校验自变量和收敛阈值
pub(crate) fn validate_inputs(paras: &[f64], conv_thresh: f64) -> Result<(), String> {
    if paras.is_empty() {
        return Err("自变量 paras 不能为空".to_string());
    }
    if paras.iter().any(|v| !v.is_finite()) {
        return Err("自变量 paras 包含非有限值".to_string());
    }
    if !conv_thresh.is_finite() || conv_thresh <= 0.0 {
        return Err(format!("收敛阈值必须为有限正数，当前为 {conv_thresh}"));
    }
    Ok(())
}

/// 获取梯度：优先使用 [`FuncEval`] 中提供的解析梯度，否则用中心差分近似。
pub(crate) fn gradient<F>(func: &F, paras: &[f64], eval: &FuncEval) -> Result<Vec<f64>, String>
where
    F: Fn(&Vec<f64>) -> FuncEval,
{
    if let Some(grad) = &eval.grad {
        if grad.len() != paras.len() {
            return Err(format!("梯度维数 {} 与自变量维数 {} 不一致", grad.len(), paras.len()));
        }
        if grad.iter().any(|v| !v.is_finite()) {
            return Err("函数返回的梯度包含非有限值".to_string());
        }
        return Ok(grad.clone());
    }
    numeric_grad(func, paras)
}

/// 中心差分数值梯度
///
/// 步长取 `EPSILON^(1/3) * max(1, |x_i|)`，这是中心差分的近似最优步长。
pub(crate) fn numeric_grad<F>(func: &F, paras: &[f64]) -> Result<Vec<f64>, String>
where
    F: Fn(&Vec<f64>) -> FuncEval,
{
    let dim = paras.len();
    let mut grad = vec![0.0; dim];
    let mut trial = paras.to_vec();
    let base_step = f64::EPSILON.cbrt();

    for i in 0..dim {
        let step = base_step * paras[i].abs().max(1.0);
        trial[i] = paras[i] + step;
        let value_p = func(&trial).value;
        trial[i] = paras[i] - step;
        let value_m = func(&trial).value;
        trial[i] = paras[i];

        if !value_p.is_finite() || !value_m.is_finite() {
            return Err(format!("数值梯度计算失败：第 {i} 个分量处的函数值非有限"));
        }
        grad[i] = (value_p - value_m) / (2.0 * step);
    }

    if grad.iter().any(|v| !v.is_finite()) {
        return Err("数值梯度包含非有限值".to_string());
    }
    Ok(grad)
}

/// 获取 Hessian：优先使用 [`FuncEval`] 中提供的解析 Hessian，否则数值近似。
///
/// 若已有解析梯度，则对梯度做中心差分（每个分量只需两次函数求值）；
/// 否则对函数值做二阶中心差分。
pub(crate) fn hessian<F>(func: &F, paras: &[f64], eval: &FuncEval) -> Result<Array2<f64>, String>
where
    F: Fn(&Vec<f64>) -> FuncEval,
{
    let dim = paras.len();
    if let Some(hessian) = &eval.hessian {
        if hessian.nrows() != dim || hessian.ncols() != dim {
            return Err(format!(
                "Hessian 形状 {:?} 与自变量维数 {} 不一致",
                hessian.shape(),
                dim
            ));
        }
        if hessian.iter().any(|v| !v.is_finite()) {
            return Err("函数返回的 Hessian 包含非有限值".to_string());
        }
        return Ok(hessian.clone());
    }

    if eval.grad.is_some() {
        numeric_hessian_from_grad(func, paras)
    } else {
        numeric_hessian_from_value(func, paras, eval.value)
    }
}

/// 对解析梯度做中心差分得到 Hessian：`H[:, j] = [g(x + h e_j) - g(x - h e_j)] / (2h)`。
fn numeric_hessian_from_grad<F>(func: &F, paras: &[f64]) -> Result<Array2<f64>, String>
where
    F: Fn(&Vec<f64>) -> FuncEval,
{
    let dim = paras.len();
    let mut hess = Array2::<f64>::zeros((dim, dim));
    let base_step = f64::EPSILON.cbrt();
    let mut trial = paras.to_vec();

    for j in 0..dim {
        let step = base_step * paras[j].abs().max(1.0);

        trial[j] = paras[j] + step;
        let eval_p = func(&trial);
        let grad_p = gradient(func, &trial, &eval_p)?;

        trial[j] = paras[j] - step;
        let eval_m = func(&trial);
        let grad_m = gradient(func, &trial, &eval_m)?;

        trial[j] = paras[j];
        for i in 0..dim {
            hess[(i, j)] = (grad_p[i] - grad_m[i]) / (2.0 * step);
        }
    }

    // 消除数值噪声造成的非对称性
    for i in 0..dim {
        for j in 0..i {
            let avg = 0.5 * (hess[(i, j)] + hess[(j, i)]);
            hess[(i, j)] = avg;
            hess[(j, i)] = avg;
        }
    }

    if hess.iter().any(|v| !v.is_finite()) {
        return Err("数值 Hessian 包含非有限值".to_string());
    }
    Ok(hess)
}

/// 对函数值做二阶中心差分得到 Hessian。
fn numeric_hessian_from_value<F>(func: &F, paras: &[f64], value0: f64) -> Result<Array2<f64>, String>
where
    F: Fn(&Vec<f64>) -> FuncEval,
{
    let dim = paras.len();
    let mut hess = Array2::<f64>::zeros((dim, dim));
    let base_step = f64::EPSILON.powf(0.25);
    let steps: Vec<f64> = paras.iter().map(|v| base_step * v.abs().max(1.0)).collect();
    let mut trial = paras.to_vec();

    // 对角元
    for i in 0..dim {
        let step = steps[i];

        trial[i] = paras[i] + step;
        let value_p = func(&trial).value;

        trial[i] = paras[i] - step;
        let value_m = func(&trial).value;

        trial[i] = paras[i];
        if !value_p.is_finite() || !value_m.is_finite() {
            return Err(format!("数值 Hessian 计算失败：第 {i} 个对角元处的函数值非有限"));
        }
        hess[(i, i)] = (value_p - 2.0 * value0 + value_m) / (step * step);
    }

    // 非对角元
    for i in 0..dim {
        let step_i = steps[i];
        for j in (i + 1)..dim {
            let step_j = steps[j];

            trial[i] = paras[i] + step_i;
            trial[j] = paras[j] + step_j;
            let value_pp = func(&trial).value;

            trial[j] = paras[j] - step_j;
            let value_pm = func(&trial).value;

            trial[i] = paras[i] - step_i;
            trial[j] = paras[j] + step_j;
            let value_mp = func(&trial).value;

            trial[j] = paras[j] - step_j;
            let value_mm = func(&trial).value;

            trial[i] = paras[i];
            trial[j] = paras[j];
            if !(value_pp.is_finite() && value_pm.is_finite() && value_mp.is_finite() && value_mm.is_finite()) {
                return Err(format!("数值 Hessian 计算失败：分量 ({i}, {j}) 处的函数值非有限"));
            }

            let value = (value_pp - value_pm - value_mp + value_mm) / (4.0 * step_i * step_j);
            hess[(i, j)] = value;
            hess[(j, i)] = value;
        }
    }

    if hess.iter().any(|v| !v.is_finite()) {
        return Err("数值 Hessian 包含非有限值".to_string());
    }
    Ok(hess)
}

/// Hessian 的特征分解结果。
pub(crate) struct HessianEig {
    pub values: Vec<f64>,
    pub vectors: Array2<f64>,
}

/// 对 Hessian 做特征分解。
pub(crate) fn eig_hessian(hessian: &Array2<f64>, dim: usize) -> Result<HessianEig, String> {
    if hessian.nrows() != dim || hessian.ncols() != dim {
        return Err(format!(
            "Hessian 形状 {:?} 与自变量维数 {} 不一致",
            hessian.shape(),
            dim
        ));
    }
    if hessian.iter().any(|v| !v.is_finite()) {
        return Err("Hessian 包含非有限值".to_string());
    }

    let (values, vectors) = maths::eigh(hessian);
    let values = values.to_vec();
    if values.iter().any(|v| !v.is_finite()) {
        return Err("Hessian 特征值包含非有限值".to_string());
    }
    Ok(HessianEig { values, vectors })
}

/// 用特征分解求解 `(H + shift I) p = -g`。
pub(crate) fn solve_shifted(eig: &HessianEig, grad: &[f64], shift: f64, singular_tol: f64) -> Result<Vec<f64>, String> {
    let dim = grad.len();
    if eig.values.len() != dim {
        return Err("Hessian 特征值个数与梯度维数不一致".to_string());
    }

    let mut step = vec![0.0; dim];
    for k in 0..dim {
        let denom = eig.values[k] + shift;
        if !denom.is_finite() || denom.abs() <= singular_tol {
            return Err(format!("修正后的 Hessian 接近奇异：第 {k} 个特征值为 {denom:.3e}"));
        }

        let mut proj = 0.0;
        for i in 0..dim {
            proj += eig.vectors[(i, k)] * grad[i];
        }
        let coeff = -proj / denom;
        for i in 0..dim {
            step[i] += coeff * eig.vectors[(i, k)];
        }
    }

    if step.iter().any(|v| !v.is_finite()) {
        return Err("牛顿步包含非有限值".to_string());
    }
    Ok(step)
}

/// 计算未加阻尼的牛顿步 `p = -H^{-1} g`。
pub(crate) fn newton_step(hessian: &Array2<f64>, grad: &[f64]) -> Result<Vec<f64>, String> {
    let eig = eig_hessian(hessian, grad.len())?;
    let scale = eig.values.iter().fold(0.0_f64, |max, v| max.max(v.abs())).max(1.0);
    solve_shifted(&eig, grad, 0.0, 1e-12 * scale)
}

/// 线搜索的返回结果
#[derive(Debug)]
pub(crate) struct LineSearchResult {
    /// 使用的步长
    pub alpha: f64,
    /// 新的自变量
    pub x: Vec<f64>,
    /// 新自变量处的函数求值结果
    pub eval: FuncEval,
    /// 新自变量处的梯度（如果线搜索过程中已经算过）
    pub grad: Option<Vec<f64>>,
}

/// 保证初始步长为有限正数
pub(crate) fn sanitize_alpha(alpha: f64) -> f64 {
    if alpha.is_finite() && alpha > 0.0 { alpha } else { 1.0 }
}

/// 根据 Hessian 估计二次模型下的初始步长。
///
/// 对搜索方向 `p`，沿 `p` 的二次模型极小点为
/// `alpha = -(g·p) / (p^T H p)`。仅当该值为有限正数且曲率为正时采用，
/// 否则回退到 `1.0`。
pub(crate) fn initial_step(grad: &[f64], dir: &[f64], hessian: Option<&Array2<f64>>) -> f64 {
    if let Some(hessian) = hessian {
        let dim = grad.len();
        if hessian.nrows() == dim && hessian.ncols() == dim {
            let hdir = mat_vec(hessian, dir);
            let curvature = dot(dir, &hdir);
            let slope = dot(grad, dir);
            if slope < 0.0 && curvature > 0.0 {
                let alpha = -slope / curvature;
                if alpha.is_finite() && alpha > 0.0 {
                    return alpha;
                }
            }
        }
    }
    1.0
}

/// Armijo 回溯线搜索。
///
/// 从 `alpha0` 开始，不断将步长乘以 `BACKTRACK_RHO`，直到满足
/// `F(x + alpha p) <= F(x) + c1 * alpha * g·p`。
pub(crate) fn backtracking<F>(
    func: &F,
    x: &[f64],
    f0: f64,
    grad: &[f64],
    dir: &[f64],
    alpha0: f64,
) -> Result<LineSearchResult, String>
where
    F: Fn(&Vec<f64>) -> FuncEval,
{
    let slope = dot(grad, dir);
    if !(slope < 0.0) {
        return Err(format!("搜索方向不是下降方向: g·p = {slope:.3e}"));
    }

    let mut alpha = sanitize_alpha(alpha0);
    for _ in 0..MAX_LINE_SEARCH_ITER {
        let trial_x = axpy(x, alpha, dir);
        let eval = func(&trial_x);
        if eval.value.is_finite() && eval.value <= f0 + ARMIJO_C1 * alpha * slope {
            return Ok(LineSearchResult {
                alpha,
                x: trial_x,
                eval,
                grad: None,
            });
        }
        alpha *= BACKTRACK_RHO;
        if !alpha.is_finite() || alpha < 1e-20 {
            break;
        }
    }

    Err("回溯线搜索失败：未能满足 Armijo 充分下降条件".to_string())
}

/// 强 Wolfe 线搜索。
///
/// 同时满足以下两个条件：
/// 1. `F(x + alpha p) <= F(x) + c1 * alpha * g·p`
/// 2. `|g(x + alpha p)·p| <= c2 * |g·p|`
///
/// 若在最大迭代次数内未找到满足条件的步长，返回 `Ok(None)`，
/// 调用方可回退到 [`backtracking`]。
pub(crate) fn strong_wolfe<F>(
    func: &F,
    x: &[f64],
    f0: f64,
    grad: &[f64],
    dir: &[f64],
    alpha0: f64,
) -> Result<Option<LineSearchResult>, String>
where
    F: Fn(&Vec<f64>) -> FuncEval,
{
    let slope0 = dot(grad, dir);
    if !(slope0 < 0.0) {
        return Ok(None);
    }
    let curvature_tol = WOLFE_C2 * slope0.abs();

    let mut alpha = sanitize_alpha(alpha0);
    let mut alpha_prev = 0.0;
    let mut phi_prev = f0;
    let mut slope_prev = slope0;

    for _ in 0..MAX_LINE_SEARCH_ITER {
        let trial_x = axpy(x, alpha, dir);
        let eval = func(&trial_x);
        let phi = eval.value;

        let armijo_failed = phi > f0 + ARMIJO_C1 * alpha * slope0;
        if !phi.is_finite() || armijo_failed || (alpha_prev > 0.0 && phi >= phi_prev) {
            return zoom(
                func,
                x,
                f0,
                dir,
                alpha_prev,
                phi_prev,
                slope_prev,
                alpha,
                phi,
                slope0,
                curvature_tol,
            );
        }

        let trial_grad = gradient(func, &trial_x, &eval)?;
        let slope = dot(&trial_grad, dir);
        if !slope.is_finite() {
            return Ok(None);
        }
        if slope.abs() <= curvature_tol {
            return Ok(Some(LineSearchResult {
                alpha,
                x: trial_x,
                eval,
                grad: Some(trial_grad),
            }));
        }
        if slope >= 0.0 {
            return zoom(
                func,
                x,
                f0,
                dir,
                alpha,
                phi,
                slope,
                alpha_prev,
                phi_prev,
                slope0,
                curvature_tol,
            );
        }

        alpha_prev = alpha;
        phi_prev = phi;
        slope_prev = slope;
        alpha *= 2.0;
        if !alpha.is_finite() || alpha > 1e20 {
            break;
        }
    }

    Ok(None)
}

/// 强 Wolfe 线搜索的区间收缩过程。
#[allow(clippy::too_many_arguments)]
fn zoom<F>(
    func: &F,
    x: &[f64],
    f0: f64,
    dir: &[f64],
    mut alpha_lo: f64,
    mut phi_lo: f64,
    mut slope_lo: f64,
    mut alpha_hi: f64,
    mut phi_hi: f64,
    slope0: f64,
    curvature_tol: f64,
) -> Result<Option<LineSearchResult>, String>
where
    F: Fn(&Vec<f64>) -> FuncEval,
{
    for _ in 0..MAX_LINE_SEARCH_ITER {
        let alpha = interpolate_alpha(alpha_lo, phi_lo, slope_lo, alpha_hi, phi_hi);
        let trial_x = axpy(x, alpha, dir);
        let eval = func(&trial_x);
        let phi = eval.value;

        if !phi.is_finite() || phi > f0 + ARMIJO_C1 * alpha * slope0 || phi >= phi_lo {
            alpha_hi = alpha;
            phi_hi = phi;
        } else {
            let trial_grad = gradient(func, &trial_x, &eval)?;
            let slope = dot(&trial_grad, dir);
            if !slope.is_finite() {
                return Ok(None);
            }
            if slope.abs() <= curvature_tol {
                return Ok(Some(LineSearchResult {
                    alpha,
                    x: trial_x,
                    eval,
                    grad: Some(trial_grad),
                }));
            }
            if slope * (alpha_hi - alpha_lo) >= 0.0 {
                alpha_hi = alpha_lo;
                phi_hi = phi_lo;
            }
            alpha_lo = alpha;
            phi_lo = phi;
            slope_lo = slope;
        }

        if (alpha_hi - alpha_lo).abs() <= 1e-16 * (1.0 + alpha.abs()) {
            break;
        }
    }

    Ok(None)
}

/// 用二次插值估计区间内的极小点，失败时退化为中点。
fn interpolate_alpha(lo: f64, phi_lo: f64, slope_lo: f64, hi: f64, phi_hi: f64) -> f64 {
    let a_min = lo.min(hi);
    let a_max = lo.max(hi);
    let width = a_max - a_min;
    if !width.is_finite() || width <= 0.0 {
        return lo;
    }

    let delta = hi - lo;
    let denom = 2.0 * (phi_hi - phi_lo - slope_lo * delta);
    let candidate = if phi_hi.is_finite() && denom.is_finite() && denom.abs() > f64::MIN_POSITIVE {
        lo - slope_lo * delta * delta / denom
    } else {
        f64::NAN
    };

    let margin = 1e-3 * width;
    if candidate.is_finite() && width > 2.0 * margin {
        candidate.clamp(a_min + margin, a_max - margin)
    } else {
        0.5 * (lo + hi)
    }
}
