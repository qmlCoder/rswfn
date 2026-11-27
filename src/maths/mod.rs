//! 数学工具：矩阵运算、向量几何、优化与等值面提取。

pub mod march;
/// 矩阵相关运算
pub mod matrix;
pub mod optim;

use crate::maths;
use crate::{base::Vec3, utils};
use ndarray::{Array1, Array2, Axis};
use std::backtrace::Backtrace;
use std::f64::consts::PI;

/// 计算(重叠)矩阵的平方根 S^{1/2}
pub fn matrix_phalf(smat: &Array2<f64>) -> Array2<f64> {
    let (eig_val, vec_mat) = maths::eigh(smat);
    let nmat = smat.shape()[0];
    let mut val_mat = Array2::<f64>::zeros((nmat, nmat));

    for i in 0..nmat {
        let eig = eig_val[i];
        let eig = if eig < 0.0 { 0.0 } else { eig };
        val_mat[(i, i)] = eig.sqrt();
    }
    vec_mat.dot(&val_mat).dot(&vec_mat.t())
}

/// 计算(重叠)矩阵的负平方根 S^{-1/2}
pub fn matrix_nhalf(smat: &Array2<f64>) -> Array2<f64> {
    let (eig_val, vec_mat) = maths::eigh(smat);
    let nmat = smat.shape()[0];
    let mut val_mat = Array2::<f64>::zeros((nmat, nmat));
    for i in 0..nmat {
        let eig = eig_val[i];
        let eig = if eig < 0.0 { 0.0 } else { eig };
        val_mat[(i, i)] = 1.0 / eig.sqrt();
    }
    vec_mat.dot(&val_mat).dot(&vec_mat.t())
}

/// 根据系数矩阵和轨道占据轨道索引计算密度矩
/// - cmat: 系数矩阵
/// - idxs: 轨道占据索引
/// - open: 轨道是否为开壳层
pub fn calc_dmat(cmat: &Array2<f64>, obt_occs: &Vec<f64>) -> Array2<f64> {
    let nato = cmat.nrows();
    let nobt = cmat.ncols();
    let mut dmat = Array2::<f64>::zeros((nato, nato));
    for i in 0..nobt {
        let occ = obt_occs[i];
        if occ == 0.0 {
            continue;
        }
        let col = cmat.column(i);
        // 外积 col * col^T 作为密度矩阵贡献
        for a in 0..nato {
            for b in 0..nato {
                dmat[(a, b)] += col[a] * col[b] * occ;
            }
        }
    }
    dmat
}

/// 计算三维空间中两点之间的距离
pub fn calc_dist(p1: &[f64; 3], p2: &[f64; 3]) -> f64 {
    let [x1, y1, z1] = p1;
    let [x2, y2, z2] = p2;
    f64::sqrt((x2 - x1).powi(2) + (y2 - y1).powi(2) + (z2 - z1).powi(2))
}

/// 拉格朗日插值，根据已知数据点计算目标点的函数值
pub fn lag_intpol(xs: &Vec<f64>, ys: &Vec<f64>, tv: f64) -> f64 {
    let nx = xs.len();
    // let mut loc;
    let mut loc = 0; // 找到离当前点最近的点
    for i in 0..nx {
        loc = i;
        if xs[i] > tv {
            break;
        }
    }
    let up = 0.max(loc as isize - 3) as usize;
    let lo = (nx).min(loc + 2);
    if tv < xs[0] {
        let slop = (ys[1] - ys[0]) / (xs[1] - xs[0]);
        ys[0] + slop * (tv - xs[0])
    } else if tv > xs[nx - 1] {
        ys[nx - 1]
    } else {
        let mut iv: f64 = 0.0;
        for i in up..lo {
            let mut pv = 1.0;
            for j in up..lo {
                if i == j {
                    continue;
                }
                pv *= (tv - xs[j]) / (xs[i] - xs[j]);
            }
            iv += pv * ys[i];
        }
        iv
    }
}

/// 将向量v1投影到v2上，v2必须是单位向量
pub fn vec_project(vec1: [f64; 3], vec2: [f64; 3]) -> Result<[f64; 3], String> {
    let [x1, y1, z1] = vec1;
    let [x2, y2, z2] = vec2;
    let len2 = vec_len(&vec2);
    if (len2 - 1.0).abs() > 1e-2 {
        let msg = format!("v2 不是单位向量{:?},{}", vec2, len2);
        log::error!("{}", msg);
        return Err(msg);
    }
    let dot = x1 * x2 + y1 * y2 + z1 * z2;
    Ok([dot * x2, dot * y2, dot * z2])
}

/// 向量的归一化
pub fn vec_norm(vec: &[f64; 3]) -> [f64; 3] {
    let [x, y, z] = vec;
    let len = f64::sqrt(x.powi(2) + y.powi(2) + z.powi(2));
    [x / len, y / len, z / len]
}

/// 向量的长度
pub fn vec_len(vec: &[f64; 3]) -> f64 {
    let [x, y, z] = vec;
    f64::sqrt(x.powi(2) + y.powi(2) + z.powi(2))
}

/// 向量的距离
pub fn vec_dist(vec1: &[f64; 3], vec2: &[f64; 3]) -> f64 {
    let [x1, y1, z1] = *vec1;
    let [x2, y2, z2] = *vec2;
    f64::sqrt((x1 - x2).powi(2) + (y1 - y2).powi(2) + (z1 - z2).powi(2))
}

/// 计算向量的叉积
pub fn vec_cross(vec1: &[f64; 3], vec2: &[f64; 3]) -> [f64; 3] {
    let [x1, y1, z1] = vec1;
    let [x2, y2, z2] = vec2;
    [y1 * z2 - z1 * y2, z1 * x2 - x1 * z2, x1 * y2 - y1 * x2]
}

/// 计算一组数的标准差 总体标准差（Population Standard Deviation）
pub fn vec_psd(nums: &Vec<f64>) -> f64 {
    let mut sum = 0.0;
    for num in nums.iter() {
        sum += *num;
    }
    let mean = sum / nums.len() as f64;
    let mut sum = 0.0;
    for i in 0..nums.len() {
        sum += (nums[i] - mean).powi(2);
    }
    (sum / nums.len() as f64).sqrt()
}

/// 通过三个不共线的点确定一个平面，返回平面方程 Ax+By+Cz+D=0 的系数
pub fn get_plan_by_3points(p1: &Vec3, p2: &Vec3, p3: &Vec3) -> [f64; 4] {
    let v21 = p1 - p2;
    let v23 = p3 - p2;
    let [a, b, c] = v23.cross(&v21).into_arr();
    let [x, y, z] = p1.into_arr();
    let d = -a * x - b * y - c * z;
    [a, b, c, d]
}

/// 返回 a 的绝对值乘以 b 的符号
pub fn sign(a: f64, b: f64) -> f64 {
    if b > 0.0 { a.abs() } else { a.abs() * (-1.0) }
}

/// 根据能量和温度计算玻尔兹曼分布概率
/// - ens: 各能级的能量
/// - t: 温度（开尔文）
pub fn bolzman_dist(ens: &Vec<f64>, t: f64) -> Vec<f64> {
    let k = 1.380649e-23; // 玻尔兹曼常数
    let mut ps = vec![0.0; ens.len()];
    for (i, e) in ens.iter().enumerate() {
        ps[i] = (-e / (k * t)).exp();
    }
    let sum = ps.iter().sum::<f64>();
    for p in ps.iter_mut() {
        *p /= sum;
    }
    ps
}

/// 获取一个数组中的最大值及其索引
pub fn arg_max(vals: &Vec<f64>) -> (usize, f64) {
    let mut idx = 0;
    let mut max = f64::MIN;
    for (i, val) in vals.iter().enumerate() {
        if *val > max {
            idx = i;
            max = *val
        }
    }
    (idx, max)
}

/// 获取一个数组中的最小值及其索引
pub fn arg_min(vals: &Vec<f64>) -> (usize, f64) {
    let mut idx = 0;
    let mut min = f64::MAX;
    for (i, val) in vals.iter().enumerate() {
        if *val < min {
            idx = i;
            min = *val
        }
    }
    (idx, min)
}

/// 计算整数的阶乘 n!
pub fn fac(n: isize) -> isize {
    if n == 0 {
        return 1;
    }
    return n * fac(n - 1);
}

/// 计算整数的双阶乘 n!!
pub fn fac2(n: isize) -> isize {
    if n % 2 == 0 {
        // 偶数情况
        let k = n / 2;
        return 2_isize.pow(k as u32) * fac(k);
    } else {
        // 奇数情况
        let k = (n + 1) / 2;
        return fac(2 * k) / (2_isize.pow(k as u32) * fac(k));
    }
}

/// 对实对称矩阵进行特征值分解，返回特征值和特征向量
/// 使用 faer 后端
pub fn eigh(mat: &Array2<f64>) -> (Array1<f64>, Array2<f64>) {
    use faer::Side;
    let mat = into_faer(mat);
    let eigen = mat.self_adjoint_eigen(Side::Lower).unwrap();
    let S = eigen.S();
    let U = eigen.U();
    let nrow = mat.nrows();
    let mut eig_val = Array1::zeros(nrow);
    let mut eig_vec = Array2::zeros((nrow, nrow));
    for i in 0..nrow {
        eig_val[i] = S[i];
        for j in 0..nrow {
            eig_vec[(i, j)] = U[(i, j)];
        }
    }
    (eig_val, eig_vec)
}

/// 计算方阵的逆矩阵
/// 使用 faer 后端
pub fn inv(mat: &Array2<f64>) -> Array2<f64> {
    use faer::linalg::solvers::DenseSolveCore;
    let mat = into_faer(mat);
    let lu = mat.partial_piv_lu();
    let mat = lu.inverse();
    from_faer(&mat)
}

/// 计算方阵的行列式
/// 使用 faer 后端
pub fn det(mat: &Array2<f64>) -> f64 {
    let mat = into_faer(mat);
    mat.determinant()
}

/// 将 ndarray 矩阵转换为 faer 矩阵
fn into_faer(src: &Array2<f64>) -> faer::Mat<f64> {
    let nrow = src.nrows();
    let ncol = src.ncols();
    let mut out = faer::Mat::zeros(nrow, ncol);
    for i in 0..nrow {
        for j in 0..ncol {
            out[(i, j)] = src[(i, j)];
        }
    }
    out
}

/// 将 faer 矩阵转换为 ndarray 矩阵
fn from_faer(mat: &faer::Mat<f64>) -> Array2<f64> {
    let nrow = mat.nrows();
    let ncol = mat.ncols();
    let mut arr = Array2::zeros((nrow, ncol));
    for i in 0..nrow {
        for j in 0..ncol {
            arr[(i, j)] = mat[(i, j)];
        }
    }
    arr
}
