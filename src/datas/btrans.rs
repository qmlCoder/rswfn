use cached::proc_macro::cached;
use ndarray::{Array2, s};
use std::sync::LazyLock;

use crate::maths;
// 数据来源于卢天博客 http://sobereva.com/97

// const D00:f64=
/// 球谐 d 函数到笛卡尔 d 函数的变换矩阵（6×5）
#[rustfmt::skip]
pub static DMAT: LazyLock<Array2<f64>> = LazyLock::new(|| {
    let mut dmat = Array2::<f64>::zeros((6,5));
    dmat[( 0, 0)]=-0.5;
    dmat[( 0, 3)]= f64::sqrt(3.0)/2.0;
    dmat[( 1, 0)]=-0.5;
    dmat[( 1, 3)]=-f64::sqrt(3.0)/2.0;
    dmat[( 2, 0)]= 1.0;
    dmat[( 3, 4)]= 1.0;
    dmat[( 4, 1)]= 1.0;
    dmat[( 5, 2)]= 1.0;
    dmat
});

/// 球谐 f 函数到笛卡尔 f 函数的变换矩阵（10×7）
#[rustfmt::skip]
pub static FMAT: LazyLock<Array2<f64>> = LazyLock::new(|| {
    let mut fmat = Array2::<f64>::zeros((10,7));
    fmat[( 0, 1)]=-f64::sqrt(3.0/8.0);
    fmat[( 0, 5)]= f64::sqrt(5.0/8.0);
    fmat[( 1, 2)]=-f64::sqrt(3.0/8.0);
    fmat[( 1, 6)]=-f64::sqrt(5.0/8.0);
    fmat[( 2, 0)]= 1.0;
    fmat[( 3, 1)]=-f64::sqrt(3.0/40.0);
    fmat[( 3, 5)]=-3.0/f64::sqrt(8.0) ;
    fmat[( 4, 2)]=-f64::sqrt(3.0/40.0);
    fmat[( 4, 6)]= 3.0/f64::sqrt(8.0);
    fmat[( 5, 0)]=-3.0/2.0/f64::sqrt(5.0);
    fmat[( 5, 3)]= f64::sqrt(3.0)/2.0;
    fmat[( 6, 1)]= f64::sqrt(6.0/5.0);
    fmat[( 7, 2)]= f64::sqrt(6.0/5.0);
    fmat[( 8, 0)]=-3.0/2.0/f64::sqrt(5.0);
    fmat[( 8, 3)]=-f64::sqrt(3.0)/2.0;
    fmat[( 9, 4)]= 1.0;
    fmat
});

/// 球谐 g 函数到笛卡尔 g 函数的变换矩阵（15×9）
#[rustfmt::skip]
pub static GMAT: LazyLock<Array2<f64>> = LazyLock::new(|| {
    let mut gmat = Array2::<f64>::zeros((15, 9));

    gmat[( 0, 0)] =  1.0;
    gmat[( 1, 2)] =  2.0 * f64::sqrt(5.0 / 14.0);
    gmat[( 2, 0)] = -3.0 * f64::sqrt(3.0 / 35.0);
    gmat[( 2, 3)] = -3.0 * f64::sqrt(3.0 / 28.0);
    gmat[( 3, 2)] = -3.0 / 2.0 * f64::sqrt(5.0 / 14.0);
    gmat[( 3, 6)] = -f64::sqrt(5.0 / 8.0);
    gmat[( 4, 0)] =  3.0 / 8.0;
    gmat[( 4, 3)] =  f64::sqrt(5.0) / 4.0;
    gmat[( 4, 7)] =  f64::sqrt(35.0) / 8.0;
    gmat[( 5, 1)] =  2.0 * f64::sqrt(5.0 / 14.0);
    gmat[( 6, 4)] =  3.0 / f64::sqrt(7.0);
    gmat[( 7, 1)] = -3.0 / 2.0 / f64::sqrt(14.0);
    gmat[( 7, 5)] = -3.0 / f64::sqrt(8.0);
    gmat[( 8, 4)] = -f64::sqrt(5.0 / 28.0);
    gmat[( 8, 8)] = -f64::sqrt(5.0)/ 2.0;
    gmat[( 9, 0)] = -3.0 * f64::sqrt(3.0 / 35.0);
    gmat[( 9, 3)] =  3.0 * f64::sqrt(3.0 / 28.0);
    gmat[(10, 2)] = -3.0 / 2.0 / f64::sqrt(14.0);
    gmat[(10, 6)] =  3.0 / f64::sqrt(8.0);
    gmat[(11, 0)] =  3.0 / 4.0 * f64::sqrt(3.0 / 35.0);
    gmat[(11, 7)] = -3.0 / 4.0 * f64::sqrt(3.0);
    gmat[(12, 1)] = -3.0 / 2.0 * f64::sqrt(5.0 / 14.0);
    gmat[(12, 5)] =  f64::sqrt(5.0 / 8.0);
    gmat[(13, 4)] = -f64::sqrt(5.0 / 28.0);
    gmat[(13, 8)] =  f64::sqrt(5.0) / 2.0;
    gmat[(14, 0)] =  3.0 / 8.0;
    gmat[(14, 3)] = -f64::sqrt(5.0) / 4.0;
    gmat[(14, 7)] =  f64::sqrt(35.0) / 8.0;

    gmat
});

/// 球谐 h 函数到笛卡尔 h 函数的变换矩阵（21×11）
#[rustfmt::skip]
pub static HMAT: LazyLock<Array2<f64>> = LazyLock::new(|| {
    let mut hmat = Array2::<f64>::zeros((21, 11));

    hmat[( 0, 0)] = 1.0;
    hmat[( 1, 2)] = f64::sqrt(5.0 / 3.0);
    hmat[( 2, 0)] = -5.0 / f64::sqrt(21.0);
    hmat[( 2, 3)] = -f64::sqrt(5.0 / 2.0);
    hmat[( 3, 2)] = -3.0 * f64::sqrt(5.0 / 28.0);
    hmat[( 3, 6)] = -f64::sqrt(5.0 / 6.0);
    hmat[( 4, 0)] = 5.0 / 8.0;
    hmat[( 4, 3)] = f64::sqrt(35.0 / 3.0) / 4.0;
    hmat[( 4, 7)] = f64::sqrt(35.0 / 8.0);
    hmat[( 5, 2)] = f64::sqrt(15.0 / 8.0);
    hmat[( 5, 6)] = f64::sqrt(35.0 / 2.0) / 8.0;
    hmat[( 5,10)] = 3.0 / 8.0 * f64::sqrt(7.0 / 2.0);
    hmat[( 6, 1)] = f64::sqrt(5.0 / 3.0);
    hmat[( 7, 4)] = f64::sqrt(5.0 / 3.0);
    hmat[( 8, 1)] = -3.0 / f64::sqrt(28.0);
    hmat[( 8, 5)] = -f64::sqrt(3.0 / 2.0);
    hmat[( 9, 4)] = -f64::sqrt(5.0 / 12.0);
    hmat[( 9, 8)] = -f64::sqrt(5.0 / 2.0);
    hmat[(10, 1)] = f64::sqrt(5.0 / 3.0) / 8.0;
    hmat[(10, 5)] = f64::sqrt(35.0 / 2.0) / 8.0;
    hmat[(10, 9)] = 5.0 / 8.0 * f64::sqrt(7.0 / 2.0);
    hmat[(11, 0)] = -5.0 / f64::sqrt(21.0);
    hmat[(11, 3)] = f64::sqrt(5.0 / 2.0);
    hmat[(12, 2)] = -3.0 / f64::sqrt(28.0);
    hmat[(12, 6)] = f64::sqrt(3.0 / 2.0);
    hmat[(13, 0)] = f64::sqrt(15.0 / 7.0) / 4.0;
    hmat[(13, 7)] = -3.0 / 4.0 * f64::sqrt(3.0);
    hmat[(14, 2)] = f64::sqrt(5.0 / 7.0) / 4.0;
    hmat[(14, 6)] = -f64::sqrt(5.0 / 6.0) / 4.0;
    hmat[(14,10)] = -5.0 / 4.0 * f64::sqrt(3.0 / 2.0);
    hmat[(15, 1)] = -3.0 * f64::sqrt(5.0 / 28.0);
    hmat[(15, 5)] = f64::sqrt(5.0 / 6.0);
    hmat[(16, 4)] = -f64::sqrt(5.0 / 12.0);
    hmat[(16, 8)] = f64::sqrt(5.0 / 2.0);
    hmat[(17, 1)] = f64::sqrt(5.0 / 7.0) / 4.0;
    hmat[(17, 5)] = f64::sqrt(5.0 / 6.0) / 4.0;
    hmat[(17, 9)] = -5.0 / 4.0 * f64::sqrt(3.0 / 2.0);
    hmat[(18, 0)] = 5.0 / 8.0;
    hmat[(18, 3)] = -f64::sqrt(35.0 / 3.0) / 4.0;
    hmat[(18, 7)] = f64::sqrt(35.0 / 8.0);
    hmat[(19, 2)] = f64::sqrt(5.0 / 3.0) / 8.0;
    hmat[(19, 6)] = -f64::sqrt(35.0 / 2.0) / 8.0;
    hmat[(19,10)] = 5.0 / 8.0 * f64::sqrt(7.0 / 2.0);
    hmat[(20, 1)] = f64::sqrt(15.0 / 8.0);
    hmat[(20, 5)] = -f64::sqrt(35.0 / 2.0) / 8.0;
    hmat[(20, 9)] = 3.0 / 8.0 * f64::sqrt(7.0 / 2.0);

    hmat
});

/// 球谐 d 到笛卡尔 d 变换矩阵的转置（5×6）
static DMAT_T: LazyLock<Array2<f64>> = LazyLock::new(|| DMAT.t().to_owned());
/// 球谐 f 到笛卡尔 f 变换矩阵的转置（7×10）
static FMAT_T: LazyLock<Array2<f64>> = LazyLock::new(|| FMAT.t().to_owned());
/// 球谐 g 到笛卡尔 g 变换矩阵的转置（9×15）
static GMAT_T: LazyLock<Array2<f64>> = LazyLock::new(|| GMAT.t().to_owned());
/// 球谐 h 到笛卡尔 h 变换矩阵的转置（11×21）
static HMAT_T: LazyLock<Array2<f64>> = LazyLock::new(|| HMAT.t().to_owned());

/// 根据角动量获取球谐到笛卡尔的变换矩阵
pub fn sph2car_tmat(ang: usize) -> Array2<f64> {
    match ang {
        0 => Array2::eye(1),
        1 => Array2::eye(3),
        2 => DMAT.to_owned(),
        3 => FMAT.to_owned(),
        4 => GMAT.to_owned(),
        5 => HMAT.to_owned(),
        _ => {
            dbg!("非法角动量:{}", ang);
            panic!("非法角动量:{}", ang)
        }
    }
}

/// 生成混合 car/sph 基组之间的系数变换矩阵。
///
/// 返回 `n_into × n_from` 的矩阵 `M`，满足 `cmat_into = M.dot(&cmat_from)`。
///
/// - 所有壳层类型不变：返回单位矩阵；
/// - 只发生收缩（`car -> sph`，可能夹杂不变壳层）：用完整 `from_smat` 做全局
///   度量修正投影 `M = (F^T S F)^{-1} F^T S`，与 `car2sph_fmat` 一致，能正确处理
///   笛卡尔 d/f 中 s 污染分量与 s 轨道的耦合；
/// - 只发生扩张（`sph -> car`，可能夹杂不变壳层）：返回嵌入矩阵的转置 `F^T`
///   （`sph -> car` 即 `U`）；
/// - 同时出现收缩和扩张：回退到逐壳层处理。
pub fn get_fmat(from_smat: &Array2<f64>, from_forms: &Vec<&str>, into_forms: &Vec<&str>, angs: &Vec<usize>) -> Array2<f64> {
    assert_eq!(from_forms.len(), into_forms.len());
    assert_eq!(from_forms.len(), angs.len());

    let mut n_from = 0;
    let mut n_into = 0;
    for (i, ang) in angs.iter().enumerate() {
        n_from += ang2syms(*ang, from_forms[i]).len();
        n_into += ang2syms(*ang, into_forms[i]).len();
    }

    let has_contract = (0..angs.len()).any(|i| from_forms[i] == "car" && into_forms[i] == "sph");
    let has_expand = (0..angs.len()).any(|i| from_forms[i] == "sph" && into_forms[i] == "car");

    // 类型完全不变
    if !has_contract && !has_expand {
        return Array2::eye(n_from);
    }

    // 构建 into -> from 的嵌入矩阵 F (n_from × n_into)
    let mut fmat_embed = Array2::zeros((n_from, n_into));
    let mut i_from = 0;
    let mut i_into = 0;
    for (i, ang) in angs.iter().enumerate() {
        let from_form = from_forms[i];
        let into_form = into_forms[i];
        let nf = ang2syms(*ang, from_form).len();
        let ni = ang2syms(*ang, into_form).len();

        let block = if from_form == into_form {
            Array2::eye(nf)
        } else if from_form == "car" && into_form == "sph" {
            sph2car_tmat(*ang) // U: car = U * sph
        } else if from_form == "sph" && into_form == "car" {
            car2sph_tmat(*ang) // U^T: sph = U^T * car
        } else {
            panic!("未知的轨道形式: {} -> {}", from_form, into_form);
        };

        fmat_embed.slice_mut(s![i_from..i_from + nf, i_into..i_into + ni]).assign(&block);
        i_from += nf;
        i_into += ni;
    }

    // 纯收缩：全局度量修正投影（等价于 car2sph_fmat）
    if has_contract && !has_expand {
        let metric = fmat_embed.t().dot(from_smat).dot(&fmat_embed);
        let metric_inv = maths::inv(&metric);
        return metric_inv.dot(&fmat_embed.t()).dot(from_smat);
    }

    // 纯扩张：直接使用嵌入转置
    if has_expand && !has_contract {
        return fmat_embed.t().to_owned();
    }

    // 混合收缩+扩张：逐壳层处理（近似）
    let mut fmat = Array2::zeros((n_into, n_from));
    let mut i_from = 0;
    let mut i_into = 0;
    for (i, ang) in angs.iter().enumerate() {
        let from_form = from_forms[i];
        let into_form = into_forms[i];
        let nf = ang2syms(*ang, from_form).len();
        let ni = ang2syms(*ang, into_form).len();

        let block = if from_form == into_form {
            Array2::eye(nf)
        } else if from_form == "car" && into_form == "sph" {
            let u = sph2car_tmat(*ang);
            let sblk = from_smat.slice(s![i_from..i_from + nf, i_from..i_from + nf]).to_owned();
            let metric = u.t().dot(&sblk).dot(&u);
            let metric_inv = maths::inv(&metric);
            metric_inv.dot(&u.t()).dot(&sblk)
        } else if from_form == "sph" && into_form == "car" {
            sph2car_tmat(*ang)
        } else {
            panic!("未知的轨道形式: {} -> {}", from_form, into_form);
        };

        fmat.slice_mut(s![i_into..i_into + ni, i_from..i_from + nf]).assign(&block);
        i_from += nf;
        i_into += ni;
    }
    fmat
}

/// 根据角动量获取笛卡尔到球谐的变换矩阵（即 `sph2car_tmat` 的转置）
///
/// # Parameters
/// - `ang`: 角动量量子数（0-5 分别对应 s、p、d、f、g、h）
///
/// # Returns
/// 对应的笛卡尔到球谐变换矩阵
///
/// 注意：该函数不能作用于分子轨道系数，只能作用于重叠积分等对称矩阵
pub fn car2sph_tmat(ang: usize) -> Array2<f64> {
    let tmat = sph2car_tmat(ang);
    tmat.t().to_owned()
}

/// 获取笛卡尔到球谐变换矩阵的静态引用（l >= 2）
///
/// 与 `car2sph_tmat` 数值相同，但不做克隆，供格点求值等热循环使用。
pub fn car2sph_tmat_ref(ang: usize) -> &'static Array2<f64> {
    match ang {
        2 => &DMAT_T,
        3 => &FMAT_T,
        4 => &GMAT_T,
        5 => &HMAT_T,
        _ => panic!("非法角动量:{}", ang),
    }
}

/// 根据角动量获取笛卡尔到球谐的变换矩阵
pub fn car2sph_fmat(smat_car: &Array2<f64>, angs: &Vec<usize>) -> Array2<f64> {
    let ncar: usize = angs.iter().map(|ang| ang2syms(*ang, "car").len()).sum();
    let nsph: usize = angs.iter().map(|ang| ang2syms(*ang, "sph").len()).sum();
    let mut irow = 0;
    let mut icol = 0;
    let mut fmat = Array2::zeros((ncar, nsph));
    for ang in angs {
        let tmat = match ang {
            0 => Array2::eye(1),
            1 => Array2::eye(3),
            2 => DMAT.to_owned(),
            3 => FMAT.to_owned(),
            4 => GMAT.to_owned(),
            5 => HMAT.to_owned(),
            _ => {
                dbg!("非法角动量:{}", ang);
                panic!("非法角动量:{}", ang)
            }
        };
        let nrow = tmat.nrows();
        let ncol = tmat.ncols();
        fmat.slice_mut(s![irow..irow + nrow, icol..icol + ncol]).assign(&tmat);
        irow += nrow;
        icol += ncol;
    }
    let inv = fmat.t().dot(smat_car).dot(&fmat);
    let inv = maths::inv(&inv);
    let fmat = inv.dot(&fmat.t()).dot(smat_car);
    fmat
}

/// s 壳层（混合类型）的轨道符号
pub const MIX_S_SYMS: [&str; 1] = ["S"];
/// s 壳层（混合类型）的角动量分量 (l, m, n)
pub const MIX_S_LMNS: [[usize; 3]; 1] = [[0, 0, 0]];

/// p 壳层（混合类型）的轨道符号
pub const MIX_P_SYMS: [&str; 3] = ["PX", "PY", "PZ"];
/// p 壳层（混合类型）的角动量分量 (l, m, n)
pub const MIX_P_LMNS: [[usize; 3]; 3] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

/// 笛卡尔 d 壳层的轨道符号（XX、YY、ZZ、XY、XZ、YZ）
pub const CAR_D_SYMS: [&str; 6] = [
    "XX", // [2,0,0]
    "YY", // [0,2,0]
    "ZZ", // [0,0,2]
    "XY", // [1,1,0]
    "XZ", // [1,0,1]
    "YZ", // [0,1,1]
];

/// 笛卡尔 d 壳层的角动量分量 (l, m, n)
pub const CAR_D_LMNS: [[usize; 3]; 6] = [[2, 0, 0], [0, 2, 0], [0, 0, 2], [1, 1, 0], [1, 0, 1], [0, 1, 1]];

/// 笛卡尔 f 壳层的轨道符号（共 10 个）
pub const CAR_F_SYMS: [&str; 10] = [
    "XXX", // [3,0,0] 0
    "YYY", // [0,3,0] 1
    "ZZZ", // [0,0,3] 2
    "XYY", // [1,2,0] 3
    "XXY", // [2,1,0] 4
    "XXZ", // [2,0,1] 5
    "XZZ", // [1,0,2] 6
    "YZZ", // [0,1,2] 7
    "YYZ", // [0,2,1] 8
    "XYZ", // [1,1,1] 9
];

/// 笛卡尔 f 壳层的角动量分量 (l, m, n)
pub const CAR_F_LMNS: [[usize; 3]; 10] = [[3, 0, 0], [0, 3, 0], [0, 0, 3], [1, 2, 0], [2, 1, 0], [2, 0, 1], [1, 0, 2], [0, 1, 2], [0, 2, 1], [1, 1, 1]];

/// 笛卡尔 g 壳层的轨道符号（共 15 个）
pub const CAR_G_SYMS: [&str; 15] = [
    "ZZZZ", // [0,0,4]
    "YZZZ", // [0,1,3]
    "YYZZ", // [0,2,2]
    "YYYZ", // [0,3,1]
    "YYYY", // [0,4,0]
    "XZZZ", // [1,0,3]
    "XYZZ", // [1,1,2]
    "XYYZ", // [1,2,1]
    "XYYY", // [1,3,0]
    "XXZZ", // [2,0,2]
    "XXYZ", // [2,1,1]
    "XXYY", // [2,2,0]
    "XXXZ", // [3,0,1]
    "XXXY", // [3,1,0]
    "XXXX", // [4,0,0]
];

/// 笛卡尔 g 壳层的角动量分量 (l, m, n)
pub const CAR_G_LMNS: [[usize; 3]; 15] = [
    [0, 0, 4],
    [0, 1, 3],
    [0, 2, 2],
    [0, 3, 1],
    [0, 4, 0],
    [1, 0, 3],
    [1, 1, 2],
    [1, 2, 1],
    [1, 3, 0],
    [2, 0, 2],
    [2, 1, 1],
    [2, 2, 0],
    [3, 0, 1],
    [3, 1, 0],
    [4, 0, 0],
];

/// 笛卡尔 h 壳层的轨道符号（共 21 个）
pub const CAR_H_SYMS: [&str; 21] = [
    "ZZZZZ", // [0,0,5]
    "YZZZZ", // [0,1,4]
    "YYZZZ", // [0,2,3]
    "YYYZZ", // [0,3,2]
    "YYYYZ", // [0,4,1]
    "YYYYY", // [0,5,0]
    "XZZZZ", // [1,0,4]
    "XYZZZ", // [1,1,3]
    "XYYZZ", // [1,2,2]
    "XYYYZ", // [1,3,1]
    "XYYYY", // [1,4,0]
    "XXZZZ", // [2,0,3]
    "XXYZZ", // [2,1,2]
    "XXYYZ", // [2,2,1]
    "XXYYY", // [2,3,0]
    "XXXZZ", // [3,0,2]
    "XXXYZ", // [3,1,1]
    "XXXYY", // [3,2,0]
    "XXXXZ", // [4,0,1]
    "XXXXY", // [4,1,0]
    "XXXXX", // [5,0,0]
];

/// 笛卡尔 h 壳层的角动量分量 (l, m, n)
pub const CAR_H_LMNS: [[usize; 3]; 21] = [
    [0, 0, 5],
    [0, 1, 4],
    [0, 2, 3],
    [0, 3, 2],
    [0, 4, 1],
    [0, 5, 0],
    [1, 0, 4],
    [1, 1, 3],
    [1, 2, 2],
    [1, 3, 1],
    [1, 4, 0],
    [2, 0, 3],
    [2, 1, 2],
    [2, 2, 1],
    [2, 3, 0],
    [3, 0, 2],
    [3, 1, 1],
    [3, 2, 0],
    [4, 0, 1],
    [4, 1, 0],
    [5, 0, 0],
];

/// 球谐 d 壳层的轨道符号（D 0、D+1、D-1、D+2、D-2）
pub const SPH_D_SYMS: [&str; 5] = [
    "D 0", // z^2
    "D+1", // xz
    "D-1", // yz
    "D+2", // x^2-y^2
    "D-2", // xy
];
/// 球谐 f 壳层的轨道符号（F 0、F±1、F±2、F±3）
#[rustfmt::skip]
pub const SPH_F_SYMS:[&str; 7]=[
    "F 0", // z^3
    "F+1", // xz^2
    "F-1", // yz^2
    "F+2", // z(x^2-y^2)
    "F-2", // xyz
    "F+3", // x(x^2-3y^2)
    "F-3", // y(3x^2-y^2)
];
/// 球谐 g 壳层的轨道符号（G 0、G±1、G±2、G±3、G±4）
#[rustfmt::skip]
pub const SPH_G_SYMS:[&str; 9]=["G 0","G+1","G-1","G+2","G-2","G+3","G-3","G+4","G-4"];
/// 球谐 h 壳层的轨道符号（H 0、H±1…H±5）
#[rustfmt::skip]
pub const SPH_H_SYMS:[&str;11]=["H 0","H+1","H-1","H+2","H-2","H+3","H-3","H+4","H-4","H+5","H-5"];

/// 所有壳层角动量分量按 s、p、d、f、g、h 顺序拼接而成的查找表
pub static LMNS: LazyLock<Vec<[usize; 3]>> = LazyLock::new(|| [&MIX_S_LMNS[..], &MIX_P_LMNS[..], &CAR_D_LMNS[..], &CAR_F_LMNS[..], &CAR_G_LMNS[..], &CAR_H_LMNS[..]].concat());

/// 所有壳层轨道符号按 s、p、d、f、g、h 顺序拼接而成的查找表
pub static SYMS: LazyLock<Vec<&str>> = LazyLock::new(|| [&MIX_S_SYMS[..], &MIX_P_SYMS[..], &CAR_D_SYMS[..], &CAR_F_SYMS[..], &CAR_G_SYMS[..], &CAR_H_SYMS[..]].concat());

/// 获取混合类型（s、p）壳层的轨道符号列表
#[cached]
pub fn get_mix_syms() -> Vec<&'static str> {
    [&MIX_S_SYMS[..], &MIX_P_SYMS[..]].concat()
}

/// 获取笛卡尔类型（d、f、g、h）壳层的轨道符号列表
#[cached]
pub fn get_car_syms() -> Vec<&'static str> {
    [&CAR_D_SYMS[..], &CAR_F_SYMS[..], &CAR_G_SYMS[..], &CAR_H_SYMS[..]].concat()
}

/// 获取球谐类型（d、f、g、h）壳层的轨道符号列表
#[cached]
pub fn get_sph_syms() -> Vec<&'static str> {
    [&SPH_D_SYMS[..], &SPH_F_SYMS[..], &SPH_G_SYMS[..], &SPH_H_SYMS[..]].concat()
}

#[rustfmt::skip]
/// 根据角动量获取对应的轨道角动量分量（l, m, n）列表
/// 仅支持0-3（s, p, d, f）
pub fn ang2lmns(ang: usize) -> Vec<[usize; 3]> {
    match ang {
        0 => MIX_S_LMNS.clone().to_vec(),
        1 => MIX_P_LMNS.clone().to_vec(),
        2 => CAR_D_LMNS.clone().to_vec(),
        3 => CAR_F_LMNS.clone().to_vec(),
        _ => panic!("未知角动量"),
    }
}
/// 根据角动量和形式获取轨道符号列表
/// - form: "car" 或 "sph"
pub fn ang2syms(ang: usize, form: &str) -> Vec<&str> {
    match (form, ang) {
        ("car", 0) => MIX_S_SYMS.clone().to_vec(),
        ("car", 1) => MIX_P_SYMS.clone().to_vec(),
        ("car", 2) => CAR_D_SYMS.clone().to_vec(),
        ("car", 3) => CAR_F_SYMS.clone().to_vec(),
        ("car", 4) => CAR_G_SYMS.clone().to_vec(),
        ("car", 5) => CAR_H_SYMS.clone().to_vec(),
        ("sph", 0) => MIX_S_SYMS.clone().to_vec(),
        ("sph", 1) => MIX_P_SYMS.clone().to_vec(),
        ("sph", 2) => SPH_D_SYMS.clone().to_vec(),
        ("sph", 3) => SPH_F_SYMS.clone().to_vec(),
        ("sph", 4) => SPH_G_SYMS.clone().to_vec(),
        ("sph", 5) => SPH_H_SYMS.clone().to_vec(),
        _ => panic!("未知角动量:{},{}", form, ang),
    }
}

/// 将角动量分量(l,m,n)转换为对应的轨道符号
pub fn lmn2sym(lmn: &[usize; 3]) -> String {
    let syms = &SYMS;
    let lmns = &LMNS;
    for i in 0..syms.len() {
        if &lmns[i] == lmn {
            return syms[i].to_string();
        }
    }
    panic!("未知角动量分量");
}

/// 将轨道符号转换为对应的角动量分量(l,m,n)
#[cached(key = "String", convert = r"{ sym.to_string() }")]
pub fn sym2lmn(sym: &str) -> [usize; 3] {
    let syms = &SYMS;
    let lmns = &LMNS;
    for i in 0..syms.len() {
        if syms[i] == sym {
            return lmns[i];
        }
    }
    panic!("未知角动量分量");
}

/// 根据壳层长度获取轨道符号列表
#[cached]
pub fn len2syms(len: usize) -> Vec<String> {
    let syms = match len {
        1 => MIX_S_SYMS.clone().to_vec(),
        3 => MIX_P_SYMS.clone().to_vec(),

        5 => SPH_D_SYMS.clone().to_vec(),
        7 => SPH_F_SYMS.clone().to_vec(),
        9 => SPH_G_SYMS.clone().to_vec(),
        11 => SPH_H_SYMS.clone().to_vec(),

        6 => CAR_D_SYMS.clone().to_vec(),
        10 => CAR_F_SYMS.clone().to_vec(),
        15 => CAR_G_SYMS.clone().to_vec(),
        21 => CAR_H_SYMS.clone().to_vec(),

        _ => panic!("未知长度"),
    };
    syms.iter().map(|sym| sym.to_string()).collect()
}

/// 将轨道符号转换为对应的角动量量子数
#[cached(key = "String", convert = r"{ sym.to_string() }")]
pub fn sym2ang(sym: &str) -> usize {
    for s_sym in MIX_S_SYMS {
        if sym == s_sym {
            return 0;
        }
    }
    for p_sym in MIX_P_SYMS {
        if sym == p_sym {
            return 1;
        }
    }
    for d_sym in CAR_D_SYMS {
        if sym == d_sym {
            return 2;
        }
    }
    for d_sym in SPH_D_SYMS {
        if sym == d_sym {
            return 2;
        }
    }
    for f_sym in SPH_F_SYMS {
        if sym == f_sym {
            return 3;
        }
    }
    for f_sym in CAR_F_SYMS {
        if sym == f_sym {
            return 3;
        }
    }
    for g_sym in CAR_G_SYMS {
        if sym == g_sym {
            return 4;
        }
    }
    for g_sym in SPH_G_SYMS {
        if sym == g_sym {
            return 4;
        }
    }
    panic!("未知符号:{}", sym);
}

/// 根据符号获取壳层类型（0:混合类型, 1:球谐类型, 2:笛卡尔类型）
pub fn sym2stp(sym: &str) -> usize {
    for mix_sym in get_mix_syms() {
        if sym == mix_sym {
            return 0;
        }
    }
    for sph_sym in get_sph_syms() {
        if sym == sph_sym {
            return 1;
        }
    }
    for car_sym in get_car_syms() {
        if sym == car_sym {
            return 2;
        }
    }
    panic!("未知符号:{}", sym);
}
