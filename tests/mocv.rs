use ndarray::Array2;
use rswfn::props::orbtprop::mocv;

/// 断言矩阵为单位矩阵（允许浮点误差）
fn assert_identity(tmat: &Array2<f64>, tol: f64) {
    let n = tmat.nrows();
    assert_eq!(tmat.ncols(), n, "变换矩阵必须是方阵");
    for i in 0..n {
        for j in 0..n {
            let expect = if i == j { 1.0 } else { 0.0 };
            assert!(
                (tmat[(i, j)] - expect).abs() <= tol,
                "tmat[{},{}] = {}，期望为 {}",
                i,
                j,
                tmat[(i, j)],
                expect
            );
        }
    }
}

#[test]
/// 无旋转（单位旋转矩阵）时，笛卡尔变换矩阵应为单位矩阵
fn test_tmat_car_identity() {
    let lmat = Array2::<f64>::eye(3);
    for num in [3, 6, 10, 15] {
        let tmat = mocv::get_tmat_car(num, &lmat);
        assert_identity(&tmat, 1e-12);
    }
}

#[test]
/// 无旋转（单位旋转矩阵）时，球谐变换矩阵应为单位矩阵
fn test_tmat_sph_identity() {
    let lmat = Array2::<f64>::eye(3);
    for num in [3, 5, 7, 9] {
        let tmat = mocv::get_tmat_sph(num, &lmat);
        assert_identity(&tmat, 1e-12);
    }
}
