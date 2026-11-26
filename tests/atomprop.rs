use rswfn::props::atomprop;
use rswfn::base::Mole;

mod utils;
use utils::allclose;

#[test]
/// 测试苯环mulliken电荷
fn test1() {
    let path = "moles/C6H6.fch";
    let mole = Mole::from_file(path);

    let caler = atomprop::charge::Calculator::new(&mole);
    let res = caler.mulliken();
    let ans = vec![
        6.12846023, 6.12846024, 6.12846023, 6.12846023, 6.12846024, 6.12846023, 0.87153977, 0.87153977, 0.87153977,
        0.87153977, 0.87153977, 0.87153977,
    ];
    assert!(allclose(res, ans, 1e-4));
}

#[test]
/// 测试苯环lowdin电荷
fn test2() {
    let path = "moles/C6H6.fch";
    let mole = Mole::from_file(path);
    let caler = atomprop::charge::Calculator::new(&mole);
    let res = caler.lowdin();
    let ans = vec![
        6.15647202, 6.15982449, 6.15982449, 6.15647202, 6.15982449, 6.15982449, 0.84289460, 0.84049220, 0.84049220,
        0.84289460, 0.84049220, 0.84049220,
    ];
    assert!(allclose(res, ans, 1e-4));
}

#[test]
/// 测试苯环hirshfeld电荷
fn test3() {
    let path = "moles/C6H6.fch";
    let mole = Mole::from_file(path);
    let caler = atomprop::charge::Calculator::new(&mole);
    let mut res = caler.hirshfeld();
    let res = (0..res.len()).map(|atm| mole.atoms.get(atm).crg - res[atm]).collect();
    let ans = vec![
        -0.04268004,
        -0.04261213,
        -0.04261213,
        -0.04268004,
        -0.04261213,
        -0.04261213,
        0.04262433,
        0.04263998,
        0.04263999,
        0.04262433,
        0.04263998,
        0.04263999,
    ];
    assert!(allclose(res, ans, 1e-4));
}

#[test]
/// 测试苯环pocv方法pi电子布居
fn test4() {
    let path = "moles/C6H6.fch";
    let mole = Mole::from_file(path);

    let caler = atomprop::charge::Calculator::new(&mole);
    let (dirs, vals) = caler.pi_pocv("mulliken").unwrap();
    let ans = vec![
        0.98514, 0.98514, 0.98514, 0.98514, 0.98514, 0.98514, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    ];
    println!("pi_pocv:{:?}", vals);
    assert!(allclose(vals, ans, 1e-4));
}

#[test]
/// 测试苯环mocv方法pi电子布居
fn test5() {
    let path = "moles/C6H6.fch";
    let mole = Mole::from_file(path);

    let caler = atomprop::charge::Calculator::new(&mole);
    let (dirs, atos, vals) = caler.pi_mocv("mulliken").unwrap();
    let ans = vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    println!("pi_mocv:{:?}", vals);
    assert!(allclose(vals, ans, 1e-4));
}

#[test]
/// 测试C6H6^- 自旋布居
fn test6() {
    let path = "moles/C6H6_N.fch";
    let mole = Mole::from_file(path);
    let caler = atomprop::charge::Calculator::new(&mole);
    let res = caler.spin("mulliken").iter().map(|row| row[2]).collect();
    let ans = vec![
        -0.11890, 0.32872, 0.32872, -0.11890, 0.32872, 0.32872, 0.00392, -0.02123, -0.02123, 0.00392, -0.02123,
        -0.02123,
    ];
    assert!(allclose(res, ans, 1e-4));
}

#[test]
/// 测试苯酚福井函数
fn test7() {
    let mole0 = Mole::from_file("moles/phenol.out");
    let molen = Mole::from_file("moles/phenol_N+1.out");
    let molep = Mole::from_file("moles/phenol_N-1.out");
    let caler = atomprop::activity::Calculator::new(&mole0);
    // 比较双描述符
    let res: Vec<f64> = caler
        .fukui(&molen, &molep, "hirshfeld")
        .iter()
        .map(|row| row[6])
        .collect();
    let ans = vec![
        -0.0935, 0.0521, 0.0429, -0.0741, 0.0429, 0.0521, -0.0129, 0.0164, 0.0151, 0.0151, 0.0164, -0.0522, -0.0202,
    ];
    println!("res:{:?}", res);
    assert!(allclose(res, ans, 1e-3));
}

#[test]
/// 测试三种轨道形式（raw/car/sph）计算的总电子数一致
fn test8() {
    let mole = Mole::from_file("moles/C2H4.out");
    // C2H4 为 6D 笛卡尔基组，car/sph 转换非平凡，适合验证总电子数守恒
    let sums: Vec<f64> = ["raw", "car", "sph"]
        .into_iter()
        .map(|form| {
            let dmat = mole.coefs().get_dmat(form, &mole);
            let smat = mole.basis().get_smat(form, &mole);
            (dmat * smat).sum()
        })
        .collect();
    println!("sums:{:?}", sums);
    // C2H4 共 8α+8β = 16 个电子
    assert!(allclose(sums.clone(), vec![16.0, 16.0, 16.0], 1e-3));
    assert!((sums[0] - sums[2]).abs() < 1e-3);
}

#[test]
/// 测试乙烯的标准两阶段 RESP 电荷
fn test10() {
    let mole = Mole::from_file("moles/C2H4.out");
    let caler = atomprop::resp::Calculator::new(&mole);
    let res = caler.resp().unwrap();
    println!("resp C2H4: {:?}", res.charges);
    println!("MK 点数: {}, RMSE: {:.3e}", res.n_points, res.rmse);

    assert!(res.n_points > 0);
    assert!(res.stage2 && res.converged);
    assert!(res.rmse.is_finite() && res.rmse > 0.0);

    // 总电荷守恒
    let total: f64 = res.charges.iter().sum();
    assert!((total - mole.charge() as f64).abs() < 1e-6);

    // 乙烯的对称等价性：两个 C 等价，四个 H 等价
    assert!((res.charges[0] - res.charges[3]).abs() < 1e-4);
    let h = res.charges[1];
    for i in [2, 4, 5] {
        assert!((res.charges[i] - h).abs() < 1e-4);
    }
}

#[test]
/// 测试苯 RESP 电荷与 Multiwfn 3.8 的一致性
fn test11() {
    let mole = Mole::from_file("moles/C6H6.fch");
    let caler = atomprop::resp::Calculator::new(&mole);
    let res = caler.resp().unwrap();
    let ans = vec![
        -0.1079158845,
        -0.1116391792,
        -0.1116391757,
        -0.1079158847,
        -0.1133629517,
        -0.1133629486,
        0.1100701395,
        0.1109257558,
        0.1109257566,
        0.1100701396,
        0.1119221161,
        0.1119221169,
    ];
    println!("RESP C6H6: {:?}", res.charges);
    // 苯没有 sp³ 碳，第二阶段应被跳过
    assert!(!res.stage2);
    // 与 Multiwfn 的 MK 采样存在极小的浮点差异，容差取 5e-4
    assert!(allclose(res.charges, ans, 5e-4));
}

#[test]
/// 测试 MOCV 方法在笛卡尔型（6D）输入基组下的 π 电子布居
fn test9() {
    let mole = Mole::from_file("moles/C2H4.out");
    // C2H4.out 的 raw 系数为笛卡尔形式，验证 MOCV 在非球谐输入下仍能正确工作
    let caler = atomprop::charge::Calculator::new(&mole);
    let (_stms, _atos, vals) = caler.pi_mocv("mulliken").unwrap();
    let ans = vec![1.0, 0.0, 0.0, 1.0, 0.0, 0.0];
    println!("pi_mocv(car):{:?}", vals);
    assert!(allclose(vals, ans, 1e-4));
}
