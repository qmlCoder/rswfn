use rswfn::base::Mole;
use rswfn::props::gridprop::{density, dftgrid};
use rswfn::maths::calc_dist;

#[test]
/// 解析电子静电势与实空间格点积分 -∫ρ(r')/|r-r'|d³r' 一致
///
/// 解析实现见 `rsint::esp`（参考 Multiwfn 的 `eleesp1`），
/// 这里用独立的数值积分交叉验证，同时验证量纲与符号。
fn test_esp_ele_vs_grid() {
    let mole = Mole::from_file("moles/C2H4.out");
    let caler = density::Calculator::new(&mole);

    // 较密的分子积分格点，保证数值积分的基本精度
    let (grids, weits) = dftgrid::Calculator::new(&mole).mol_grids(75, 302);
    let dens: Vec<f64> = grids.iter().map(|g| caler.mol_rho_dm(g, 0).val0).collect();

    for rc in [[1.2, 0.7, -0.4], [-2.5, 1.8, 0.9], [0.0, 0.0, 3.5]] {
        let ana = caler.ele_esp(&rc);
        let mut num = 0.0;
        for ((g, w), d) in grids.iter().zip(weits.iter()).zip(dens.iter()) {
            let dist = calc_dist(&rc, g);
            // 查询点与积分格点重合时 1/|r-r'| 有奇点，直接跳过（其真实贡献很小）
            if dist < 1e-3 {
                continue;
            }
            num -= w * d / dist;
        }
        println!("rc={rc:?} 解析={ana:.6}  格点={num:.6}  差={:.2e}", (ana - num).abs());
        assert!(ana < 0.0, "电子静电势应为负值: {ana}");
        assert!((ana - num).abs() < 2e-2, "rc {rc:?}: 解析 {ana} vs 格点 {num}");
    }
}

#[test]
/// 中性分子远端的静电势应趋于零（不会出现虚假的单极势）
fn test_esp_neutral_far_field() {
    let mole = Mole::from_file("moles/C2H4.out");
    let caler = density::Calculator::new(&mole);
    // C2H4 无偶极矩，远端以四极矩 1/r³ 衰减
    for r in [20.0, 40.0] {
        let esp = caler.tot_esp(&[r, 0.0, 0.0]);
        println!("r={r:>4} esp={esp:.3e}");
        assert!(esp.abs() < 1e-3, "r={r}: esp={esp}");
    }
}

#[test]
/// 原子核及其附近的核势按 Multiwfn 约定取 1e3，而不是发散
fn test_nuc_esp_on_nucleus() {
    let mole = Mole::from_file("moles/C2H4.out");
    let caler = density::Calculator::new(&mole);
    let xyz = mole.atoms().xyzs()[0];
    assert_eq!(caler.nuc_esp(&xyz), 1e3);

    // 立方格点常常与核只差浮点误差，这时也应返回有限值（否则 cub 出现 1e15 的尖峰）
    let near = [xyz[0] + 1e-9, xyz[1] - 1e-9, xyz[2]];
    assert_eq!(caler.nuc_esp(&near), 1e3);

    // 离核稍远时按 Z/r 计算
    let far = [xyz[0] + 0.5, xyz[1], xyz[2]];
    let val = caler.nuc_esp(&far);
    assert!(val > 0.0 && val < 1e3, "far = {val}");
}
