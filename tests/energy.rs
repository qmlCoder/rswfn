use rswfn::base::Mole;
use rswfn::moleprop::energy;

/// 用苯（C6H6，B3LYP/6-31G(d)）的 log 文件做 π 能量分解（POCV / MOCV）集成测试
///
/// 验证：
/// 1. `E_π + E_σ + E_πσ = E_elec`（分解闭合）；
/// 2. `E_total = E_elec + E_nuc` 与 log 中的 `SCF Done` 总能量一致；
/// 3. π 贡献 `E_π + E_πσ` 有限且为负。
#[test]
fn test_pi_energy_c6h6() {
    // log 中的 E(RB3LYP)（含核排斥），单位 Hartree
    const E_TOTAL_LOG: f64 = -232.247385084;

    let mole = Mole::from_file("moles/C6H6.out");
    let caler = energy::Calculator::new(&mole);

    let pocv = caler.pi_pocv_parts().expect("POCV π 能量分解失败");
    let mocv = caler.pi_mocv_parts().expect("MOCV π 能量分解失败");

    for (name, p) in [("POCV", &pocv), ("MOCV", &mocv)] {
        println!(
            "{name}: E_pi={:.6}  E_sigma={:.6}  E_cross={:.6} (库仑={:.6} 交换={:.6} XC={:.6})  E_elec={:.6}  E_nuc={:.6}  E_total={:.6}  pi={:.6}",
            p.e_pi,
            p.e_sigma,
            p.e_cross,
            p.e_cross_coul,
            p.e_cross_exch,
            p.e_cross_xc,
            p.e_elec,
            p.e_nuc,
            p.e_total,
            p.pi()
        );

        // 1. 分解闭合：E_π + E_σ + E_πσ == E_elec
        let sum = p.e_pi + p.e_sigma + p.e_cross;
        assert!(
            (sum - p.e_elec).abs() < 1.0e-6,
            "{name}: 分解不闭合 {sum} vs E_elec={}",
            p.e_elec
        );

        // 2. 总能量与 log 的 SCF Done 对照（用同一密度，网格/泛函微差，放宽到 0.02 Ha）
        assert!(
            (p.e_total - E_TOTAL_LOG).abs() < 0.02,
            "{name}: E_total={} 与 log {E_TOTAL_LOG} 不符",
            p.e_total
        );

        // 3. π 贡献有限且为负（电子能量份额）
        assert!(p.pi().is_finite(), "{name}: π 贡献非有限值");
        assert!(p.pi() < 0.0, "{name}: π 贡献应为负，实际 {}", p.pi());
    }
}
