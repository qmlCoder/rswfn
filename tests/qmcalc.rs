use rswfn::base::{Atoms, Basis, Mole, Shl};
use rswfn::qmcal::hf::RHF;
use rswfn::qmcal::ks::{RKS, UKS};
use rswfn::rsxcs::XcFunc;
use rswfn::writer::{FchWriter, Writer};

/// H2 / STO-3G（R = 1.4 bohr）的参考 HF 总能量约 -1.1167 Hartree
/// 这个用例很快，用已知解析参考值校验积分与 SCF 主流程
#[test]
fn test_rhf_h2_sto3g() {
    let alps = vec![3.42525091, 0.62391373, 0.16885540];
    let coes = vec![0.15432897, 0.53532814, 0.44463500];
    let shls = vec![
        Shl {
            atm: 0,
            shl: 0,
            ang: 0,
            form: "sph".to_string(),
            alps: alps.clone(),
            coes: coes.clone(),
        },
        Shl {
            atm: 1,
            shl: 0,
            ang: 0,
            form: "sph".to_string(),
            alps,
            coes,
        },
    ];
    let atoms = Atoms::new(&vec![1, 1], &vec![[0.0, 0.0, 0.0], [0.0, 0.0, 1.4]], &vec![1.0, 1.0]);
    let mole = Mole::new(atoms, Some(Basis::new(shls, None)), None);

    let res = RHF::new(&mole).run(None);
    assert!(res.converged, "SCF 未收敛");
    assert!(
        (res.total_eng() + 1.1167).abs() < 1.0e-3,
        "H2/STO-3G 能量异常: {}",
        res.total_eng()
    );
}

/// 乙烯 RHF/6-31G(d)：检查自洽收敛、电子数和总能量
#[test]
fn test_rhf_c2h4() {
    let mole = Mole::from_file("moles/C2H4.out");
    assert_eq!(mole.atoms.neles, [8, 8]);

    let res = RHF::new(&mole).run(None);
    assert!(res.converged, "SCF 未收敛（迭代 {} 次）", res.niter);
    // DIIS 正常工作时乙烯应 10 轮以内收敛（无 DIIS 时约 19 轮）
    assert!(res.niter <= 12, "SCF 迭代次数过多: {}（DIIS 可能未生效）", res.niter);

    // 乙烯 RHF/6-31G(d) 总能量约 -78.03 Hartree
    // （同一日志里的 B3LYP 能量是 -78.587，HF 因缺少相关能而高出约 0.55 Hartree）
    let total = res.total_eng();
    assert!(
        (total + 78.031).abs() < 0.05,
        "总能量异常: {} Hartree（电子 {}, 核排斥 {}）",
        total,
        res.ele_eng,
        res.nuc_eng
    );

    // 密度矩阵应给出 16 个电子
    let calced = Mole::new(mole.atoms.clone(), Some(res.basis), Some(res.coefs));
    assert_eq!(calced.coefs().obt_eles(&calced), [8, 8]);
}

/// 乙烯 B3LYP/6-31G(d)：同一日志里 Gaussian 的总能量为 -78.5875 Hartree。
/// 本实现与 Gaussian 的差异主要来自数值积分格点（约 1e-3 Hartree 量级）。
#[test]
fn test_rks_c2h4() {
    let mole = Mole::from_file("moles/C2H4.out");
    let res = RKS::with_grid(&mole, 50, 74, XcFunc::B3LYP).run();
    assert!(res.converged, "SCF 未收敛（迭代 {} 次）", res.niter);
    assert!(res.niter <= 12, "SCF 迭代次数过多: {}（DIIS 可能未生效）", res.niter);
    assert!(res.xc_eng < 0.0, "交换相关能应为负值: {}", res.xc_eng);

    let total = res.total_eng();
    println!("B3LYP/6-31G(d) 乙烯总能量 = {total}（Gaussian -78.5875）");
    assert!(
        (total + 78.5875).abs() < 0.01,
        "B3LYP 总能量异常: {total} Hartree（Gaussian 参考 -78.5875）"
    );

    let calced = Mole::new(mole.atoms.clone(), Some(res.basis), Some(res.coefs));
    assert_eq!(calced.coefs().obt_eles(&calced), [8, 8]);
}

/// 闭壳层极限：同一格点上 UKS 必须与 RKS 给出完全一致的能量
///
/// 当 ρ_α = ρ_β 时，自旋极化泛函退化为未极化泛函，
/// 且 UKS 的 Fock/能量表达式也与 RKS 代数等价，因此两者差异只来自浮点误差。
#[test]
fn test_uks_closed_shell_matches_rks() {
    let mole = Mole::from_file("moles/C2H4.out");
    let rks = RKS::with_grid(&mole, 50, 74, XcFunc::B3LYP).run();
    let uks = UKS::with_grid(&mole, 50, 74, XcFunc::B3LYP).run();

    assert!(uks.converged, "UKS 未收敛（迭代 {} 次）", uks.niter);
    println!("RKS = {}, UKS = {}", rks.total_eng(), uks.total_eng());
    assert!(
        (uks.total_eng() - rks.total_eng()).abs() < 1.0e-6,
        "UKS({}) 与 RKS({}) 总能量不一致",
        uks.total_eng(),
        rks.total_eng()
    );
    assert!(
        uks.spin_sq.abs() < 1.0e-6,
        "闭壳层的 <S^2> 应为 0，实际 {}",
        uks.spin_sq
    );

    let calced = Mole::new(mole.atoms.clone(), Some(uks.basis), Some(uks.coefs));
    assert_eq!(calced.coefs().obt_eles(&calced), [8, 8]);
}

/// 开壳层（C6H6^-）UKS：能收敛、电子数正确、<S^2> 接近 0.75
///
/// C6H6^- 的 ROB3LYP 参考总能量约 -232.1631 Hartree；
/// UKS 与 ROKS 略有差异，这里只做量级检查。
#[test]
fn test_uks_open_shell() {
    let mole = Mole::from_file("moles/C6H6_N.fch");
    assert_eq!(mole.atoms.neles, [22, 21]);

    let res = UKS::with_grid(&mole, 40, 50, XcFunc::B3LYP).run();
    let total = res.total_eng();
    assert!(res.converged, "UKS 未收敛（迭代 {} 次）", res.niter);
    println!("C6H6^- UKS 总能量 = {total}（ROB3LYP 参考 -232.1631）");
    assert!(
        (res.spin_sq - 0.75).abs() < 0.02,
        "开壳层 <S^2> 偏离 0.75 过多: {}",
        res.spin_sq
    );
    assert!((total + 232.1631).abs() < 0.05, "UKS 总能量异常: {total} Hartree");

    let calced = Mole::new(mole.atoms.clone(), Some(res.basis), Some(res.coefs));
    assert_eq!(calced.coefs().obt_eles(&calced), [22, 21]);

    // 开壳层 fch 应能完整回读（α/β 两套系数与轨道能）
    let path = std::env::temp_dir().join("rswfn_uks_c6h6n.fch");
    let path_str = path.to_str().expect("路径转字符串失败");
    let mut writer = FchWriter::new();
    writer.read_mole(&calced);
    writer.charge = calced.charge();
    writer.multip = calced.multip();
    writer.energy = total;
    writer.save(path_str);

    let back = Mole::from_file(path_str);
    assert_eq!(back.atoms.neles, [22, 21], "回读的电子数不一致");
    assert_eq!(back.coefs().obt_eles(&back), [22, 21], "回读的占据数不一致");
    std::fs::remove_file(&path).ok();
}

/// HF 结果写成 fch 后应能被重新读取，且电子数与系数一致
#[test]
fn test_rhf_fch_roundtrip() {
    let mole = Mole::from_file("moles/C2H4.out");
    let res = RHF::new(&mole).run(None);
    assert!(res.converged, "SCF 未收敛");

    let total_eng = res.total_eng();
    let calced = Mole::new(mole.atoms.clone(), Some(res.basis), Some(res.coefs));

    let path = std::env::temp_dir().join("rswfn_rhf_c2h4.fch");
    let path_str = path.to_str().expect("路径转字符串失败");
    let mut writer = FchWriter::new();
    writer.read_mole(&calced);
    writer.energy = total_eng;
    writer.save(path_str);

    let back = Mole::from_file(path_str);
    assert_eq!(back.atoms.neles, [8, 8], "回读的电子数不一致");
    assert_eq!(back.coefs().obt_eles(&back), [8, 8], "回读的占据数不一致");
    // 输入是 6-31G(d) 的 6D 形式，因此笛卡尔基函数数为 38
    assert_eq!(back.basis().int_form(), "car", "回读的基组形式不一致");
    assert_eq!(back.basis().shls.nato("car"), 38, "回读的基函数数量不一致");

    std::fs::remove_file(&path).ok();
}
