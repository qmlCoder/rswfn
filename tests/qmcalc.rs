use rswfn::base::{Atoms, Basis, Mole, Shl};
use rswfn::datas::btrans;
use rswfn::props::gridprop::wfnfunc;
use rswfn::qmcal::hf::RHF;
use rswfn::qmcal::ks::{RKS, UKS};
use rswfn::rsxcs::XcFunc;
use rswfn::writer::{FchWriter, MfnWriter, Writer};

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

    let res = RHF::new(&mole).run();
    assert!(res.converged, "SCF 未收敛");
    assert!((res.total_eng() + 1.1167).abs() < 1.0e-3, "H2/STO-3G 能量异常: {}", res.total_eng());
}

/// 乙烯 RHF/6-31G(d)：检查自洽收敛、电子数和总能量
#[test]
fn test_rhf_c2h4() {
    let mole = Mole::from_file("moles/C2H4.out");
    assert_eq!(mole.atoms.neles, [8, 8]);

    let res = RHF::new(&mole).run();
    assert!(res.converged, "SCF 未收敛（迭代 {} 次）", res.niter);
    // DIIS 正常工作时乙烯应 10 轮以内收敛（无 DIIS 时约 19 轮）
    assert!(res.niter <= 12, "SCF 迭代次数过多: {}（DIIS 可能未生效）", res.niter);

    // 乙烯 RHF/6-31G(d) 总能量约 -78.03 Hartree
    // （同一日志里的 B3LYP 能量是 -78.587，HF 因缺少相关能而高出约 0.55 Hartree）
    let total = res.total_eng();
    assert!((total + 78.031).abs() < 0.05, "总能量异常: {} Hartree（电子 {}, 核排斥 {}）", total, res.ele_eng, res.nuc_eng);

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
    assert!((total + 78.5875).abs() < 0.01, "B3LYP 总能量异常: {total} Hartree（Gaussian 参考 -78.5875）");

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
    assert!((uks.total_eng() - rks.total_eng()).abs() < 1.0e-6, "UKS({}) 与 RKS({}) 总能量不一致", uks.total_eng(), rks.total_eng());
    assert!(uks.spin_sq.abs() < 1.0e-6, "闭壳层的 <S^2> 应为 0，实际 {}", uks.spin_sq);

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
    assert!((res.spin_sq - 0.75).abs() < 0.02, "开壳层 <S^2> 偏离 0.75 过多: {}", res.spin_sq);
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
    let res = RHF::new(&mole).run();
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

/// 局部坐标系 π-only KS（LKS）：乙烯 B3LYP/6-31G(d)
///
/// 验证 π 约化基远小于全基组、π 电子数正确、SCF 收敛、局部系数表标签只含 π 成分，
/// 且系数可展开回全基组并导出 fch / mwfn 后回读一致。
#[test]
fn test_lks_c2h4() {
    use rswfn::qmcal::lks::LKS;

    let mole = Mole::from_file("moles/C2H4.out");
    let res = LKS::with_grid(&mole, 30, 50, XcFunc::B3LYP, None).run();

    assert!(res.converged, "LKS 未收敛（迭代 {} 次）", res.niter);
    assert_eq!(res.npi_elec, 2, "乙烯 π 电子数应为 2");
    // 6-31G(d)：每个 C 有 2 条 p 壳层（各 1 个 π 分量）与 1 条 d 壳层（2 个 π 分量），
    // 因此 π 约化基 = 2 × (2 + 2) = 8
    assert_eq!(res.npi, 8, "乙烯 π 约化基大小异常: {}", res.npi);
    assert!(res.npi < mole.basis().shls.nato("car"), "约化基应小于全基组");
    assert_eq!(res.nocc(), 1, "乙烯 π 占据轨道数应为 1");
    assert_eq!(res.obt_engs.len(), res.npi);
    assert_eq!(res.cmat.shape(), [res.npi, res.npi]);

    // 系数可展开回完整笛卡尔基
    let full = res.expand_coefs(&mole);
    assert_eq!(full.shape(), [mole.basis().shls.nato("car"), res.npi]);
    let nz = (0..full.nrows()).filter(|&i| full[(i, 0)].abs() > 1.0e-10).count();
    assert!(nz > 0, "占据 π 轨道展开到全局基后为空");

    // 局部 π 系数表：行标签只含 π 成分
    let table = res.mo_table_local(&mole);
    assert!(table.contains(":PX"), "局部系数表缺少 π 标签");
    assert_eq!(res.local_labels(&mole).len(), res.npi, "π 行标签数量应与 π 轨道数一致");

    // 导出 fch 并回读
    let out_mole = res.to_mole(&mole);
    let fch_path = std::env::temp_dir().join("rswfn_lks_c2h4.fch");
    let fch_str = fch_path.to_str().expect("路径转字符串失败");
    let mut fw = FchWriter::from_mole(&out_mole);
    fw.energy = res.ele_eng;
    fw.charge = out_mole.charge();
    fw.multip = out_mole.multip();
    fw.save(fch_str);
    let back = Mole::from_file(fch_str);
    assert_eq!(back.basis().shls.nato("car"), mole.basis().shls.nato("car"), "fch 基函数数不一致");
    assert_eq!(back.coefs().get_cmat("car", &back).ncols(), res.npi, "fch π 轨道数不一致");
    assert_eq!(back.coefs().obt_eles(&back), [res.nocc(), res.nocc()], "fch π 电子数不一致");
    std::fs::remove_file(&fch_path).ok();

    // 导出 mwfn 并回读
    let mwfn_path = std::env::temp_dir().join("rswfn_lks_c2h4.mwfn");
    let mwfn_str = mwfn_path.to_str().expect("路径转字符串失败");
    let mut mw = MfnWriter::from_mole(&out_mole);
    mw.energy = res.ele_eng;
    mw.save(mwfn_str);
    let back = Mole::from_file(mwfn_str);
    assert_eq!(back.basis().shls.nato("car"), mole.basis().shls.nato("car"), "mwfn 基函数数不一致");
    assert_eq!(back.coefs().get_cmat("car", &back).ncols(), res.npi, "mwfn π 轨道数不一致");
    assert_eq!(back.coefs().obt_eles(&back), [res.nocc(), res.nocc()], "mwfn π 电子数不一致");
    std::fs::remove_file(&mwfn_path).ok();
}

/// 球谐（5D/7F）基组下 KS 格点求值：
/// 球谐壳层没有单独的笛卡尔 (l, m, n)，格点上的 χ 必须先算该壳层的笛卡尔分量，
/// 再按 `φ_sph = car2sph_tmat · φ_car` 组合。
#[test]
fn test_spherical_basis_ks_grid() {
    let mole = Mole::from_file("moles/C2H4.out");
    assert_eq!(mole.basis().int_form(), "car", "fixture 应为 6D 基组");

    // 把 6D 基组统一成球谐形式，构造一个不含系数的球谐分子
    let sph_basis = mole.basis().unify_form("sph");
    let mole_sph = Mole::new(mole.atoms.clone(), Some(sph_basis), None);
    assert_eq!(mole_sph.basis().int_form(), "sph");

    // 修复前这里会在 sym2lmn("D 0") 处 panic
    let rks = RKS::with_grid(&mole_sph, 20, 50, XcFunc::B3LYP);
    let nato = rks.nato;
    assert_eq!(nato, mole_sph.basis().shls.nato("sph"), "球谐基函数数量不一致");
    assert_eq!(rks.chi0.len(), rks.grids.len() * nato);

    // 取一个 d 壳层，校验格点值等于笛卡尔分量的 car2sph 组合
    let basis = mole_sph.basis();
    let atos = basis.atos();
    let shls = &basis.shls;
    let uplo_sph = shls.uplo("sph");
    let uplo_car = shls.uplo("car");
    let s = (0..shls.len()).find(|&s| shls.get(s).ang == 2).expect("基线应包含 d 壳层");
    let ang = shls.get(s).ang;
    let [su, sl] = uplo_sph[s];
    let [cu, cl] = uplo_car[s];
    // 球谐基函数应为归一化的：S_sph 的 d 块对角应接近 1，
    // 否则用日志里的球谐 MO 系数构造的密度会失真
    let smat_sph = mole_sph.basis().get_smat("sph", &mole_sph);
    for i in su..sl {
        assert!((smat_sph[(i, i)] - 1.0).abs() < 1.0e-6, "sph 基函数未归一化: S[{i},{i}] = {}", smat_sph[(i, i)]);
    }

    let tmat = btrans::car2sph_tmat(ang);
    let grid = rks.grids[0];
    let atoms = mole_sph.atoms();
    let ncar = cl - cu;
    let mut car0 = vec![0.0; ncar];
    let mut car1 = vec![[0.0; 3]; ncar];
    for (j, idx) in (cu..cl).enumerate() {
        let ato = &atos.cars[idx];
        let [ax, ay, az] = atoms.get(ato.atm).xyz;
        let xyz = [grid[0] - ax, grid[1] - ay, grid[2] - az];
        let (val0, val1, _) = wfnfunc::cgf(&xyz, ato.lmn(), &ato.alps, ato.coes_norm(), 1);
        car0[j] = val0;
        car1[j] = val1;
    }
    for k in 0..(sl - su) {
        let mut expected0 = 0.0;
        let mut expected1 = [0.0; 3];
        for j in 0..ncar {
            let c = tmat[(k, j)];
            expected0 += c * car0[j];
            for r in 0..3 {
                expected1[r] += c * car1[j][r];
            }
        }
        assert!((rks.chi0[su + k] - expected0).abs() < 1.0e-10, "sph χ[{k}] 与笛卡尔组合不一致");
        for r in 0..3 {
            assert!((rks.chi1[r * nato + su + k] - expected1[r]).abs() < 1.0e-10, "sph ∂χ/∂{r}[{k}] 与笛卡尔组合不一致");
        }
    }
}
