//! 量子化学计算任务（HF / DFT）。

use std::path::{Path, absolute};
use std::time::{Duration, Instant};

use crate::base::Mole;
use crate::qmcal;
use crate::report::{Field, Report};
use crate::rsbse;
use crate::rsxcs::{XcFunc, supported_names};
use crate::writer::{FchWriter, MfnWriter, Writer};

use super::{ServiceResult, ext_of};

/// `qmcal hf` / `qmcal dft` 的 SCF 数值参数
#[derive(Clone, Copy, Debug)]
pub struct ScfOptions {
    /// 双电子积分 Schwarz 筛选阈值
    pub screen: f64,
    /// 双电子积分块内存上限（MB），超过后自动退化为直接 SCF
    pub max_mem_mb: usize,
    /// SCF 能量收敛阈值（Hartree）
    pub conv: f64,
    /// 最大 SCF 迭代次数
    pub max_iter: usize,
}

impl Default for ScfOptions {
    fn default() -> Self {
        Self {
            screen: 1.0e-11,
            max_mem_mb: 8000,
            conv: 1.0e-8,
            max_iter: 100,
        }
    }
}

/// `qmcal dft` 的数值积分格点参数
#[derive(Clone, Copy, Debug)]
pub struct GridOptions {
    /// 径向格点数（默认 99，对应 Gaussian 的 UltraFine）
    pub radial: usize,
    /// 球面（Lebedev）格点数（默认 590，对应 Gaussian 的 UltraFine）
    pub angular: usize,
    /// 是否做径向-角向剪枝（默认开启，格点数约为不剪枝的 1/4）
    pub prune: bool,
    /// 格点积分中显著原子轨道的距离截断（Bohr；0 表示关闭）
    pub ao_cut: f64,
}

impl Default for GridOptions {
    fn default() -> Self {
        Self {
            radial: 99,
            angular: 590,
            prune: true,
            ao_cut: 12.0,
        }
    }
}

/// 读取 qmcal 的输入文件，并按需使用程序提供的基组数据覆盖文件中的基组
fn load_calc_mole(input: &Path, basis: Option<&str>) -> ServiceResult<(Mole, String)> {
    let from = absolute(input).map_err(|err| format!("无法解析输入路径 {}: {err}", input.display()))?;
    let from_path = from.to_str().ok_or_else(|| format!("路径不是有效的 UTF-8: {}", from.display()))?.to_string();
    let mut mole = super::load_any_mole(input)?;
    match basis {
        Some(name) => {
            let basis = rsbse::from_name(name, &mole.atoms).ok_or_else(|| format!("不支持的基组: {name}"))?;
            mole.basis = Some(basis);
            // 基组已更换，文件中原有的轨道系数不再对应
            mole.coefs = None;
        }
        None => {
            if mole.basis.is_none() {
                return Err("输入文件不含基组信息，请提供 .fch/.fchk、.log/.out 等包含基组的文件，或用 -b 指定基组".to_string());
            }
        }
    }
    Ok((mole, from_path))
}

/// 计算公共信息字段
fn calc_header_fields(path: &str, mole: &Mole, method: &str) -> Vec<Field> {
    let [nela, nelb] = mole.atoms.neles;
    let form = mole.basis().int_form();
    vec![
        Field::new("方法", method.to_string()),
        Field::new("输入文件", path.to_string()),
        Field::new("分子式", mole.formula()),
        Field::new("电子数", format!("{}（α {nela}，β {nelb}）", nela + nelb)),
        Field::new("基函数数", mole.basis().shls.nato(form).to_string()),
        Field::new("积分形式", form.to_string()),
    ]
}

/// SCF 收敛情况字段
fn calc_status_fields(converged: bool, niter: usize, delta_eng: f64) -> Vec<Field> {
    vec![
        Field::new("SCF 状态", if converged { "收敛" } else { "未收敛" }),
        Field::new("迭代次数", niter.to_string()),
        Field::new("ΔE", format!("{delta_eng:.2e}")),
    ]
}

/// 各阶段耗时（备注行）
fn timing_notes(t_load: Duration, t_grid: Option<Duration>, t_int1e: Duration, t_int2e: Duration, t_scf: Duration, t_export: Duration, t_total: Duration) -> Vec<String> {
    let mut notes = vec!["耗时统计".to_string()];
    notes.push(format!("  读取输入  : {:>9.3} s", t_load.as_secs_f64()));
    if let Some(t_grid) = t_grid {
        notes.push(format!("  格点初始化: {:>9.3} s", t_grid.as_secs_f64()));
    }
    notes.push(format!("  单电子积分: {:>9.3} s", t_int1e.as_secs_f64()));
    notes.push(format!("  双电子积分: {:>9.3} s", t_int2e.as_secs_f64()));
    notes.push(format!("  SCF 迭代  : {:>9.3} s", t_scf.as_secs_f64()));
    notes.push(format!("  导出波函数: {:>9.3} s", t_export.as_secs_f64()));
    notes.push(format!("  总耗时    : {:>9.3} s", t_total.as_secs_f64()));
    notes
}

/// 将计算得到的波函数导出为 .fch 文件，返回实际写入路径
fn write_fch(mole: &Mole, output: Option<&Path>, basis: crate::base::Basis, coefs: crate::base::Coefs, energy: f64) -> ServiceResult<Option<String>> {
    let Some(output) = output else {
        return Ok(None);
    };
    let dest = absolute(output).map_err(|err| format!("无法解析输出路径 {}: {err}", output.display()))?;
    if ext_of(&dest).as_deref() != Some("fch") {
        return Err(format!("计算结果目前仅支持导出 .fch 文件，收到: {}", dest.display()));
    }
    let dest_path = dest.to_str().ok_or_else(|| format!("输出路径不是有效的 UTF-8: {}", dest.display()))?;
    let calc_mole = Mole::new(mole.atoms.clone(), Some(basis), Some(coefs));
    let mut writer = FchWriter::new();
    writer.read_mole(&calc_mole);
    writer.charge = calc_mole.charge();
    writer.multip = calc_mole.multip();
    writer.energy = energy;
    writer.save(dest_path);
    Ok(Some(dest_path.to_string()))
}

/// 闭壳层 Hartree-Fock（RHF）计算
pub fn hf(input: &Path, basis: Option<&str>, output: Option<&Path>) -> ServiceResult<Report> {
    hf_with(input, basis, output, &ScfOptions::default())
}

/// 带 SCF 数值参数的闭壳层 Hartree-Fock（RHF）计算
pub fn hf_with(input: &Path, basis: Option<&str>, output: Option<&Path>, opts: &ScfOptions) -> ServiceResult<Report> {
    let t_total = Instant::now();
    let t_load = Instant::now();
    let (mole, from_path) = load_calc_mole(input, basis)?;
    let t_load = t_load.elapsed();

    let [nela, nelb] = mole.atoms.neles;
    if nela != nelb {
        return Err(format!("qmcal hf 目前仅支持闭壳层体系（α {nela} ≠ β {nelb}）；开壳层请使用 qmcal dft（UKS）"));
    }

    let mut report = Report::new("Hartree-Fock（RHF）计算");
    report.fields.extend(calc_header_fields(&from_path, &mole, "RHF"));

    let mut caler = qmcal::hf::RHF::new(&mole);
    caler.screen = opts.screen;
    caler.max_mem = opts.max_mem_mb * 1024 * 1024;
    caler.tol = opts.conv;
    caler.max_iter = opts.max_iter;
    let res = caler.run();
    let total_eng = res.total_eng();
    let (ele_eng, nuc_eng, converged, niter, delta_eng) = (res.ele_eng, res.nuc_eng, res.converged, res.niter, res.delta_eng);
    let (time_int1e, time_int2e, time_scf, time_jk) = (res.time_int1e, res.time_int2e, res.time_scf, res.time_jk);

    report.fields.extend(calc_status_fields(converged, niter, delta_eng));
    report.fields.push(Field::new("电子能量", format!("{ele_eng:>20.10} Hartree")));
    report.fields.push(Field::new("核排斥能", format!("{nuc_eng:>20.10} Hartree")));
    report.fields.push(Field::new("总能量", format!("{total_eng:>20.10} Hartree")));

    let t_export = Instant::now();
    let saved = write_fch(&mole, output, res.basis, res.coefs, total_eng)?;
    let t_export = t_export.elapsed();

    report.notes.extend(timing_notes(t_load, None, time_int1e, time_int2e, time_scf, t_export, t_total.elapsed()));
    report.notes.push(format!("  SCF 细分  : J/K 构造 {:>7.3} s", time_jk.as_secs_f64()));
    if let Some(path) = saved {
        report.notes.push(format!("波函数已导出到: {path}"));
    }
    if !converged {
        report.notes.push(format!("警告: SCF 在 {niter} 次迭代后仍未收敛，结果可能不可靠"));
    }
    Ok(report)
}

/// 密度泛函计算（闭壳层 RKS / 开壳层 UKS 自动选择）
pub fn dft(input: &Path, basis: Option<&str>, func: &str, output: Option<&Path>) -> ServiceResult<Report> {
    dft_with(input, basis, func, output, &ScfOptions::default(), &GridOptions::default())
}

/// 带 SCF 与格点参数的密度泛函计算
pub fn dft_with(input: &Path, basis: Option<&str>, func: &str, output: Option<&Path>, opts: &ScfOptions, grid: &GridOptions) -> ServiceResult<Report> {
    let t_total = Instant::now();
    let t_load = Instant::now();
    let (mole, from_path) = load_calc_mole(input, basis)?;
    let t_load = t_load.elapsed();

    let xc = XcFunc::from_name(func).ok_or_else(|| format!("未知的交换相关泛函: {func}（可用: {}）", supported_names().join(", ")))?;
    let func_name = xc.name().to_uppercase();
    let [nela, nelb] = mole.atoms.neles;

    if nela == nelb {
        let mut report = Report::new(format!("Kohn-Sham（RKS，{func_name}）计算"));
        report.fields.extend(calc_header_fields(&from_path, &mole, "RKS"));

        let t_grid = Instant::now();
        let mut caler = qmcal::ks::RKS::with_grid_prune(&mole, grid.radial, grid.angular, xc, grid.prune);
        caler.screen = opts.screen;
        caler.max_mem = opts.max_mem_mb * 1024 * 1024;
        caler.tol = opts.conv;
        caler.max_iter = opts.max_iter;
        caler.ao_cut = grid.ao_cut;
        let t_grid = t_grid.elapsed();
        report.fields.push(Field::new("格点", format!("{} 径向 × {} 球面（共 {} 个）", caler.radial, caler.angular, caler.grids.len())));
        if caler.ao_cut > 0.0 {
            let (nblk, mean_ao, max_ao) = caler.ao_stats();
            report.fields.push(Field::new(
                "显著 AO",
                format!("截断 {:.2} Bohr，{} 个格点块，平均 {:.0}/{:.0}（最多 {}）", caler.ao_cut, nblk, mean_ao, caler.nato as f64, max_ao),
            ));
        }

        let res = caler.run();
        let total_eng = res.total_eng();
        let (ele_eng, xc_eng, nuc_eng) = (res.ele_eng, res.xc_eng, res.nuc_eng);
        let (converged, niter, delta_eng) = (res.converged, res.niter, res.delta_eng);
        let (time_int1e, time_int2e, time_scf, time_jk, time_grid) = (res.time_int1e, res.time_int2e, res.time_scf, res.time_jk, res.time_grid);

        report.fields.extend(calc_status_fields(converged, niter, delta_eng));
        report.fields.push(Field::new("电子能量", format!("{ele_eng:>20.10} Hartree")));
        report.fields.push(Field::new("交换相关能", format!("{xc_eng:>19.10} Hartree")));
        report.fields.push(Field::new("核排斥能", format!("{nuc_eng:>20.10} Hartree")));
        report.fields.push(Field::new("总能量", format!("{total_eng:>20.10} Hartree")));
        report.notes.push(format!("  SCF 细分  : J/K 构造 {:>7.3} s，格点部分 {:>7.3} s", time_jk.as_secs_f64(), time_grid.as_secs_f64()));

        let t_export = Instant::now();
        let saved = write_fch(&mole, output, res.basis, res.coefs, total_eng)?;
        let t_export = t_export.elapsed();

        report.notes.extend(timing_notes(t_load, Some(t_grid), time_int1e, time_int2e, time_scf, t_export, t_total.elapsed()));
        if let Some(path) = saved {
            report.notes.push(format!("波函数已导出到: {path}"));
        }
        if !converged {
            report.notes.push(format!("警告: SCF 在 {niter} 次迭代后仍未收敛，结果可能不可靠"));
        }
        Ok(report)
    } else {
        let mut report = Report::new(format!("Kohn-Sham（UKS，{func_name}）计算"));
        report.fields.extend(calc_header_fields(&from_path, &mole, "UKS"));

        let t_grid = Instant::now();
        let mut caler = qmcal::ks::UKS::with_grid_prune(&mole, grid.radial, grid.angular, xc, grid.prune);
        caler.base.screen = opts.screen;
        caler.base.max_mem = opts.max_mem_mb * 1024 * 1024;
        caler.base.tol = opts.conv;
        caler.base.max_iter = opts.max_iter;
        caler.base.ao_cut = grid.ao_cut;
        let t_grid = t_grid.elapsed();
        report
            .fields
            .push(Field::new("格点", format!("{} 径向 × {} 球面（共 {} 个）", caler.base.radial, caler.base.angular, caler.base.grids.len())));
        if caler.base.ao_cut > 0.0 {
            let (nblk, mean_ao, max_ao) = caler.base.ao_stats();
            report.fields.push(Field::new(
                "显著 AO",
                format!("截断 {:.2} Bohr，{} 个格点块，平均 {:.0}/{:.0}（最多 {}）", caler.base.ao_cut, nblk, mean_ao, caler.base.nato as f64, max_ao),
            ));
        }

        let res = caler.run();
        let total_eng = res.total_eng();
        let (ele_eng, xc_eng, nuc_eng) = (res.ele_eng, res.xc_eng, res.nuc_eng);
        let (converged, niter, delta_eng) = (res.converged, res.niter, res.delta_eng);
        let (time_int1e, time_int2e, time_scf) = (res.time_int1e, res.time_int2e, res.time_scf);
        let (time_jk, time_grid) = (res.time_jk, res.time_grid);

        report.fields.extend(calc_status_fields(converged, niter, delta_eng));
        report.fields.push(Field::new("电子能量", format!("{ele_eng:>20.10} Hartree")));
        report.fields.push(Field::new("交换相关能", format!("{xc_eng:>19.10} Hartree")));
        report.fields.push(Field::new("核排斥能", format!("{nuc_eng:>20.10} Hartree")));
        report.fields.push(Field::new("总能量", format!("{total_eng:>20.10} Hartree")));
        report.fields.push(Field::new("<S^2>", format!("{:.6}（多重度 {}）", res.spin_sq, res.multip())));
        report.notes.push(format!("  SCF 细分  : J/K 构造 {:>7.3} s，格点部分 {:>7.3} s", time_jk.as_secs_f64(), time_grid.as_secs_f64()));

        let t_export = Instant::now();
        let saved = write_fch(&mole, output, res.basis, res.coefs, total_eng)?;
        let t_export = t_export.elapsed();

        report.notes.extend(timing_notes(t_load, Some(t_grid), time_int1e, time_int2e, time_scf, t_export, t_total.elapsed()));
        if let Some(path) = saved {
            report.notes.push(format!("波函数已导出到: {path}"));
        }
        if !converged {
            report.notes.push(format!("警告: SCF 在 {niter} 次迭代后仍未收敛，结果可能不可靠"));
        }
        Ok(report)
    }
}

/// 局部坐标系下的 π-only Kohn-Sham（DFT）计算（MOCV，π 约化基）
///
/// 只保留局部坐标系中 π 对称的原子轨道，在远小于全基组的约化基上做自洽场。
/// 因丢掉 σ 与核内层，其能量是 π 子系统模型值，**不是**分子总能量（参考
/// `docs/pi电子能量.md`）；主要用途是低成本获得 π 分子轨道与轨道能。
pub fn lks(input: &Path, basis: Option<&str>, func: &str, npi_elec: Option<usize>, output: Option<&Path>, print_mo: bool) -> ServiceResult<Report> {
    lks_with(input, basis, func, npi_elec, output, print_mo, &ScfOptions::default())
}

/// 带 SCF 数值参数的 π-only Kohn-Sham（LKS）计算
pub fn lks_with(input: &Path, basis: Option<&str>, func: &str, npi_elec: Option<usize>, output: Option<&Path>, print_mo: bool, opts: &ScfOptions) -> ServiceResult<Report> {
    let t_total = Instant::now();
    let t_load = Instant::now();
    let (mole, from_path) = load_calc_mole(input, basis)?;
    let t_load = t_load.elapsed();

    let xc = XcFunc::from_name(func).ok_or_else(|| format!("未知的交换相关泛函: {func}（可用: {}）", supported_names().join(", ")))?;
    let func_name = xc.name().to_uppercase();
    let [nela, nelb] = mole.atoms.neles;
    if nela != nelb {
        return Err(format!("qmcal lks 目前仅支持闭壳层体系（α {nela} ≠ β {nelb}）"));
    }

    let mut report = Report::new(format!("局部坐标系 π-only Kohn-Sham（LKS，{func_name}）计算"));
    report.fields.extend(calc_header_fields(&from_path, &mole, "LKS"));

    let t_grid = Instant::now();
    let mut caler = qmcal::lks::LKS::with_grid(&mole, 100, 74, xc, npi_elec);
    caler.screen = opts.screen;
    caler.max_mem = opts.max_mem_mb * 1024 * 1024;
    caler.tol = opts.conv;
    caler.max_iter = opts.max_iter;
    let t_grid = t_grid.elapsed();
    let nato_car = mole.basis().shls.nato("car");
    report.fields.push(Field::new("格点", format!("{} 径向 × {} 球面（共 {} 个）", caler.radial, caler.angular, caler.grids.len())));
    report.fields.push(Field::new("π 约化基", format!("{} / {}（笛卡尔原子轨道）", caler.pi_atos.len(), nato_car)));
    report.fields.push(Field::new("π 积分形式", "car（局部坐标系）"));

    let res = caler.run();
    let (ele_eng, xc_eng, nuc_eng) = (res.ele_eng, res.xc_eng, res.nuc_eng);
    let (converged, niter, delta_eng) = (res.converged, res.niter, res.delta_eng);
    let (time_int1e, time_int2e, time_scf) = (res.time_int1e, res.time_int2e, res.time_scf);

    report.fields.extend(calc_status_fields(converged, niter, delta_eng));
    report.fields.push(Field::new("π 电子数", res.npi_elec.to_string()));
    report.fields.push(Field::new("π 电子能量（模型）", format!("{ele_eng:>20.10} Hartree")));
    report.fields.push(Field::new("交换相关能", format!("{xc_eng:>19.10} Hartree")));
    report.fields.push(Field::new("核排斥能（完整分子）", format!("{nuc_eng:>18.10} Hartree")));

    // π 分子轨道能量
    let mut mo_notes = vec!["π 分子轨道能量（Hartree）".to_string()];
    let nocc = res.nocc();
    for (i, eng) in res.obt_engs.iter().enumerate() {
        let occ = if i < nocc { "occ" } else { "vir" };
        mo_notes.push(format!("  #{:<3} {:<4} {eng:>18.10}", i + 1, occ));
    }
    report.notes.extend(mo_notes);

    // 局部坐标系下的 π 分子轨道系数（可选，行标签只含 π 成分）
    if print_mo {
        report.notes.push("局部坐标系下的 π 分子轨道系数（行 = π 原子轨道，列 = π 分子轨道）".to_string());
        for line in res.mo_table_local(&mole).lines() {
            report.notes.push(line.to_string());
        }
    }

    // 导出为 fch / mwfn（系数已展开回全局笛卡尔基）
    let t_export = Instant::now();
    let saved = match output {
        Some(output) => {
            let dest = absolute(output).map_err(|err| format!("无法解析输出路径 {}: {err}", output.display()))?;
            let ext = ext_of(&dest).unwrap_or_default();
            let dest_path = dest.to_str().ok_or_else(|| format!("输出路径不是有效的 UTF-8: {}", dest.display()))?.to_string();
            let out_mole = res.to_mole(&mole);
            match ext.as_str() {
                "fch" | "fchk" => {
                    let mut writer = FchWriter::from_mole(&out_mole);
                    writer.energy = res.ele_eng;
                    writer.charge = out_mole.charge();
                    writer.multip = out_mole.multip();
                    writer.save(&dest_path);
                }
                "mwfn" | "mfn" => {
                    let mut writer = MfnWriter::from_mole(&out_mole);
                    writer.energy = res.ele_eng;
                    writer.save(&dest_path);
                }
                _ => return Err(format!("LKS 仅支持导出 .fch/.fchk 或 .mwfn，收到: {}", dest.display())),
            }
            Some(dest_path)
        }
        None => None,
    };
    let t_export = t_export.elapsed();

    report.notes.extend(timing_notes(t_load, Some(t_grid), time_int1e, time_int2e, time_scf, t_export, t_total.elapsed()));
    if let Some(path) = saved {
        report.notes.push(format!("已导出 π 波函数（系数已展开到全局基，电子数=π 电子数）到: {path}"));
    }
    report.notes.push("提示: LKS 是 π 子系统模型，能量不含 σ 电子与核内层，不是分子总能量".to_string());
    if !converged {
        report.notes.push(format!("警告: SCF 在 {niter} 次迭代后仍未收敛，结果可能不可靠"));
    }
    Ok(report)
}
