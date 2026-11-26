//! 量子化学计算任务（HF / DFT）。

use std::path::{Path, absolute};
use std::time::{Duration, Instant};

use crate::base::Mole;
use crate::qmcal;
use crate::report::{Field, Report};
use crate::rsbse;
use crate::rsxcs::{XcFunc, supported_names};
use crate::writer::{FchWriter, Writer};

use super::{ServiceResult, ext_of};

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

    let caler = qmcal::hf::RHF::new(&mole);
    let res = caler.run(None);
    let total_eng = res.total_eng();
    let (ele_eng, nuc_eng, converged, niter, delta_eng) = (res.ele_eng, res.nuc_eng, res.converged, res.niter, res.delta_eng);
    let (time_int1e, time_int2e, time_scf) = (res.time_int1e, res.time_int2e, res.time_scf);

    report.fields.extend(calc_status_fields(converged, niter, delta_eng));
    report.fields.push(Field::new("电子能量", format!("{ele_eng:>20.10} Hartree")));
    report.fields.push(Field::new("核排斥能", format!("{nuc_eng:>20.10} Hartree")));
    report.fields.push(Field::new("总能量", format!("{total_eng:>20.10} Hartree")));

    let t_export = Instant::now();
    let saved = write_fch(&mole, output, res.basis, res.coefs, total_eng)?;
    let t_export = t_export.elapsed();

    report.notes.extend(timing_notes(t_load, None, time_int1e, time_int2e, time_scf, t_export, t_total.elapsed()));
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
        let caler = qmcal::ks::RKS::new(&mole, xc);
        let t_grid = t_grid.elapsed();
        report.fields.push(Field::new("格点", format!("{} 径向 × {} 球面（共 {} 个）", caler.radial, caler.angular, caler.grids.len())));

        let res = caler.run();
        let total_eng = res.total_eng();
        let (ele_eng, xc_eng, nuc_eng) = (res.ele_eng, res.xc_eng, res.nuc_eng);
        let (converged, niter, delta_eng) = (res.converged, res.niter, res.delta_eng);
        let (time_int1e, time_int2e, time_scf) = (res.time_int1e, res.time_int2e, res.time_scf);

        report.fields.extend(calc_status_fields(converged, niter, delta_eng));
        report.fields.push(Field::new("电子能量", format!("{ele_eng:>20.10} Hartree")));
        report.fields.push(Field::new("交换相关能", format!("{xc_eng:>19.10} Hartree")));
        report.fields.push(Field::new("核排斥能", format!("{nuc_eng:>20.10} Hartree")));
        report.fields.push(Field::new("总能量", format!("{total_eng:>20.10} Hartree")));

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
        let caler = qmcal::ks::UKS::new(&mole, xc);
        let t_grid = t_grid.elapsed();
        report
            .fields
            .push(Field::new("格点", format!("{} 径向 × {} 球面（共 {} 个）", caler.base.radial, caler.base.angular, caler.base.grids.len())));

        let res = caler.run();
        let total_eng = res.total_eng();
        let (ele_eng, xc_eng, nuc_eng) = (res.ele_eng, res.xc_eng, res.nuc_eng);
        let (converged, niter, delta_eng) = (res.converged, res.niter, res.delta_eng);
        let (time_int1e, time_int2e, time_scf) = (res.time_int1e, res.time_int2e, res.time_scf);

        report.fields.extend(calc_status_fields(converged, niter, delta_eng));
        report.fields.push(Field::new("电子能量", format!("{ele_eng:>20.10} Hartree")));
        report.fields.push(Field::new("交换相关能", format!("{xc_eng:>19.10} Hartree")));
        report.fields.push(Field::new("核排斥能", format!("{nuc_eng:>20.10} Hartree")));
        report.fields.push(Field::new("总能量", format!("{total_eng:>20.10} Hartree")));
        report.fields.push(Field::new("<S^2>", format!("{:.6}（多重度 {}）", res.spin_sq, res.multip())));

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
