//! 格点属性任务：实空间函数导出、静电势、拓扑分析。

use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use rayon::prelude::*;

use crate::base::Mole;
use crate::props::gridprop::{CubeGrid, density, topolog};
use crate::report::{Column, Field, Report, Table};
use crate::utils;
use crate::writer::{CubWriter, Writer};

use super::{ServiceResult, report_for, require_coefs};

/// cub 写入结果摘要
struct CubeInfo {
    npts: usize,
    shape: [usize; 3],
    p0: [f64; 3],
    p1: [f64; 3],
    elapsed: Duration,
}

/// 实空间函数的中文/英文标题
fn density_title(func: &str) -> &'static str {
    match func {
        "density" => "Electron Density",
        "pro_rho" => "Promolecular Density",
        "spin" => "Spin Density (alpha-beta)",
        "rho_a" => "Alpha Electron Density",
        "rho_b" => "Beta Electron Density",
        "del_rho" => "Deformation Density (Delta rho)",
        "pi_pocv" => "Pi Electron Density (POCV)",
        "pi_mocv" => "Pi Electron Density (MOCV)",
        "rdg" => "Reduced Density Gradient (RDG)",
        "iri" => "Interaction Region Indicator (IRI)",
        "laplace" => "Laplacian of Electron Density",
        "signl2rho" => "sign(lambda2)*rho",
        "ked" => "Kinetic Energy Density",
        "lagkin" => "Lagrangian Kinetic Energy Density",
        "hamkin" => "Hamiltonian Kinetic Energy Density",
        "elf" => "Electron Localization Function (ELF)",
        "lol" => "Localized Orbital Locator (LOL)",
        "infoentro" => "Information Entropy Density",
        "activate" => "Activation Density",
        "curvature" => "Electron Density Curvature",
        "fukui_plus" => "Fukui Function f+ = rho(N+1) - rho(N)",
        "fukui_minus" => "Fukui Function f- = rho(N) - rho(N-1)",
        "fukui_zero" => "Fukui Function f0",
        "fukui_df" => "Dual Descriptor Df = f+ - f-",
        _ => "Scalar Field",
    }
}

/// 根据函数名生成格点标量场计算闭包
///
/// fukui_* 系列函数需要额外提供 N+1（阴离子）与 N-1（阳离子）波函数文件。
fn density_field<'a>(mole: &'a Mole, func: &str, mole_n: Option<&'a Mole>, mole_p: Option<&'a Mole>) -> ServiceResult<Box<dyn Fn(&[f64; 3]) -> f64 + Send + Sync + 'a>> {
    let caler = density::Calculator::new(mole);
    let field: Box<dyn Fn(&[f64; 3]) -> f64 + Send + Sync + 'a> = match func {
        "density" => Box::new(move |grid| caler.mol_rho_dm(grid, 0).val0),
        "pro_rho" => Box::new(move |grid| caler.pro_rho(grid)),
        "spin" => Box::new(move |grid| {
            let [rho_a, rho_b] = caler.mol_rho_spin(grid, 0);
            rho_a.val0 - rho_b.val0
        }),
        "rho_a" => Box::new(move |grid| caler.mol_rho_spin(grid, 0)[0].val0),
        "rho_b" => Box::new(move |grid| caler.mol_rho_spin(grid, 0)[1].val0),
        "pi_pocv" => Box::new(move |grid| caler.pi_rho_pocv(grid, 0).val0),
        "pi_mocv" => Box::new(move |grid| caler.pi_rho_mocv(grid, 0).val0),
        "del_rho" => Box::new(move |grid| caler.del_rho(grid)),
        "rdg" => Box::new(move |grid| caler.RDG(grid)),
        "iri" => Box::new(move |grid| caler.IRI(grid, None)),
        "laplace" => Box::new(move |grid| caler.laplace(grid)),
        "signl2rho" => Box::new(move |grid| caler.signl2rho(grid)),
        "ked" => Box::new(move |grid| caler.KED(grid)),
        "lagkin" => Box::new(move |grid| caler.lagkin(grid)[0]),
        "hamkin" => Box::new(move |grid| caler.hamkin(grid)),
        "elf" => Box::new(move |grid| caler.ELF(grid)),
        "lol" => Box::new(move |grid| caler.LOL(grid)),
        "infoentro" => Box::new(move |grid| caler.infoentro(grid)),
        "activate" => Box::new(move |grid| caler.activate(grid)),
        "curvature" => Box::new(move |grid| caler.curvature(grid)),
        // 网格版福井函数：Δf 与原子级分析保持一致（f⁺ - f⁻）
        "fukui_plus" => {
            let mole_n = mole_n.ok_or_else(|| "fukui_plus 需要 --mole-n 指定阴离子（N+1）波函数文件".to_string())?;
            let caler_n = density::Calculator::new(mole_n);
            Box::new(move |grid| caler_n.mol_rho_dm(grid, 0).val0 - caler.mol_rho_dm(grid, 0).val0)
        }
        "fukui_minus" => {
            let mole_p = mole_p.ok_or_else(|| "fukui_minus 需要 --mole-p 指定阳离子（N-1）波函数文件".to_string())?;
            let caler_p = density::Calculator::new(mole_p);
            Box::new(move |grid| caler.mol_rho_dm(grid, 0).val0 - caler_p.mol_rho_dm(grid, 0).val0)
        }
        "fukui_zero" => {
            let mole_n = mole_n.ok_or_else(|| "fukui_zero 需要 --mole-n 指定阴离子（N+1）波函数文件".to_string())?;
            let mole_p = mole_p.ok_or_else(|| "fukui_zero 需要 --mole-p 指定阳离子（N-1）波函数文件".to_string())?;
            let caler_n = density::Calculator::new(mole_n);
            let caler_p = density::Calculator::new(mole_p);
            Box::new(move |grid| (caler_n.mol_rho_dm(grid, 0).val0 - caler_p.mol_rho_dm(grid, 0).val0) / 2.0)
        }
        "fukui_df" => {
            let mole_n = mole_n.ok_or_else(|| "fukui_df 需要 --mole-n 指定阴离子（N+1）波函数文件".to_string())?;
            let mole_p = mole_p.ok_or_else(|| "fukui_df 需要 --mole-p 指定阳离子（N-1）波函数文件".to_string())?;
            let caler_n = density::Calculator::new(mole_n);
            let caler_p = density::Calculator::new(mole_p);
            Box::new(move |grid| {
                let rho0 = caler.mol_rho_dm(grid, 0).val0;
                let rho_n = caler_n.mol_rho_dm(grid, 0).val0;
                let rho_p = caler_p.mol_rho_dm(grid, 0).val0;
                rho_n + rho_p - 2.0 * rho0
            })
        }
        other => {
            return Err(format!(
                "未知的实空间函数: {other}（可用: density/pro_rho/spin/rho_a/rho_b/del_rho/pi_pocv/pi_mocv/rdg/iri/laplace/signl2rho/ked/lagkin/hamkin/elf/lol/infoentro/activate/curvature/fukui_plus/fukui_minus/fukui_zero/fukui_df）"
            ));
        }
    };
    Ok(field)
}

/// 计算分子边界 + 边框
fn molecule_bounding_box(mole: &Mole, border: f64) -> ([f64; 3], [f64; 3]) {
    let xyzs = mole.atoms().xyzs();
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for xyz in xyzs.iter() {
        for i in 0..3 {
            if xyz[i] < min[i] {
                min[i] = xyz[i];
            }
            if xyz[i] > max[i] {
                max[i] = xyz[i];
            }
        }
    }
    (min.map(|v| v - border), max.map(|v| v + border))
}

/// 生成默认的 cub 输出路径: 输入文件名_{func}.cub
fn default_cub_path(mole: &Mole, func: &str) -> ServiceResult<PathBuf> {
    let path = mole.path.as_deref().ok_or_else(|| "分子没有来源文件路径，无法生成默认 cub 文件名".to_string())?;
    Ok(PathBuf::from(format!("{}_{}.cub", utils::file_name(path), func)))
}

/// 在立方体格点上计算标量场值并写入 cub 文件
fn write_cube<F>(mole: &Mole, output: &Path, step: f64, title: &str, field: F) -> ServiceResult<CubeInfo>
where
    F: Fn(&[f64; 3]) -> f64 + Send + Sync,
{
    let started = Instant::now();
    let border = 4.0;
    let (p0, p1) = molecule_bounding_box(mole, border);

    let mut cube_grid = CubeGrid::new();
    cube_grid.set_v1(&p0, &p1, step, 0.0);
    let (shape, grids) = cube_grid.get();
    let start = cube_grid.start;
    let steps = cube_grid.steps;
    let [nx, ny, nz] = shape;
    if nx == 0 || ny == 0 || nz == 0 {
        return Err(format!("格点数为零，请减小步长或检查分子坐标（shape = {nx} x {ny} x {nz}）"));
    }

    // 使用 rayon 并行计算每个格点的值
    let vals: Vec<f64> = grids.par_iter().map(|grid| field(grid)).collect();

    let mut writer = CubWriter::new();
    writer.read_mole(mole);
    writer.title = title.to_string();
    writer.pos0 = start;
    writer.size = shape;
    writer.step = steps;
    writer.obts = vec![0];
    writer.vals = vec![vals];
    let output_str = output.to_str().ok_or_else(|| format!("输出路径不是有效的 UTF-8: {}", output.display()))?;
    writer.save(output_str);

    Ok(CubeInfo {
        npts: grids.len(),
        shape,
        p0,
        p1,
        elapsed: started.elapsed(),
    })
}

fn cube_report(mole: &Mole, func: &str, output: &Path, step: f64, info: CubeInfo) -> Report {
    let [nx, ny, nz] = info.shape;
    let [x0, y0, z0] = info.p0;
    let [x1, y1, z1] = info.p1;

    let mut report = report_for(format!("实空间函数计算: {} ({func})", density_title(func)), mole);
    report.fields.push(Field::new("函数", func.to_string()));
    report.fields.push(Field::new("输出文件", output.display().to_string()));
    report.fields.push(Field::new("格点数", format!("{}（{} × {} × {}）", info.npts, nx, ny, nz)));
    report.fields.push(Field::new("步长", format!("{step:.3} Bohr")));
    report.fields.push(Field::new("空间范围", format!("[{x0:.3}, {x1:.3}] × [{y0:.3}, {y1:.3}] × [{z0:.3}, {z1:.3}] Bohr")));
    report.fields.push(Field::new("耗时", format!("{:.3} s", info.elapsed.as_secs_f64())));
    report.notes.push("cub 文件已保存".to_string());
    report
}

/// 实空间函数（含 fukui_* 系列）导出 cub
pub fn density(mole: &Mole, output: Option<&Path>, step: f64, func: &str, mole_n: Option<&Mole>, mole_p: Option<&Mole>) -> ServiceResult<Report> {
    require_coefs(mole)?;
    if !(step.is_finite() && step > 0.0) {
        return Err(format!("格点步长必须为正数，收到: {step}"));
    }
    let field = density_field(mole, func, mole_n, mole_p)?;
    let output = match output {
        Some(path) => path.to_path_buf(),
        None => default_cub_path(mole, func)?,
    };
    let title = format!("{} ({func})", density_title(func));
    let info = write_cube(mole, &output, step, &title, field)?;
    Ok(cube_report(mole, func, &output, step, info))
}

/// 静电势（原子单位，核项 + 电子项）导出 cub
pub fn potential(mole: &Mole, output: Option<&Path>, step: f64) -> ServiceResult<Report> {
    require_coefs(mole)?;
    if !(step.is_finite() && step > 0.0) {
        return Err(format!("格点步长必须为正数，收到: {step}"));
    }
    let caler = density::Calculator::new(mole);
    let output = match output {
        Some(path) => path.to_path_buf(),
        None => default_cub_path(mole, "potential")?,
    };
    let info = write_cube(mole, &output, step, "Electrostatic Potential", |grid| caler.tot_esp(grid))?;
    Ok(cube_report(mole, "potential", &output, step, info))
}

/// 拓扑分析（寻找临界点）
pub fn topology(mole: &Mole, func: &str) -> ServiceResult<Report> {
    require_coefs(mole)?;
    if !matches!(func, "mol_rho_dm" | "pi_rho_pocv" | "pi_rho_mocv") {
        return Err(format!("未知的拓扑分析函数: {func}（可用: mol_rho_dm | pi_rho_pocv | pi_rho_mocv）"));
    }
    let caler = topolog::Calculator::new(mole);
    let cpts = caler.analysis(func);

    let mut report = report_for(format!("拓扑分析（函数: {func}）"), mole);
    report.fields.push(Field::new("临界点数", cpts.len().to_string()));

    let mut table = Table::new(vec![
        Column::right("序号"),
        Column::left("临界点类型"),
        Column::right("X (Bohr)"),
        Column::right("Y (Bohr)"),
        Column::right("Z (Bohr)"),
        Column::right("函数值"),
    ]);
    for (i, (pos, val, cpt)) in cpts.iter().enumerate() {
        let cpt_name = match cpt {
            -3 => "(3, -3) 核临界点",
            -1 => "(3, -1) 键临界点",
            1 => "(3, +1) 环临界点",
            3 => "(3, +3) 笼临界点",
            _ => "其他",
        };
        table.push_row(vec![
            (i + 1).to_string(),
            cpt_name.to_string(),
            format!("{:.12}", pos[0]),
            format!("{:.12}", pos[1]),
            format!("{:.12}", pos[2]),
            format!("{:.8}", val),
        ]);
    }
    report.tables.push(table);
    if cpts.is_empty() {
        report.notes.push("未找到临界点".to_string());
    }
    Ok(report)
}
