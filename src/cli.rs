use clap::{Parser, Subcommand, command};
use rayon::prelude::*;
use std::collections::HashMap;
use std::path::{PathBuf, absolute};

use crate::atomprop;
use crate::atomprop::activity;
use crate::base::{Mole, Vec3};
use crate::bondprop;
use crate::dev;
use crate::gridprop::{CubeGrid, density, potential, topolog};
use crate::moleprop;
use crate::orbtprop;
use crate::reader::{FchReader, GjfReader, LogReader};
use crate::tools::log::Tool;
use crate::writer::{CubWriter, DmoWriter, GjfWriter, Writer, XyzWriter};
use crate::{qmcalc, utils};

#[derive(Parser, Debug)]
#[command(
    author = "Xiaofei Shi",
    version,
    about = "Weave Function Analysis Tool"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// 原子属性分析
    #[command(subcommand)]
    Atomprop(AtompropCmds),
    /// 键属性分析
    #[command(subcommand)]
    Bondprop(BondpropCmds),
    /// 分子属性分析
    #[command(subcommand)]
    Moleprop(MolepropCmds),
    /// 轨道属性分析
    #[command(subcommand)]
    Orbtprop(OrbtpropCmds),
    /// 格点属性分析，导出 cub 文件
    #[command(subcommand)]
    Gridprop(GridpropCmds),

    /// 提取分子坐标（打印到屏幕，或导出 .xyz 文件）
    Coord {
        /// 输入文件路径（支持 log/out/fch/fchk/gjf/xyz/mol/sdf/molden/dmo/inp 等格式）
        input: PathBuf,
        /// 输出 .xyz 文件路径（可选，默认只打印坐标）
        output: Option<PathBuf>,
        /// 打印坐标单位: ang（埃） | bohr（原子单位），默认 ang
        #[arg(long, default_value = "ang", value_parser = ["ang", "bohr"])]
        unit: String,
    },

    /// 文件格式转换
    Trans { input: PathBuf, output: PathBuf },
    /// DHMO 计算
    Dmo { input: PathBuf, output: Option<PathBuf> },

    /// 提取支持信息（Supporting Information），汇总能量、频率与结构
    Si {
        /// 输入 Gaussian 输出文件路径（.log/.out），可指定多个
        inputs: Vec<PathBuf>,
        /// 输出文本文件路径（可选，默认打印到屏幕）
        #[arg(short = 'o', long = "output")]
        output: Option<PathBuf>,
    },

    /// 开发时测试
    Dev,
}

// ===== atomprop 子命令 =====

#[derive(Subcommand, Debug)]
enum AtompropCmds {
    /// 原子电荷/布居分析
    #[command(subcommand)]
    Charge(ChargeCmds),
    /// 原子π电子布居（POCV/MOCV 轨道改造方法）
    PiPop {
        /// 波函数文件路径（支持 .log/.out/.fch/.fchk/.gjf）
        filepath: PathBuf,
        /// 布居方法: mulliken | lowdin | hirshfeld
        #[arg(long, default_value = "mulliken", value_parser = ["mulliken", "lowdin", "hirshfeld"])]
        ctype: String,
        /// 轨道改造方法: pocv | mocv
        #[arg(long, default_value = "pocv", value_parser = ["pocv", "mocv"])]
        method: String,
    },
    /// 福井函数与双描述符（需要 N、N+1、N-1 三个波函数文件）
    Fukui {
        /// 中性分子文件（N 电子）
        filepath: PathBuf,
        /// 阴离子文件（N+1 电子）
        anion: PathBuf,
        /// 阳离子文件（N-1 电子）
        cation: PathBuf,
        /// 布居方法: mulliken | lowdin | hirshfeld
        #[arg(long, default_value = "mulliken", value_parser = ["mulliken", "lowdin", "hirshfeld"])]
        ctype: String,
    },
    /// π轨道的福井函数与双描述符（POCV/MOCV 方法）
    FukuiPi {
        /// 中性分子文件（N 电子）
        filepath: PathBuf,
        /// 阴离子文件（N+1 电子）
        anion: PathBuf,
        /// 阳离子文件（N-1 电子）
        cation: PathBuf,
        /// 布居方法: mulliken | lowdin | hirshfeld
        #[arg(long, default_value = "mulliken", value_parser = ["mulliken", "lowdin", "hirshfeld"])]
        ctype: String,
        /// 轨道改造方法: pocv | mocv
        #[arg(long, default_value = "pocv", value_parser = ["pocv", "mocv"])]
        method: String,
    },
    /// 自由价分析（free valence）
    Freev {
        /// 波函数文件路径
        filepath: PathBuf,
        /// 键级类型: mayer | wiberg
        #[arg(long, default_value = "mayer", value_parser = ["mayer", "wiberg"])]
        otype: String,
    },
    /// π自由价分析（POCV 方法）
    FreevPi {
        /// 波函数文件路径
        filepath: PathBuf,
        /// 键级类型: mayer | wiberg
        #[arg(long, default_value = "mayer", value_parser = ["mayer", "wiberg"])]
        otype: String,
    },
}

#[derive(Subcommand, Debug)]
enum ChargeCmds {
    /// Mulliken 布居分析
    Mulliken {
        /// 波函数文件路径（支持 .log/.out/.fch/.fchk/.gjf）
        filepath: PathBuf,
    },
    /// Löwdin 布居分析
    Lowdin {
        /// 波函数文件路径
        filepath: PathBuf,
    },
    /// Hirshfeld 电荷分析
    Hirshfeld {
        /// 波函数文件路径
        filepath: PathBuf,
    },
    /// 自旋密度布居分析
    Spin {
        /// 波函数文件路径
        filepath: PathBuf,
        /// 电荷类型: mulliken | lowdin | hirshfeld
        #[arg(long, default_value = "mulliken", value_parser = ["mulliken", "lowdin", "hirshfeld"])]
        ctype: String,
    },
}

// ===== bondprop 子命令 =====

#[derive(Subcommand, Debug)]
enum BondpropCmds {
    /// 键级分析
    #[command(subcommand)]
    Order(OrderCmds),
    /// π键级分析（POCV/MOCV 方法）
    PiOrder {
        /// 波函数文件路径
        filepath: PathBuf,
        /// 轨道改造方法: pocv | mocv
        #[arg(long, default_value = "pocv", value_parser = ["pocv", "mocv"])]
        method: String,
        /// 键级类型: mayer | wiberg
        #[arg(long, default_value = "mayer", value_parser = ["mayer", "wiberg"])]
        otype: String,
    },
    /// 键级分解，将两个原子间的键级分解为 σ、π、δ、φ 分量
    Decom {
        /// 原子1 序号（从1开始）
        atm1: usize,
        /// 原子2 序号（从1开始）
        atm2: usize,
        /// 波函数文件路径
        filepath: PathBuf,
    },
}

#[derive(Subcommand, Debug)]
enum OrderCmds {
    /// Mayer 键级
    Mayer {
        /// 波函数文件路径
        filepath: PathBuf,
    },
    /// Wiberg 键级
    Wiberg {
        /// 波函数文件路径
        filepath: PathBuf,
    },
    /// 指定分子轨道对键级的贡献
    Obt {
        /// 分子轨道序号（从1开始）
        obt: usize,
        /// 波函数文件路径
        filepath: PathBuf,
        /// 键级类型: mayer | wiberg
        #[arg(long, default_value = "mayer", value_parser = ["mayer", "wiberg"])]
        otype: String,
    },
}

// ===== moleprop 子命令 =====

#[derive(Subcommand, Debug)]
enum MolepropCmds {
    /// 芳香性分析（基于π键级，PISD 方法）
    Aromacity {
        /// 波函数文件路径
        filepath: PathBuf,
        /// 环原子列表（序号从1开始，如 "1-6" 或 "1,2,3,4,5,6"），可重复指定；未指定时自动检测分子中的环
        #[arg(long = "ring")]
        rings: Vec<String>,
    },
}

// ===== orbtprop 子命令 =====

#[derive(Subcommand, Debug)]
enum OrbtpropCmds {
    /// 每个分子轨道的π电子数（POCV 方法）
    PiElesPocv {
        /// 波函数文件路径
        filepath: PathBuf,
    },
    /// 每个分子轨道的π电子数（MOCV 方法）
    PiElesMocv {
        /// 波函数文件路径
        filepath: PathBuf,
    },
    /// 指定分子轨道的原子电子布居（轨道成分分析）
    Composition {
        /// 分子轨道序号（从1开始）
        obt: usize,
        /// 波函数文件路径
        filepath: PathBuf,
        /// 布居方法: mulliken | lowdin | hirshfeld
        #[arg(long, default_value = "mulliken", value_parser = ["mulliken", "lowdin", "hirshfeld"])]
        ctype: String,
    },
}

// ===== gridprop 子命令 =====

#[derive(Subcommand, Debug)]
enum GridpropCmds {
    /// 计算实空间函数格点值，输出 cub 文件
    Density {
        /// 波函数文件路径（支持 .log/.out/.fch/.fchk）
        filepath: PathBuf,
        /// 输出 cub 文件路径（默认: 输入文件名_{func}.cub）
        output: Option<PathBuf>,
        /// 格点步长，单位 Bohr（默认 0.2）
        #[arg(default_value = "0.2", long = "step")]
        step: f64,
        /// 实空间函数名: density | pro_rho | spin | rho_a | rho_b | del_rho | pi_pocv | pi_mocv
        /// | rdg | iri | laplace | signl2rho | ked | lagkin | hamkin | elf | lol
        /// | infoentro | activate | curvature | fukui_plus | fukui_minus | fukui_zero | fukui_df
        #[arg(default_value = "density", long = "func")]
        func: String,
        /// fukui_* 函数需要: 阴离子(N+1)波函数文件
        #[arg(long = "mole-n")]
        mole_n: Option<PathBuf>,
        /// fukui_* 函数需要: 阳离子(N-1)波函数文件
        #[arg(long = "mole-p")]
        mole_p: Option<PathBuf>,
    },
    /// 静电势，输出 cub 文件
    Potential {
        /// 波函数文件路径（支持 .log/.out/.fch/.fchk）
        filepath: PathBuf,
        /// 输出 cub 文件路径（默认: 输入文件名_potential.cub）
        output: Option<PathBuf>,
        /// 格点步长，单位 Bohr（默认 0.2）
        #[arg(default_value = "0.2", long = "step")]
        step: f64,
    },
    /// 拓扑分析，寻找临界点
    Topology {
        /// 函数名: mol_rho_dm | pi_rho_pocv | pi_rho_mocv
        #[arg(short = 'f')]
        func: String,
        /// 波函数文件路径（支持 .log/.out/.fch/.fchk）
        filepath: PathBuf,
    },
}

// ===== 读取分子文件 =====

/// 根据文件路径自动选择合适的读取器加载分子
fn load_mole(filepath: &PathBuf) -> Mole {
    let path = absolute(filepath).expect("无法获取文件绝对路径");
    let path_str = path.to_str().expect("路径转字符串失败");
    let ftype = utils::file_type(path_str);
    match ftype.as_str() {
        "log" | "out" => Mole::from_read(&LogReader::new(path_str)),
        "gjf" => Mole::from_read(&GjfReader::new(path_str)),
        "fch" | "fchk" => Mole::from_read(&FchReader::new(path_str)),
        _ => panic!("不支持的波函数文件类型: {ftype}，支持: log/out/gjf/fch/fchk"),
    }
}

// ===== 打印辅助函数 =====

/// 打印原子布居结果表
fn print_atom_population(mole: &Mole, pops: &[f64], title: &str) {
    let atoms = mole.atoms();
    println!("\n{}", title);
    println!("{}", "=".repeat(title.len()));
    println!("{:>4}  {:>6}  {:>14}", "原子", "符号", "电子布居");
    println!("{:-<28}", "");
    for (i, &pop) in pops.iter().enumerate() {
        println!("{:>4}  {:>6}  {:>14.8}", i + 1, atoms.atoms[i].sym, pop);
    }
    let total: f64 = pops.iter().sum();
    println!("{:-<28}", "");
    println!("{:>12}  {:>14.8}", "总电子数", total);
}

/// 打印键级矩阵（仅上三角，筛选显著键级）
fn print_bond_order(mole: &Mole, omat: &ndarray::Array2<f64>, title: &str) {
    let natm = mole.atoms().len();
    let atoms = mole.atoms();
    println!("\n{}", title);
    println!("{}", "=".repeat(title.len()));
    println!("{:>6}  {:>6}  {:>12}", "原子对", "序号", "键级");
    println!("{:-<30}", "");
    for i in 0..natm {
        for j in (i + 1)..natm {
            let order = omat[(i, j)];
            if order.abs() > 0.01 {
                println!(
                    "{:>6}  {:>3}-{:<3}  {:>12.6}",
                    format!("{}-{}", atoms.atoms[i].sym, atoms.atoms[j].sym),
                    i + 1,
                    j + 1,
                    order
                );
            }
        }
    }
}

/// 打印福井函数结果表（7 列: e0, e(N+1), e(N-1), f⁻, f⁺, f⁰, Δf）
fn print_fukui_table(mole: &Mole, vals: &[[f64; 7]], title: &str) {
    let atoms = mole.atoms();
    println!("\n{}", title);
    println!("{}", "=".repeat(title.len()));
    println!(
        "{:>4}  {:>6}  {:>12}  {:>12}  {:>12}  {:>10}  {:>10}  {:>10}  {:>10}",
        "原子", "符号", "e0(N)", "e(N+1)", "e(N-1)", "f⁻", "f⁺", "f⁰", "Δf"
    );
    println!("{:-<96}", "");
    for (i, v) in vals.iter().enumerate() {
        println!(
            "{:>4}  {:>6}  {:>12.8}  {:>12.8}  {:>12.8}  {:>10.6}  {:>10.6}  {:>10.6}  {:>10.6}",
            i + 1,
            atoms.atoms[i].sym,
            v[0],
            v[1],
            v[2],
            v[3],
            v[4],
            v[5],
            v[6]
        );
    }
}

/// 打印自由价结果表（3 列: 标准价, 键级和, 自由价）
fn print_freev_table(mole: &Mole, vals: &[[f64; 3]], title: &str) {
    let atoms = mole.atoms();
    println!("\n{}", title);
    println!("{}", "=".repeat(title.len()));
    println!(
        "{:>4}  {:>6}  {:>14}  {:>14}  {:>14}",
        "原子", "符号", "标准价", "键级和", "自由价"
    );
    println!("{:-<56}", "");
    for (i, v) in vals.iter().enumerate() {
        println!(
            "{:>4}  {:>6}  {:>14.6}  {:>14.6}  {:>14.6}",
            i + 1,
            atoms.atoms[i].sym,
            v[0],
            v[1],
            v[2]
        );
    }
}

/// 打印π方向向量表
fn print_pi_dirs(dirs: &HashMap<usize, Vec3>, mole: &Mole) {
    let syms = mole.atoms().syms();
    let mut idxs: Vec<usize> = dirs.keys().cloned().collect();
    idxs.sort_unstable();
    println!("{:>4}  {:>6}  {:>18}  {:>18}  {:>18}", "原子", "符号", "X", "Y", "Z");
    println!("{:-<70}", "");
    for idx in idxs {
        let [x, y, z] = dirs[&idx].into_arr();
        println!(
            "{:>4}  {:>6}  {:>18.10}  {:>18.10}  {:>18.10}",
            idx + 1,
            syms[idx],
            x,
            y,
            z
        );
    }
}

/// 解析形如 "1-3,5,7-9" 的原子列表（1基、含逗号/分号/空格分隔与区间），返回 0 基索引
fn parse_atom_list(spec: &str, natm: usize) -> Vec<usize> {
    let mut idxs = vec![];
    for token in spec.split(|c: char| c == ',' || c == ';' || c.is_whitespace()) {
        let token = token.trim();
        if token.is_empty() {
            continue;
        }
        let nums: Vec<usize> = token
            .split('-')
            .map(|s| s.trim().parse().unwrap_or_else(|_| panic!("无法解析原子序号: {spec}")))
            .collect();
        let (a, b) = match nums[..] {
            [a] => (a, a),
            [a, b] => (a, b),
            _ => panic!("无法解析原子序号区间: {spec}"),
        };
        if a == 0 || b == 0 || a > b || b > natm {
            panic!("原子序号越界: {spec}（分子共有 {natm} 个原子，序号从 1 开始）");
        }
        for i in a..=b {
            idxs.push(i - 1);
        }
    }
    idxs
}

/// 计算分子边界 + 4 Bohr 边框
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
    let p0 = min.map(|v| v - border);
    let p1 = max.map(|v| v + border);
    (p0, p1)
}

/// 生成默认的 cub 输出路径: 输入文件名_{func}.cub
fn default_cub_path(filepath: &PathBuf, func: &str) -> String {
    let path = absolute(filepath).expect("无法获取文件绝对路径");
    let path_str = path.to_str().expect("路径转字符串失败");
    let name = utils::file_name(path_str);
    format!("{}_{}.cub", name, func)
}

/// 在立方体格点上计算标量场值，并写入 cub 文件
fn compute_and_write_cube<F>(mole: &Mole, output_path: &str, step: f64, title: &str, field_calc: F)
where
    F: Fn(&[f64; 3]) -> f64 + Send + Sync,
{
    // 计算分子边界 + 4 Bohr 边框
    let border = 4.0;
    let (p0, p1) = molecule_bounding_box(mole, border);

    // 创建立方体格点
    let mut cube_grid = CubeGrid::new();
    cube_grid.set_v1(&p0, &p1, step, 0.0);
    let (shape, grids) = cube_grid.get();
    let [nx, ny, nz] = shape;
    let start = cube_grid.start;
    let steps = cube_grid.steps;

    let npts = grids.len();
    println!("格点数: {} ({} x {} x {}), 步长: {:.3} Bohr", npts, nx, ny, nz, step);
    println!(
        "空间范围: [{:.3}, {:.3}] x [{:.3}, {:.3}] x [{:.3}, {:.3}] (Bohr)",
        p0[0], p1[0], p0[1], p1[1], p0[2], p1[2]
    );
    println!("正在计算 {} ...", title);

    // 使用 rayon 并行计算每个格点的值
    let vals: Vec<f64> = grids.par_iter().map(|grid| field_calc(grid)).collect();

    // 写入 cub 文件
    let mut writer = CubWriter::new();
    writer.read_mole(mole);
    writer.title = title.to_string();
    writer.pos0 = start;
    writer.size = [nx, ny, nz];
    writer.step = steps;
    writer.obts = vec![0];
    writer.vals = vec![vals];
    writer.save(output_path);

    println!("cub 文件已保存: {}", output_path);
}

// ===== 实空间函数 =====

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
fn density_field<'a>(
    mole: &'a Mole,
    func: &str,
    mole_n: Option<&'a Mole>,
    mole_p: Option<&'a Mole>,
) -> Box<dyn Fn(&[f64; 3]) -> f64 + Send + Sync + 'a> {
    let caler = density::Calculator::new(mole);
    match func {
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
            let caler_n =
                density::Calculator::new(mole_n.expect("fukui_plus 需要 --mole-n 指定阴离子（N+1）波函数文件"));
            Box::new(move |grid| caler_n.mol_rho_dm(grid, 0).val0 - caler.mol_rho_dm(grid, 0).val0)
        }
        "fukui_minus" => {
            let caler_p =
                density::Calculator::new(mole_p.expect("fukui_minus 需要 --mole-p 指定阳离子（N-1）波函数文件"));
            Box::new(move |grid| caler.mol_rho_dm(grid, 0).val0 - caler_p.mol_rho_dm(grid, 0).val0)
        }
        "fukui_zero" => {
            let caler_n =
                density::Calculator::new(mole_n.expect("fukui_zero 需要 --mole-n 指定阴离子（N+1）波函数文件"));
            let caler_p =
                density::Calculator::new(mole_p.expect("fukui_zero 需要 --mole-p 指定阳离子（N-1）波函数文件"));
            Box::new(move |grid| (caler_n.mol_rho_dm(grid, 0).val0 - caler_p.mol_rho_dm(grid, 0).val0) / 2.0)
        }
        "fukui_df" => {
            let caler_n = density::Calculator::new(mole_n.expect("fukui_df 需要 --mole-n 指定阴离子（N+1）波函数文件"));
            let caler_p = density::Calculator::new(mole_p.expect("fukui_df 需要 --mole-p 指定阳离子（N-1）波函数文件"));
            Box::new(move |grid| {
                let rho0 = caler.mol_rho_dm(grid, 0).val0;
                let rho_n = caler_n.mol_rho_dm(grid, 0).val0;
                let rho_p = caler_p.mol_rho_dm(grid, 0).val0;
                rho_n + rho_p - 2.0 * rho0
            })
        }
        _ => panic!("未知的实空间函数: {func}，请查看 --help 中的函数列表"),
    }
}

// ===== 主入口 =====

/// CLI 模式入口函数，解析命令行参数执行对应操作
pub fn main() {
    let cli = Cli::parse();
    match &cli.command {
        // ==================== Atomprop ====================
        Commands::Atomprop(cmd) => match cmd {
            AtompropCmds::Charge(cmd) => match cmd {
                ChargeCmds::Mulliken { filepath } => {
                    let mole = load_mole(filepath);
                    let caler = atomprop::charge::Calculator::new(&mole);
                    let pops = caler.mulliken();
                    print_atom_population(&mole, &pops, "Mulliken 电子布居分析");
                }
                ChargeCmds::Lowdin { filepath } => {
                    let mole = load_mole(filepath);
                    let caler = atomprop::charge::Calculator::new(&mole);
                    let pops = caler.lowdin();
                    print_atom_population(&mole, &pops, "Löwdin 电子布居分析");
                }
                ChargeCmds::Hirshfeld { filepath } => {
                    let mole = load_mole(filepath);
                    let caler = atomprop::charge::Calculator::new(&mole);
                    let pops = caler.hirshfeld();
                    print_atom_population(&mole, &pops, "Hirshfeld 电子布居分析");
                }
                ChargeCmds::Spin { ctype, filepath } => {
                    let mole = load_mole(filepath);
                    if !mole.is_open() {
                        panic!("非开壳层分子（多重度=1）无法计算自旋布居");
                    }
                    let caler = atomprop::charge::Calculator::new(&mole);
                    let spin = caler.spin(ctype);
                    println!("\n自旋密度布居分析（ctype: {ctype}）");
                    println!("{}", "=".repeat(58));
                    println!(
                        "{:>4}  {:>6}  {:>14}  {:>14}  {:>14}",
                        "原子", "符号", "α电子", "β电子", "α - β"
                    );
                    println!("{:-<58}", "");
                    for (i, &[a, b, d]) in spin.iter().enumerate() {
                        println!(
                            "{:>4}  {:>6}  {:>14.8}  {:>14.8}  {:>14.8}",
                            i + 1,
                            mole.atoms().atoms[i].sym,
                            a,
                            b,
                            d
                        );
                    }
                }
            },
            AtompropCmds::PiPop {
                filepath,
                ctype,
                method,
            } => {
                let mole = load_mole(filepath);
                let caler = atomprop::charge::Calculator::new(&mole);
                match method.as_str() {
                    "pocv" => {
                        let (_dirs, eles) = caler.pi_pocv(ctype).expect("π-POCV 投影失败");
                        print_atom_population(&mole, &eles, &format!("原子π电子布居（POCV，{ctype}）"));
                    }
                    "mocv" => {
                        let (_stms, _atos, eles) = caler.pi_mocv(ctype).expect("π-MOCV 投影失败");
                        print_atom_population(&mole, &eles, &format!("原子π电子布居（MOCV，{ctype}）"));
                    }
                    _ => unreachable!(),
                }
            }
            AtompropCmds::Fukui {
                filepath,
                anion,
                cation,
                ctype,
            } => {
                let mole = load_mole(filepath);
                let mole_n = load_mole(anion);
                let mole_p = load_mole(cation);
                let caler = activity::Calculator::new(&mole);
                let vals = caler.fukui(&mole_n, &mole_p, ctype);
                print_fukui_table(&mole, &vals, &format!("福井函数与双描述符分析（{ctype}）"));
            }
            AtompropCmds::FukuiPi {
                filepath,
                anion,
                cation,
                ctype,
                method,
            } => {
                let mole = load_mole(filepath);
                let mole_n = load_mole(anion);
                let mole_p = load_mole(cation);
                let caler = activity::Calculator::new(&mole);
                match method.as_str() {
                    "pocv" => {
                        let (dirs, vals) = caler
                            .fukui_pi_pocv(&mole_n, &mole_p, ctype)
                            .expect("π福井函数（POCV）计算失败");
                        print_fukui_table(&mole, &vals, &format!("π福井函数与双描述符分析（POCV，{ctype}）"));
                        print_pi_dirs(&dirs, &mole);
                    }
                    "mocv" => {
                        let (_stms, _atos, vals) = caler
                            .fukui_pi_mocv(&mole_n, &mole_p, ctype)
                            .expect("π福井函数（MOCV）计算失败");
                        print_fukui_table(&mole, &vals, &format!("π福井函数与双描述符分析（MOCV，{ctype}）"));
                    }
                    _ => unreachable!(),
                }
            }
            AtompropCmds::Freev { filepath, otype } => {
                let mole = load_mole(filepath);
                let caler = activity::Calculator::new(&mole);
                let vals = caler.freev(otype);
                print_freev_table(&mole, &vals, &format!("自由价分析（{otype}）"));
            }
            AtompropCmds::FreevPi { filepath, otype } => {
                let mole = load_mole(filepath);
                let caler = activity::Calculator::new(&mole);
                let (_dirs, vals) = caler.freev_pi(otype).expect("π自由价计算失败");
                print_freev_table(&mole, &vals, &format!("π自由价分析（{otype}）"));
            }
        },

        // ==================== Bondprop ====================
        Commands::Bondprop(cmd) => match cmd {
            BondpropCmds::Order(cmd) => match cmd {
                OrderCmds::Mayer { filepath } => {
                    let mole = load_mole(filepath);
                    let caler = bondprop::order::Calculator::new(&mole);
                    let omat = caler.mayer();
                    print_bond_order(&mole, &omat, "Mayer 键级分析");
                }
                OrderCmds::Wiberg { filepath } => {
                    let mole = load_mole(filepath);
                    let caler = bondprop::order::Calculator::new(&mole);
                    let omat = caler.wiberg();
                    print_bond_order(&mole, &omat, "Wiberg 键级分析");
                }
                OrderCmds::Obt { obt, filepath, otype } => {
                    let mole = load_mole(filepath);
                    let nobt = mole.coefs().obt_occs.len();
                    if *obt == 0 || *obt > nobt {
                        panic!("轨道序号越界: {obt}（共 {nobt} 个轨道，序号从 1 开始）");
                    }
                    let caler = bondprop::order::Calculator::new(&mole);
                    let omat = caler.order_obt(otype, *obt - 1);
                    print_bond_order(&mole, &omat, &format!("分子轨道 {obt} 对{otype}键级的贡献"));
                }
            },
            BondpropCmds::PiOrder {
                filepath,
                method,
                otype,
            } => {
                let mole = load_mole(filepath);
                let caler = bondprop::order::Calculator::new(&mole);
                match method.as_str() {
                    "pocv" => {
                        let (_dirs, omat) = caler.pi_pocv(otype).expect("π键级（POCV）计算失败");
                        print_bond_order(&mole, &omat, &format!("π键级分析（POCV，{otype}）"));
                    }
                    "mocv" => {
                        let (_stms, _atos, omat) = caler.pi_mocv(otype);
                        print_bond_order(&mole, &omat, &format!("π键级分析（MOCV，{otype}）"));
                    }
                    _ => unreachable!(),
                }
            }
            BondpropCmds::Decom { atm1, atm2, filepath } => {
                let mole = load_mole(filepath);
                let natm = mole.atoms().len();
                if *atm1 == 0 || *atm1 > natm || *atm2 == 0 || *atm2 > natm {
                    panic!("原子序号越界: {atm1}, {atm2}（分子共有 {natm} 个原子）");
                }
                if atm1 == atm2 {
                    panic!("键级分解需要两个不同的原子");
                }
                let caler = bondprop::order::Calculator::new(&mole);
                let vals = caler.decom(*atm1 - 1, *atm2 - 1, "mayer");
                let syms = mole.atoms().syms();
                println!(
                    "\nMayer 键级分解: 原子{}({}) - 原子{}({})",
                    atm1,
                    syms[*atm1 - 1],
                    atm2,
                    syms[*atm2 - 1]
                );
                println!("{}", "=".repeat(40));
                for (i, name) in ["σ键级", "π键级", "δ键级", "φ键级"].iter().enumerate() {
                    println!("{:>8}: {:>14.6}", name, vals[i]);
                }
            }
        },

        // ==================== Moleprop ====================
        Commands::Moleprop(cmd) => match cmd {
            MolepropCmds::Aromacity { filepath, rings } => {
                let mole = load_mole(filepath);
                let natm = mole.atoms().len();
                let caler = moleprop::aromacity::Calculator::new(&mole);
                // 未指定环时自动检测分子中的最小环（SSSR）
                let rings: Vec<Vec<usize>> = if rings.is_empty() {
                    mole.rings()
                } else {
                    rings.iter().map(|r| parse_atom_list(r, natm)).collect()
                };
                if rings.is_empty() {
                    // 无环结构：回退为计算全分子非氢键的π键级和
                    let vals = caler.PISD(None).expect("PISD 计算失败");
                    println!("\n芳香性分析（π键级方法，未检测到环，按全分子统计）");
                    println!("{}", "=".repeat(50));
                    println!("全分子π键级和: {:>14.6}", vals[0]);
                    return;
                }
                let vals = caler.PISD(Some(&rings)).expect("PISD 计算失败");
                let syms = mole.atoms().syms();
                println!("\n芳香性分析（π键级方法）");
                println!("{}", "=".repeat(50));
                println!("{:>4}  {:>36}  {:>14}", "环", "原子列表", "PISD 值");
                println!("{:-<60}", "");
                for (i, (ring, &val)) in rings.iter().zip(vals.iter()).enumerate() {
                    let desc: Vec<String> = ring.iter().map(|&a| format!("{}({})", a + 1, syms[a])).collect();
                    println!("{:>4}  {:>36}  {:>14.6}", i + 1, desc.join(" "), val);
                }
            }
        },

        // ==================== Orbtprop ====================
        Commands::Orbtprop(cmd) => match cmd {
            OrbtpropCmds::PiElesPocv { filepath } => {
                let mole = load_mole(filepath);
                let caler = orbtprop::obtval::Calculator::new(&mole);
                let pi_eles = caler.pi_eles_pocv();
                println!("\n分子轨道 π 电子数（POCV 方法）");
                println!("{}", "=".repeat(40));
                for (i, &val) in pi_eles.iter().enumerate() {
                    println!("轨道 {:>4}: {:>12.8}", i + 1, val);
                }
            }
            OrbtpropCmds::PiElesMocv { filepath } => {
                let mole = load_mole(filepath);
                let caler = orbtprop::obtval::Calculator::new(&mole);
                let pi_eles = caler.pi_eles_mocv();
                println!("\n分子轨道 π 电子数（MOCV 方法）");
                println!("{}", "=".repeat(40));
                for (i, &val) in pi_eles.iter().enumerate() {
                    println!("轨道 {:>4}: {:>12.8}", i + 1, val);
                }
            }
            OrbtpropCmds::Composition { obt, filepath, ctype } => {
                let mole = load_mole(filepath);
                let nobt = mole.coefs().obt_occs.len();
                if *obt == 0 || *obt > nobt {
                    panic!("轨道序号越界: {obt}（共 {nobt} 个轨道，序号从 1 开始）");
                }
                let caler = atomprop::charge::Calculator::new(&mole);
                let pops = caler.eles_obt(ctype, *obt - 1);
                print_atom_population(&mole, &pops, &format!("第 {obt} 个分子轨道的 {ctype} 电子布居"));
            }
        },

        // ==================== Gridprop ====================
        Commands::Gridprop(cmd) => match cmd {
            GridpropCmds::Density {
                filepath,
                output,
                step,
                func,
                mole_n,
                mole_p,
            } => {
                let mole = load_mole(filepath);
                // 加载 fukui_* 函数所需的 N+1 / N-1 分子文件
                let mole_n = match mole_n {
                    Some(path) => Some(load_mole(path)),
                    None => None,
                };
                let mole_p = match mole_p {
                    Some(path) => Some(load_mole(path)),
                    None => None,
                };
                let output = match output {
                    Some(path) => path.to_str().expect("输出路径转换失败").to_owned(),
                    None => default_cub_path(filepath, func),
                };
                let title = format!("{} ({})", density_title(func), func);
                compute_and_write_cube(
                    &mole,
                    &output,
                    *step,
                    &title,
                    density_field(&mole, func, mole_n.as_ref(), mole_p.as_ref()),
                );
            }
            GridpropCmds::Potential { filepath, output, step } => {
                let mole = load_mole(filepath);
                let caler = potential::Calculator::new(&mole);
                let output = match output {
                    Some(path) => path.to_str().expect("输出路径转换失败").to_owned(),
                    None => default_cub_path(filepath, "potential"),
                };
                compute_and_write_cube(&mole, &output, *step, "Electrostatic Potential", |grid| {
                    let nuc = caler.nuc_potential(grid);
                    let ele = caler.ele_potential(grid);
                    ele - nuc
                });
            }
            GridpropCmds::Topology { func, filepath } => {
                let mole = load_mole(filepath);
                let caler = topolog::Calculator::new(&mole);
                let cpts = caler.analysis(func);
                println!("\n拓扑分析结果（函数: {func}）");
                println!("{:=<90}", "");
                println!(
                    "{:>4}  {:>20}  {:>20}  {:>20}  {:>12}  {:>8}",
                    "序号", "X (Bohr)", "Y (Bohr)", "Z (Bohr)", "函数值", "临界点类型"
                );
                println!("{:-<90}", "");
                for (i, (pos, val, cpt)) in cpts.iter().enumerate() {
                    let cpt_name = match cpt {
                        -3 => "(3, -3) 核临界点",
                        -1 => "(3, -1) 键临界点",
                        1 => "(3, +1) 环临界点",
                        3 => "(3, +3) 笼临界点",
                        _ => "其他",
                    };
                    println!(
                        "{:>4}  {:>20.12}  {:>20.12}  {:>20.12}  {:>12.8}  {:>8}",
                        i + 1,
                        pos[0],
                        pos[1],
                        pos[2],
                        val,
                        cpt_name,
                    );
                }
                println!("{:=<90}", "");
            }
        },

        // ==================== 提取分子坐标 ====================
        Commands::Coord { input, output, unit } => {
            let path = absolute(input).expect("无法获取文件绝对路径");
            let path_str = path.to_str().expect("路径转字符串失败");
            let mole = Mole::from_file(path_str);
            let atoms = mole.atoms();
            let natm = atoms.len();
            // 程序内部坐标为 Bohr，打印时可选择埃（ang）或玻尔（bohr）
            let xyzs: Vec<[f64; 3]> = match unit.as_str() {
                "bohr" => atoms.xyzs().clone(),
                _ => atoms.xyzs_ang(),
            };
            let unit_name = if unit == "bohr" { "Bohr" } else { "Å" };
            println!("\n分子坐标提取: {}", path_str);
            println!("{}", "=".repeat(50));
            println!(
                "分子式: {}    电荷: {}    多重度: {}    原子数: {}",
                mole.formula(),
                mole.charge(),
                mole.multip(),
                natm
            );
            println!("{:-<70}", "");
            println!(
                "{:>4}  {:>6}  {:>20}  {:>20}  {:>20}",
                "序号",
                "符号",
                format!("X ({})", unit_name),
                format!("Y ({})", unit_name),
                format!("Z ({})", unit_name)
            );
            println!("{:-<70}", "");
            let syms = atoms.syms();
            for i in 0..natm {
                let [x, y, z] = xyzs[i];
                println!("{:>4}  {:>6}  {:>20.10}  {:>20.10}  {:>20.10}", i + 1, syms[i], x, y, z);
            }
            // 可选：导出到 .xyz 文件（xyz 格式固定使用埃）
            if let Some(output) = output {
                let dest = absolute(output).expect("输出文件路径转换失败");
                let dest_path = dest.to_str().expect("路径转字符串失败");
                let ftype = utils::file_type(dest_path);
                if ftype != "xyz" {
                    panic!("坐标导出目前仅支持 .xyz 文件格式，收到: .{ftype}");
                }
                let mut writer = XyzWriter::new();
                writer.title = format!("{} extracted by rswfn", mole.formula());
                writer.read_mole(&mole);
                writer.save(dest_path);
                println!("{:-<70}", "");
                println!("坐标已导出到: {}", dest_path);
            }
        }

        // ==================== 原有命令 ====================
        Commands::Trans { input, output } => {
            let from = absolute(input).expect("输入文件转换失败");
            let dest = absolute(output).expect("输出文件转换失败");
            let from_path = from.to_str().expect("准换字符串失败");
            let dest_path = dest.to_str().expect("准换字符串失败");
            let from_type = utils::file_type(from_path);
            let dest_type = utils::file_type(dest_path);
            let dest_name = utils::file_name(dest_path);
            println!("{} -> {}", from_type, dest_type);
            match (from_type.as_str(), dest_type.as_str()) {
                ("log" | "out", "gjf") => {
                    let reader = LogReader::new(&from_path);
                    let mole = Mole::from_read(&reader);
                    let job = reader.read_job();
                    let mut writer = GjfWriter::new();
                    writer.read_mole(&mole);
                    writer.job = job.to_owned();
                    writer.chk = dest_name;
                    writer.charge = reader.charge as isize;
                    writer.multip = reader.multip as usize;
                    writer.save(&dest_path);
                    println!("文件转换成功 >_< ");
                }
                ("gjf", "dmi") => {
                    let reader = GjfReader::new(&from_path);
                    let mole = Mole::from_read(&reader);
                    let mut writer = DmoWriter::blank();
                    writer.read_mole(&mole);
                    writer.save(&dest_path);
                    println!("文件转换成功 >_< ");
                }
                _ => panic!("不支持的转换类型!!"),
            }
        }
        Commands::Dmo { input, output } => {
            let from = absolute(input)
                .expect("输入文件转换失败")
                .to_str()
                .expect("转换字符串失败")
                .to_owned();
            let dest = match output {
                Some(output) => {
                    let dest = absolute(output).expect("输出文件转换失败");
                    let dest = dest.to_str().expect("转换字符串失败").to_owned();
                    Some(dest)
                }
                None => None,
            };
            qmcalc::dhmo::output(&from, dest);
            println!("DHMO计算完成 >_< ");
        }
        // ==================== 提取支持信息 ====================
        Commands::Si { inputs, output } => {
            if inputs.is_empty() {
                panic!("至少需要指定一个输入文件");
            }
            let paths: Vec<String> = inputs
                .iter()
                .map(|path| {
                    absolute(path)
                        .expect("无法获取文件绝对路径")
                        .to_str()
                        .expect("路径转字符串失败")
                        .to_owned()
                })
                .collect();
            let tool = Tool::new();
            let text = tool.get_SI(&paths);
            match output {
                Some(path) => {
                    let dest = absolute(path).expect("输出文件路径转换失败");
                    let dest_path = dest.to_str().expect("路径转字符串失败");
                    std::fs::write(dest_path, text).expect("写入支持信息文件失败");
                    println!("支持信息已导出到: {}", dest_path);
                }
                None => print!("{}", text),
            }
        }

        Commands::Dev => {
            dev::test();
        }
    }
}
