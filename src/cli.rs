use std::path::PathBuf;

use clap::{Parser, Subcommand, command};

use crate::dev;
use crate::report::Report;
use crate::service::{self, atomprop, bondprop, file, gridprop, moleprop, orbtprop, qmcal};
use crate::viewer;

#[derive(Parser, Debug)]
#[command(author = "Xiaofei Shi", version, about = "Weave Function Analysis Tool")]
struct Cli {
    /// 以 JSON 输出结构化结果（便于脚本处理）
    #[arg(long, global = true)]
    json: bool,
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

    /// 量子化学计算（根据坐标与基组求解波函数）
    #[command(subcommand)]
    Qmcal(QmcalCmds),

    /// 提取支持信息（Supporting Information），汇总能量、频率与结构
    Si {
        /// 输入 Gaussian 输出文件路径（.log/.out），可指定多个
        inputs: Vec<PathBuf>,
        /// 输出文本文件路径（可选，默认打印到屏幕）
        #[arg(short = 'o', long = "output")]
        output: Option<PathBuf>,
    },

    /// 在终端预览分子结构（球棍 / 空间填充 / 线框）
    View {
        /// 输入文件路径（支持 xyz/mol/sdf/gjf/log/out/fch/fchk/mwfn/dmo 等格式）
        input: PathBuf,
        /// 画布尺寸，如 100x40；默认取当前终端大小
        #[arg(long, value_parser = parse_size)]
        size: Option<(u16, u16)>,
        /// 绘制模型: ball（球棍） | space（空间填充） | wire（线框）
        #[arg(long, default_value = "ball", value_parser = ["ball", "space", "wire"])]
        style: String,
        /// 初始视角（度）: yaw,pitch，如 20,12
        #[arg(long, default_value = "20,12", value_parser = parse_rotate)]
        rotate: (f64, f64),
        /// 缩放倍数（1.0 为自动铺满画布）
        #[arg(long, default_value_t = 1.0)]
        zoom: f64,
        /// 不显示元素符号
        #[arg(long)]
        no_labels: bool,
        /// 关闭景深雾化
        #[arg(long)]
        no_depth: bool,
        /// 颜色输出: auto（仅终端上色） | always | never
        #[arg(long, default_value = "auto", value_parser = ["auto", "always", "never"])]
        color: String,
    },

    /// 开发时测试
    Dev,
}

/// 解析 `100x40` 形式的画布尺寸
fn parse_size(value: &str) -> Result<(u16, u16), String> {
    let (w, h) = value.split_once(['x', 'X']).ok_or_else(|| format!("尺寸格式应为 宽x高，例如 100x40：{value}"))?;
    let w: u16 = w.trim().parse().map_err(|_| format!("宽度不是有效数字: {w}"))?;
    let h: u16 = h.trim().parse().map_err(|_| format!("高度不是有效数字: {h}"))?;
    if w < 4 || h < 2 {
        return Err("画布至少需要 4 列 × 2 行".to_string());
    }
    Ok((w, h))
}

/// 解析 `34,-18` 形式的初始视角（度）
fn parse_rotate(value: &str) -> Result<(f64, f64), String> {
    let (yaw, pitch) = value.split_once(',').ok_or_else(|| format!("角度格式应为 yaw,pitch，例如 34,-18：{value}"))?;
    let yaw: f64 = yaw.trim().parse().map_err(|_| format!("yaw 不是有效数字: {yaw}"))?;
    let pitch: f64 = pitch.trim().parse().map_err(|_| format!("pitch 不是有效数字: {pitch}"))?;
    Ok((yaw, pitch))
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
    /// RESP（限制性静电势）电荷分析
    Resp {
        /// 波函数文件路径（支持 .log/.out/.fch/.fchk）
        filepath: PathBuf,
        /// MK 拟合点密度（点/Å²）
        #[arg(long, default_value_t = 6.0)]
        pt_density: f64,
        /// 只执行第一阶段（弱限制）拟合
        #[arg(long, default_value_t = false)]
        one_stage: bool,
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
    /// 键能量分解（bond / rest / cross；仅支持闭壳层杂化泛函的 log/out 文件）
    Energy {
        /// 原子1 序号（从1开始）
        atm1: usize,
        /// 原子2 序号（从1开始）
        atm2: usize,
        /// 高斯 log/out 文件路径（需含泛函与基组信息）
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
    /// π 电子能量分解（POCV/MOCV；仅支持闭壳层杂化泛函的 log/out 文件）
    PiEnergy {
        /// 高斯 log/out 文件路径（需含泛函与基组信息）
        filepath: PathBuf,
        /// 轨道改造方法: pocv | mocv
        #[arg(long, default_value = "pocv", value_parser = ["pocv", "mocv"])]
        method: String,
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

// ===== qmcal 子命令 =====

#[derive(Subcommand, Debug)]
enum QmcalCmds {
    /// 闭壳层 Hartree-Fock（RHF）计算
    Hf {
        /// 输入文件路径（需包含坐标；基组可由文件提供，或用 -b 指定）
        #[arg(short, long)]
        input: PathBuf,
        /// 基组名（使用程序提供的基组数据，如 6-31G(d,p)、def2-SVP），指定后忽略文件中的基组
        #[arg(short = 'b', long = "basis")]
        basis: Option<String>,
        /// 输出 .fch 波函数文件路径（可选）
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// 双电子积分 Schwarz 筛选阈值（越大越快、精度越低）
        #[arg(long = "screen", default_value = "1e-11")]
        screen: f64,
        /// 双电子积分块内存上限（MB），超过后自动退化为直接 SCF
        #[arg(long = "max-mem", default_value = "8000")]
        max_mem: usize,
        /// SCF 能量收敛阈值（Hartree）
        #[arg(long = "conv", default_value = "1e-8")]
        conv: f64,
        /// 最大 SCF 迭代次数
        #[arg(long = "maxiter", default_value = "100")]
        maxiter: usize,
    },
    /// 密度泛函计算（默认 B3LYP，局部相关用 VWN3，与 Gaussian 一致）
    ///
    /// 闭壳层体系用 RKS，开壳层体系（α ≠ β 电子数）自动切换到 UKS。
    Dft {
        /// 输入文件路径（需包含坐标；基组可由文件提供，或用 -b 指定）
        #[arg(short, long)]
        input: PathBuf,
        /// 基组名（使用程序提供的基组数据，如 6-31G(d,p)、def2-SVP），指定后忽略文件中的基组
        #[arg(short = 'b', long = "basis")]
        basis: Option<String>,
        /// 交换相关泛函：内置名（b3lyp/pbe0/pbe/blyp/bp86/svwn5，lda 为 svwn5 别名）
        /// 或外部 rsxcs 支持的 libxc 短名（如 gga_x_pw91，需 use-rsxcs 特性）
        #[arg(short = 'f', long = "func", default_value = "b3lyp")]
        func: String,
        /// 输出 .fch 波函数文件路径（可选）
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// 径向格点数（默认 99，对应 Gaussian 的 UltraFine）
        #[arg(long = "radial", default_value = "99")]
        radial: usize,
        /// 球面（Lebedev）格点数（默认 590，对应 Gaussian 的 UltraFine）
        #[arg(long = "angular", default_value = "590")]
        angular: usize,
        /// 关闭径向-角向剪枝（默认开启：格点数约为不剪枝的 1/4）
        #[arg(long = "no-prune", default_value = "false")]
        no_prune: bool,
        /// 格点积分中显著原子轨道的距离截断（Bohr；默认 12，设为 0 关闭）
        #[arg(long = "ao-cut", default_value = "12")]
        ao_cut: f64,
        /// 双电子积分 Schwarz 筛选阈值（越大越快、精度越低）
        #[arg(long = "screen", default_value = "1e-11")]
        screen: f64,
        /// 双电子积分块内存上限（MB），超过后自动退化为直接 SCF
        #[arg(long = "max-mem", default_value = "8000")]
        max_mem: usize,
        /// SCF 能量收敛阈值（Hartree）
        #[arg(long = "conv", default_value = "1e-8")]
        conv: f64,
        /// 最大 SCF 迭代次数
        #[arg(long = "maxiter", default_value = "60")]
        maxiter: usize,
    },
    /// 局部坐标系下的 π-only 密度泛函计算（MOCV，π 约化基，矩阵显著变小）
    ///
    /// 仅保留局部坐标系中 π 对称的原子轨道，得到的能量是 π 子系统模型值，
    /// 不是分子总能量。
    Lks {
        /// 输入文件路径（需包含坐标；基组可由文件提供，或用 -b 指定）
        #[arg(short, long)]
        input: PathBuf,
        /// 基组名（使用程序提供的基组数据，如 6-31G(d,p)、def2-SVP），指定后忽略文件中的基组
        #[arg(short = 'b', long = "basis")]
        basis: Option<String>,
        /// 交换相关泛函：内置名（b3lyp/pbe0/pbe/blyp/bp86/svwn5，lda 为 svwn5 别名）
        #[arg(short = 'f', long = "func", default_value = "b3lyp")]
        func: String,
        /// π 电子数（缺省时按参与共轭的 p 轨道数量估计）
        #[arg(long = "npi-elec")]
        npi_elec: Option<usize>,
        /// 输出波函数文件路径（.fch/.fchk 或 .mwfn；系数已展开回全局基）
        #[arg(short, long)]
        output: Option<PathBuf>,
        /// 打印局部坐标系下的 π 分子轨道系数表
        #[arg(long = "print-mo")]
        print_mo: bool,
        /// 双电子积分 Schwarz 筛选阈值（越大越快、精度越低）
        #[arg(long = "screen", default_value = "1e-11")]
        screen: f64,
        /// 双电子积分块内存上限（MB），超过后自动退化为直接模式
        #[arg(long = "max-mem", default_value = "8000")]
        max_mem: usize,
        /// SCF 能量收敛阈值（Hartree）
        #[arg(long = "conv", default_value = "1e-8")]
        conv: f64,
        /// 最大 SCF 迭代次数
        #[arg(long = "maxiter", default_value = "60")]
        maxiter: usize,
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
    /// 静电势（原子单位，核项 + 电子项），输出 cub 文件
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

/// CLI 模式入口函数，解析命令行参数执行对应操作
pub fn main() {
    crate::loger::init_stderr_log();
    let cli = Cli::parse();
    // 结构预览输出的是终端字符画而不是分析报告，因此不经过 Report
    if let Commands::View {
        input,
        size,
        style,
        rotate,
        zoom,
        no_labels,
        no_depth,
        color,
    } = &cli.command
    {
        let opts = viewer::RenderOptions {
            size: *size,
            style: viewer::Style::parse(style).unwrap_or_default(),
            rotate: *rotate,
            zoom: *zoom,
            labels: !*no_labels,
            depth_cue: !*no_depth,
            color: viewer::ColorMode::parse(color).unwrap_or_default(),
            json: cli.json,
        };
        if let Err(err) = viewer::run(input, &opts) {
            eprintln!("错误: {err}");
            std::process::exit(1);
        }
        return;
    }
    match run(&cli.command) {
        Ok(report) => {
            if cli.json {
                println!("{}", report.to_json());
            } else {
                print!("{report}");
            }
        }
        Err(err) => {
            eprintln!("错误: {err}");
            std::process::exit(1);
        }
    }
}

/// 把命令行参数映射为服务层任务，文本/JSON 输出由前端统一处理
fn run(command: &Commands) -> Result<Report, String> {
    match command {
        Commands::Atomprop(cmd) => run_atomprop(cmd),
        Commands::Bondprop(cmd) => run_bondprop(cmd),
        Commands::Moleprop(cmd) => run_moleprop(cmd),
        Commands::Orbtprop(cmd) => run_orbtprop(cmd),
        Commands::Gridprop(cmd) => run_gridprop(cmd),
        Commands::Coord { input, output, unit } => file::coord(input, output.as_deref(), unit),
        Commands::Trans { input, output } => file::trans(input, output),
        Commands::Dmo { input, output } => file::dmo(input, output.as_deref()),
        Commands::Qmcal(cmd) => run_qmcal(cmd),
        Commands::Si { inputs, output } => file::si(inputs, output.as_deref()),
        Commands::View { .. } => Err("view 由 CLI 入口直接处理".to_string()),
        Commands::Dev => {
            dev::test();
            Ok(Report::new("开发测试").note("dev 测试已执行"))
        }
    }
}

fn run_atomprop(cmd: &AtompropCmds) -> Result<Report, String> {
    match cmd {
        AtompropCmds::Charge(cmd) => match cmd {
            ChargeCmds::Mulliken { filepath } => atomprop::mulliken(&service::load_mole(filepath)?),
            ChargeCmds::Lowdin { filepath } => atomprop::lowdin(&service::load_mole(filepath)?),
            ChargeCmds::Hirshfeld { filepath } => atomprop::hirshfeld(&service::load_mole(filepath)?),
            ChargeCmds::Spin { filepath, ctype } => atomprop::spin(&service::load_mole(filepath)?, ctype),
            ChargeCmds::Resp { filepath, pt_density, one_stage } => atomprop::resp(&service::load_mole(filepath)?, *pt_density, *one_stage),
        },
        AtompropCmds::PiPop { filepath, ctype, method } => atomprop::pi_pop(&service::load_mole(filepath)?, ctype, method),
        AtompropCmds::Fukui { filepath, anion, cation, ctype } => {
            let mole = service::load_mole(filepath)?;
            let mole_n = service::load_mole(anion)?;
            let mole_p = service::load_mole(cation)?;
            atomprop::fukui(&mole, &mole_n, &mole_p, ctype)
        }
        AtompropCmds::FukuiPi { filepath, anion, cation, ctype, method } => {
            let mole = service::load_mole(filepath)?;
            let mole_n = service::load_mole(anion)?;
            let mole_p = service::load_mole(cation)?;
            atomprop::fukui_pi(&mole, &mole_n, &mole_p, ctype, method)
        }
        AtompropCmds::Freev { filepath, otype } => atomprop::freev(&service::load_mole(filepath)?, otype),
        AtompropCmds::FreevPi { filepath, otype } => atomprop::freev_pi(&service::load_mole(filepath)?, otype),
    }
}

fn run_bondprop(cmd: &BondpropCmds) -> Result<Report, String> {
    match cmd {
        BondpropCmds::Order(cmd) => match cmd {
            OrderCmds::Mayer { filepath } => bondprop::mayer(&service::load_mole(filepath)?),
            OrderCmds::Wiberg { filepath } => bondprop::wiberg(&service::load_mole(filepath)?),
            OrderCmds::Obt { obt, filepath, otype } => bondprop::order_obt(&service::load_mole(filepath)?, *obt, otype),
        },
        BondpropCmds::PiOrder { filepath, method, otype } => bondprop::pi_order(&service::load_mole(filepath)?, method, otype),
        BondpropCmds::Decom { atm1, atm2, filepath } => bondprop::decom(&service::load_mole(filepath)?, *atm1, *atm2),
        BondpropCmds::Energy { atm1, atm2, filepath } => bondprop::bond_energy(&service::load_mole(filepath)?, *atm1, *atm2),
    }
}

fn run_moleprop(cmd: &MolepropCmds) -> Result<Report, String> {
    match cmd {
        MolepropCmds::Aromacity { filepath, rings } => moleprop::aromacity(&service::load_mole(filepath)?, rings),
        MolepropCmds::PiEnergy { filepath, method } => moleprop::pi_energy(&service::load_mole(filepath)?, method),
    }
}

fn run_orbtprop(cmd: &OrbtpropCmds) -> Result<Report, String> {
    match cmd {
        OrbtpropCmds::PiElesPocv { filepath } => orbtprop::pi_eles_pocv(&service::load_mole(filepath)?),
        OrbtpropCmds::PiElesMocv { filepath } => orbtprop::pi_eles_mocv(&service::load_mole(filepath)?),
        OrbtpropCmds::Composition { obt, filepath, ctype } => orbtprop::composition(&service::load_mole(filepath)?, *obt, ctype),
    }
}

fn run_gridprop(cmd: &GridpropCmds) -> Result<Report, String> {
    match cmd {
        GridpropCmds::Density {
            filepath,
            output,
            step,
            func,
            mole_n,
            mole_p,
        } => {
            let mole = service::load_mole(filepath)?;
            let mole_n = match mole_n {
                Some(path) => Some(service::load_mole(path)?),
                None => None,
            };
            let mole_p = match mole_p {
                Some(path) => Some(service::load_mole(path)?),
                None => None,
            };
            gridprop::density(&mole, output.as_deref(), *step, func, mole_n.as_ref(), mole_p.as_ref())
        }
        GridpropCmds::Potential { filepath, output, step } => gridprop::potential(&service::load_mole(filepath)?, output.as_deref(), *step),
        GridpropCmds::Topology { func, filepath } => gridprop::topology(&service::load_mole(filepath)?, func),
    }
}

fn run_qmcal(cmd: &QmcalCmds) -> Result<Report, String> {
    match cmd {
        QmcalCmds::Hf {
            input,
            basis,
            output,
            screen,
            max_mem,
            conv,
            maxiter,
        } => {
            let opts = service::qmcal::ScfOptions {
                screen: *screen,
                max_mem_mb: *max_mem,
                conv: *conv,
                max_iter: *maxiter,
            };
            qmcal::hf_with(input, basis.as_deref(), output.as_deref(), &opts)
        }
        QmcalCmds::Dft {
            input,
            basis,
            func,
            output,
            radial,
            angular,
            no_prune,
            ao_cut,
            screen,
            max_mem,
            conv,
            maxiter,
        } => {
            let opts = service::qmcal::ScfOptions {
                screen: *screen,
                max_mem_mb: *max_mem,
                conv: *conv,
                max_iter: *maxiter,
            };
            let grid = service::qmcal::GridOptions {
                radial: *radial,
                angular: *angular,
                prune: !*no_prune,
                ao_cut: *ao_cut,
            };
            qmcal::dft_with(input, basis.as_deref(), func, output.as_deref(), &opts, &grid)
        }
        QmcalCmds::Lks {
            input,
            basis,
            func,
            npi_elec,
            output,
            print_mo,
            screen,
            max_mem,
            conv,
            maxiter,
        } => {
            let opts = service::qmcal::ScfOptions {
                screen: *screen,
                max_mem_mb: *max_mem,
                conv: *conv,
                max_iter: *maxiter,
            };
            qmcal::lks_with(input, basis.as_deref(), func, *npi_elec, output.as_deref(), *print_mo, &opts)
        }
    }
}
