use std::io::Write;

use simplelog;

/// 仅写 stderr 的轻量日志器。
///
/// CLI 用它输出进度（SCF 迭代等），stdout 留给分析结果/JSON。
struct StderrLogger;

static STDERR_LOGGER: StderrLogger = StderrLogger;

impl log::Log for StderrLogger {
    fn enabled(&self, metadata: &log::Metadata) -> bool {
        metadata.level() <= log::Level::Info
    }

    fn log(&self, record: &log::Record) {
        if self.enabled(record.metadata()) {
            let _ = writeln!(std::io::stderr(), "[{}] {}", record.level(), record.args());
        }
    }

    fn flush(&self) {
        let _ = std::io::stderr().flush();
    }
}

/// 初始化 CLI 日志：只写 stderr，不污染 stdout。
pub fn init_stderr_log() {
    let _ = log::set_logger(&STDERR_LOGGER);
    log::set_max_level(log::LevelFilter::Info);
}

#[rustfmt::skip]
/// 初始化调试日志系统，同时输出到终端和日志文件（仅 dev 测试使用）。
pub fn init_log() {
    let file = std::fs::File::create("rswfn日志.txt").ok();
    let mut loggers: Vec<Box<dyn simplelog::SharedLogger>> = vec![
        simplelog::TermLogger::new(
            simplelog::LevelFilter::Info,
            simplelog::Config::default(),
            simplelog::TerminalMode::Mixed,
            simplelog::ColorChoice::Auto,
        ),
    ];
    if let Some(file) = file {
        loggers.push(simplelog::WriteLogger::new(
            simplelog::LevelFilter::Info,
            simplelog::Config::default(),
            file,
        ));
    }
    if simplelog::CombinedLogger::init(loggers).is_ok() {
        log::info!("日志系统已启动");
    }
}
