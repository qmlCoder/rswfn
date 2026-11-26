mod base;
mod datas;
mod editor;
mod loger;
mod maths;
mod qmcal;
mod reader;
mod rsbse;
mod rsint;
mod rsxcs;
mod tools;
mod utils;
mod writer;

mod props;

mod cli;
mod dev;
mod report;
mod service;
mod tui;
mod viewer;
use std;

/// 程序入口，根据命令行参数决定启动 CLI 模式或 TUI 模式
fn main() {
    let args: Vec<String> = std::env::args().collect();

    if args.len() > 1 {
        // rswfn --tui [路径]
        if args[1] == "--tui" {
            let initial = args.get(2).map(|path| std::path::PathBuf::from(path.as_str()));
            tui::run(initial);
            return;
        }
        // rswfn <已存在的文件或目录> 直接进入 TUI
        if args.len() == 2 && !args[1].starts_with('-') {
            let path = std::path::PathBuf::from(args[1].as_str());
            if path.exists() {
                tui::run(Some(path));
                return;
            }
        }
        cli::main();
    } else {
        tui::run(None);
    }
}
