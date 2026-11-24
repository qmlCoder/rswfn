mod base;
mod datas;
mod editor;
mod loger;
mod march;
mod maths;
mod qmcal;
mod reader;
mod rsbse;
mod rsint;
mod rsxcs;
mod tools;
mod utils;
mod writer;

mod atomprop;
mod bondprop;
mod gridprop;
mod moleprop;
mod orbtprop;

mod cli;
mod dev;
mod tui;
use std;

/// 程序入口，根据命令行参数决定启动 CLI 模式或 TUI 模式
fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() > 1 {
        cli::main();
    } else {
        tui::Main::new();
    }
}
