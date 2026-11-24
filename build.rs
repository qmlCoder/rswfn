//! 编译期扫描 `src/rsbse/txts` 下的基组文件, 生成内嵌到程序中的基组表.
//!
//! 往该目录新增 `.txt` 文件后重新编译即可通过 `rsbse::from_name` 使用,
//! 基组名即文件名(不含扩展名), 匹配时忽略大小写和空白字符.

use std::env;
use std::fmt::Write as _;
use std::fs;
use std::path::Path;

fn main() {
    let manifest_dir = env::var("CARGO_MANIFEST_DIR").expect("读取CARGO_MANIFEST_DIR失败");
    let txts_dir = Path::new(&manifest_dir).join("src").join("rsbse").join("txts");
    // 目录本身用于感知新增/删除文件, 下面的每个文件用于感知内容修改
    println!("cargo:rerun-if-changed={}", txts_dir.display());

    let mut bases = Vec::new();
    if txts_dir.is_dir() {
        for entry in fs::read_dir(&txts_dir)
            .expect("读取基组目录失败")
            .filter_map(|entry| entry.ok())
        {
            let path = entry.path();
            if path.extension().and_then(|ext| ext.to_str()) != Some("txt") {
                continue;
            }
            let (Some(stem), Some(file)) = (
                path.file_stem().and_then(|stem| stem.to_str()),
                path.file_name().and_then(|file| file.to_str()),
            ) else {
                continue;
            };
            let name = stem
                .chars()
                .filter(|char| !char.is_whitespace())
                .collect::<String>()
                .to_lowercase();
            println!("cargo:rerun-if-changed={}", path.display());
            bases.push((name, file.to_string()));
        }
    } else {
        println!("cargo:warning=基组目录不存在: {}", txts_dir.display());
    }
    // 排序保证生成结果稳定
    bases.sort();

    let mut code = String::new();
    code.push_str("// 该文件由 build.rs 自动生成, 请勿手动修改\n");
    code.push_str("pub static BASIS_SETS: &[(&str, &str)] = &[\n");
    for (name, file) in bases.iter() {
        let _ = writeln!(
            code,
            "    ({name:?}, include_str!(concat!(env!(\"CARGO_MANIFEST_DIR\"), \"/src/rsbse/txts/\", {file:?}))),"
        );
    }
    code.push_str("];\n");

    let out_dir = env::var("OUT_DIR").expect("读取OUT_DIR失败");
    let out_path = Path::new(&out_dir).join("basis_sets.rs");
    fs::write(&out_path, code).expect("写入basis_sets.rs失败");
}
