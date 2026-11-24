#![allow(non_snake_case)]
mod cub;
mod dmo;
mod fch;
mod gjf;
mod inp;
mod log;
mod mde;
mod mfn;
mod mol;
mod non;
mod out;
mod sdf;
mod xyz;

use crate::base::Atoms;
use crate::base::Basis;
use crate::base::Coefs;
use crate::base::Mole;
use crate::datas::btrans::sym2stp;
use crate::maths;

pub use cub::CubReader;
pub use dmo::DmoReader;
pub use fch::FchReader;
pub use gjf::{GauJob, GjfReader, parse_job};
pub use inp::InpReader;
pub use log::LogReader;
pub use mde::MdeReader;
pub use mfn::MfnReader;
pub use mol::MolReader;
pub use non::NonReader;
pub use out::OutReader;
pub use sdf::SdfReader;
pub use xyz::XyzReader;

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

/// 用于在文件中标记关键数据位置的索引集合
pub struct Marks {
    /// 标记名称到行号的映射
    pub marks: HashMap<String, usize>,
}

impl Marks {
    /// 创建新的Marks实例
    pub fn new() -> Marks {
        Self { marks: HashMap::new() }
    }

    /// 添加一个标记
    pub fn add(&mut self, marker: String, number: usize) {
        self.marks.insert(marker, number);
    }
}

/// 文件读取器trait
///
/// 定义从各种量子化学程序输出文件中读取分子信息的统一接口。
pub trait Reader {
    /// 获取分子中的原子信息
    fn get_atoms(&self) -> Atoms;
    /// 获取基组信息(如果文件包含)
    fn get_basis(&self) -> Option<Basis>;
    /// 获取波函数系数(如果文件包含)
    fn get_coefs(&self) -> Option<Coefs>;
}

/// 读取文件的全部行
///
/// 打开指定路径的文件,逐行读取并返回一个字符串向量。
pub fn read_lines(path: &str) -> Vec<String> {
    let file = File::open(path).expect(format!("无法打开文件:{}\n", path).as_ref());
    let read = BufReader::new(file);
    let mut texts = Vec::new();
    for line in read.lines() {
        let line = line.expect("line read error");
        texts.push(line);
    }
    texts
}

/// 根据电荷和自旋多重度计算α和β电子数
pub fn charge_number(charge: isize, multip: usize, total: usize) -> [usize; 2] {
    let N = total as isize;
    let C = charge as isize;
    let S = multip as isize;
    let nela = (N - C + S - 1) / 2;
    let nelb = (N - C - S + 1) / 2;
    [nela as usize, nelb as usize]
}

#[derive(Debug)]
pub enum FileType {
    Gjf, // 高斯输入文件
    Log, // 高斯输出文件
    Fch, // 高斯检查点文件
    Mol, // 分子文件
    Sdf,
    Xyz,
    Mde, // molden文件
    Mfn, // Multiwfn波函数文件
    Dmo, // 自定义的包含原子局部坐标系的文件
    Inp, // Orca输入文件
    Cub, // 格点文件
}

/// 判断文件类型
///
/// 根据文件扩展名返回对应的`FileType`枚举值。
pub fn get_fileType(path: &str) -> FileType {
    if path.ends_with(".gjf") {
        return FileType::Gjf;
    }
    if path.ends_with(".log") {
        return FileType::Log;
    }
    if path.ends_with(".out") {
        return FileType::Log;
    }
    if path.ends_with(".fch") || path.ends_with(".fchk") {
        return FileType::Fch;
    }
    if path.ends_with(".mol") {
        return FileType::Mol;
    }
    if path.ends_with(".sdf") {
        return FileType::Sdf;
    }
    if path.ends_with(".xyz") {
        return FileType::Xyz;
    }
    if path.ends_with(".molden") {
        return FileType::Mde;
    }
    if path.ends_with(".mwfn") {
        return FileType::Mfn;
    }
    if path.ends_with(".dmi") || path.ends_with(".dmo") {
        return FileType::Dmo;
    }
    if path.ends_with(".inp") {
        return FileType::Inp;
    }
    if path.ends_with(".cub") {
        return FileType::Cub;
    }
    panic!("未知文件类型:{path}");
}

use crate::datas::btrans;

/// 确定每个壳层的基函数形式(球谐或笛卡尔)
pub fn shl_forms(ato_atms: &Vec<usize>, ato_shls: &Vec<usize>, ato_syms: &Vec<String>) -> Vec<String> {
    let mut shl_form = vec![];
    let nato = ato_atms.len();
    let mut up = 0;
    let mut lo = 0;
    let mut keys = vec![];
    for i in 0..nato {
        let atm = ato_atms[i];
        let shl = ato_shls[i];
        let ang = btrans::sym2ang(&ato_syms[i]);
        let key = (atm, shl, ang);
        if !keys.contains(&key) {
            keys.push(key);
            let stp = btrans::sym2stp(&ato_syms[i]);
            let form = match stp {
                0 | 1 => "sph",
                2 => "car",
                _ => panic!("未知的壳层类型"),
            };
            shl_form.push(form.to_owned());
        }
    }
    shl_form
}

use crate::base::Shl;
pub type BasData = (usize, usize, usize, f64, f64); // atm, shl, ang, alp, coe
/// 从原始基组数据构建Basis对象
pub fn shl_basis(datas: &Vec<BasData>, forms: &Vec<String>) -> Basis {
    let mut keys = vec![];
    let mut shls = Vec::new();
    let mut ishl = 0;
    for data in datas.iter() {
        let (atm, shl, ang, alp, coe) = *data;
        let key = (atm, shl, ang);
        if !keys.contains(&key) {
            keys.push(key);
            let form = forms[ishl].clone();
            shls.push(Shl {
                atm,
                shl,
                ang,
                form,
                alps: vec![],
                coes: vec![],
            });
            ishl += 1;
        }
        let nshl = shls.len();
        shls[nshl - 1].alps.push(alp);
        shls[nshl - 1].coes.push(coe);
    }
    Basis::new(shls, None)
}

/// 振动频率信息
///
/// 包含频率、约化质量、力常数、红外强度及振动模式向量。
pub struct Vibe {
    pub freq: f64,
    /// 约化质量
    pub redm: f64,
    /// 力常数
    pub frcc: f64,
    /// 红外强度
    pub irin: f64,
    /// 原子振动位移向量列表
    pub vibs: Vec<[f64; 3]>,
}

impl Vibe {
    /// 获取振动幅度最大的原子及其位移方向
    pub fn get_max(&self) -> (usize, [f64; 3], f64) {
        let mut max_idx = 0;
        let mut max_vib = [0.0; 3];
        let mut max_len = 0.0;
        for (i, vib) in self.vibs.iter().enumerate() {
            let li = maths::vec_len(vib);
            if li > max_len {
                max_len = li;
                max_idx = i;
                max_vib = vib.clone();
            }
        }
        (max_idx, max_vib, max_len)
    }
}

#[derive(Debug)]
pub struct Nics {
    pub iso: f64,
    pub ani: f64,
    pub mat: [f64; 9],
    pub eig: [f64; 3],
}
