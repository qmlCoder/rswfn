use std::collections::HashMap;
use std::f64::consts::PI;

use log::Level;
use ndarray::{Array2, s};
use ndarray::{Axis, Shape};
use rayon::iter::IntoParallelIterator;
use rayon::iter::IntoParallelRefIterator;
use rayon::prelude::*;

use crate::base::{Mole, Stm, Vec3};
use crate::bondprop::direction;
use crate::datas::consts::BOHR;
use crate::datas::{self, btrans};
use crate::datas::{consts, roots};
use crate::editor::Editor;
use crate::gridprop;
use crate::gridprop::{density, wfnfunc};
use crate::loger;
use crate::march;
use crate::maths;
use crate::moleprop;
use crate::orbtprop;
use crate::orbtprop::mocv;
use crate::reader::{GjfReader, LogReader, MdeReader, Reader, XyzReader};
use crate::rsint::get_2eid;
use crate::tools;
use crate::utils::{self, comp_matrix};
use crate::writer;
use crate::writer::{FchWriter, GjfWriter, Writer};
use crate::{atomprop, qmcal};
use crate::{bondprop, editor};

/// 测试读取 log 文件并解析分子信息
pub fn test_reader() {
    let path = r"C:\Users\11032\Desktop\gfile\errors\CO-BP86-opt-1.fchk";
    let mole = Mole::from_file(path);
    println!("mole: {:?}", mole);
}

pub fn test_charge() {
    let path = String::from(r"E:\code\rswfn\moles\C6H6.fch");
    let mole = Mole::from_file(&path);
    let caler = atomprop::charge::Calculator::new(&mole);
    let res = caler.hirshfeld();
    println!("res:{:?}", res);
}

pub fn test_topology() {
    let path = String::from(r"E:\code\rswfn\moles\C6H6.fch");
    let mole = Mole::from_file(&path);
    let caler = gridprop::topolog::Calculator::new(&mole);
    let res = caler.analysis("mol_rho_dm");
    for (pos, val, cpt) in res {
        let [px, py, pz] = pos;
        println!(
            "pos: ({:>10.4}, {:>10.4}, {:>10.4}), val: {}, cpt: {}",
            px, py, pz, val, cpt
        );
    }
}

pub fn test_order() {
    let path = r"E:\code\articles\MOCV\moles\bond-order-decom\CO2.log";
    let mole = Mole::from_file(path);

    // ① 几何：确认线性/取向
    for i in 0..mole.atoms().len() {
        let atom = mole.atoms().get(i);
        println!("atom {}: sym={} xyz={:?}", i, atom.sym, atom.xyz);
    }
    let bond_vec = mole.atoms().vec(0, 2);
    println!("bond 0->2 向量: {:?}", bond_vec.into_arr());

    // 生成 72 个绕 z 轴旋转的局部坐标系（endpoint=False，步长 2π/72）
    let n_ang = 72;
    let stml: Vec<Stm> = (0..n_ang)
        .map(|i| {
            let ang = 2.0 * PI * (i as f64) / (n_ang as f64);
            let xx = ang.cos();
            let xy = ang.sin();
            Stm::new([xx, xy, 0.0], [-xy, xx, 0.0], [0.0, 0.0, 1.0])
        })
        .collect();

    // 两组要保留的轨道类型：π(PX) 和 δ(D-2)
    let syml: [&[&str]; 2] = [&["PX"], &["D-2"]];

    let ato_syms = mole.basis().atos().syms("sph");

    // MOCV：data[i, j] 第 i 个朝向、第 j 种轨道类型下的 Mayer 键级 omat[0,2]
    // POCV-z：data_pocv[i] 方向取局部坐标系 x 轴（绕全局 z 轴，应与 MOCV-PX 逐角度相等）
    let mut data = Array2::<f64>::zeros((stml.len(), syml.len()));
    let mut data_pocv = vec![0.0; stml.len()];
    for (i, stm) in stml.iter().enumerate() {
        for (j, syms) in syml.iter().enumerate() {
            let stms: HashMap<usize, Stm> = HashMap::from([(0usize, stm.clone()), (2usize, stm.clone())]);
            let atos: Vec<usize> = ato_syms
                .iter()
                .enumerate()
                .filter(|(_, s)| syms.contains(&s.as_str()))
                .map(|(idx, _)| idx)
                .collect();
            let caler = bondprop::order::Calculator::new(&mole);
            let omat = caler.mocv(&stms, &atos, "mayer");
            data[(i, j)] = omat[(0, 2)];
        }

        let dir = Vec3::from_arr(stm.ex);
        let dirs: HashMap<usize, Vec3> = HashMap::from([(0usize, dir.clone()), (2usize, dir)]);
        let caler = bondprop::order::Calculator::new(&mole);
        let omat = caler.pocv(&dirs, false, false, "mayer").expect("pocv 计算失败");
        data_pocv[i] = omat[(0, 2)];
    }

    // ② POCV-verts：与 wmview bondfunc.rs 完全一致，
    //    绕 0->2 实际键轴转 36 个方向（每 10° 一个），取 val.max(0).sqrt()
    let dir_caler = bondprop::direction::Calculator::new(&mole);
    let mut data_verts: Vec<(usize, f64)> = Vec::new();
    for (i, dir) in dir_caler.verts(0, 2).iter().enumerate() {
        let dirs: HashMap<usize, Vec3> = HashMap::from([(0usize, dir.clone()), (2usize, dir.clone())]);
        let caler = bondprop::order::Calculator::new(&mole);
        let omat = caler.pocv(&dirs, false, false, "mayer").expect("pocv 计算失败");
        data_verts.push((i, omat[(0, 2)].max(0.0).sqrt()));
    }

    // ③ 打印每个朝向的三种键级
    println!("--- 绕全局 z 轴扫描(72) ---");
    for i in 0..data.nrows() {
        println!(
            "{:>3}: MOCV-PX {:>10.6} MOCV-D2 {:>10.6} POCV-z {:>10.6}",
            i,
            data[(i, 0)],
            data[(i, 1)],
            data_pocv[i]
        );
    }

    // ④ 打印 wmview 方式的 36 个方向
    println!("--- 绕实际键轴 0->2 扫描(36, wmview方式) ---");
    for (i, val) in data_verts.iter() {
        println!("{:>3}: {:>10.6}", i, val);
    }

    // ⑤ 统计每列最小/最大，检查是否随朝向变化
    for (j, name) in ["PX", "D-2"].iter().enumerate() {
        let col = data.column(j);
        let mn = col.iter().cloned().fold(f64::INFINITY, f64::min);
        let mx = col.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        println!("MOCV {}: min={:.6} max={:.6} diff={:.2e}", name, mn, mx, mx - mn);
    }
    let mn = data_pocv.iter().cloned().fold(f64::INFINITY, f64::min);
    let mx = data_pocv.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    println!("POCV-z: min={:.6} max={:.6} diff={:.2e}", mn, mx, mx - mn);
    let mn = data_verts.iter().map(|x| x.1).fold(f64::INFINITY, f64::min);
    let mx = data_verts.iter().map(|x| x.1).fold(f64::NEG_INFINITY, f64::max);
    println!("POCV-verts: min={:.6} max={:.6} diff={:.2e}", mn, mx, mx - mn);

    // ⑥ wmview「方向键级.mayer.pocv.环轴方向.pi」输出的 36 组 (dx,dy,dz,val)，
    //    用相同方向在本地算 POCV 键级逐行对比（wm 的 val 是 raw.max(0).sqrt()）
    let wm: [(f64, f64, f64, f64); 36] = [
        (-0.4172, 0.9088, 0.0000, 0.4688),
        (-0.5686, 0.8226, 0.0000, 0.5869),
        (-0.7028, 0.7114, 0.0000, 0.6793),
        (-0.8157, 0.5785, 0.0000, 0.7349),
        (-0.9037, 0.4281, 0.0000, 0.7470),
        (-0.9643, 0.2646, 0.0000, 0.7141),
        (-0.9957, 0.0932, 0.0000, 0.6402),
        (-0.9967, -0.0812, 0.0000, 0.5342),
        (-0.9675, -0.2530, 0.0000, 0.4089),
        (-0.9088, -0.4172, 0.0000, 0.2794),
        (-0.8226, -0.5686, 0.0000, 0.1613),
        (-0.7114, -0.7028, 0.0000, 0.0689),
        (-0.5785, -0.8157, 0.0000, 0.0133),
        (-0.4281, -0.9037, 0.0000, 0.0012),
        (-0.2646, -0.9643, 0.0000, 0.0341),
        (-0.0932, -0.9957, 0.0000, 0.1080),
        (0.0812, -0.9967, 0.0000, 0.2140),
        (0.2530, -0.9675, 0.0000, 0.3393),
        (0.4172, -0.9088, 0.0000, 0.4688),
        (0.5686, -0.8226, 0.0000, 0.5869),
        (0.7028, -0.7114, 0.0000, 0.6793),
        (0.8157, -0.5785, 0.0000, 0.7349),
        (0.9037, -0.4281, 0.0000, 0.7470),
        (0.9643, -0.2646, 0.0000, 0.7141),
        (0.9957, -0.0932, 0.0000, 0.6402),
        (0.9967, 0.0812, 0.0000, 0.5342),
        (0.9675, 0.2530, 0.0000, 0.4089),
        (0.9088, 0.4172, 0.0000, 0.2794),
        (0.8226, 0.5686, 0.0000, 0.1613),
        (0.7114, 0.7028, 0.0000, 0.0689),
        (0.5785, 0.8157, 0.0000, 0.0133),
        (0.4281, 0.9037, 0.0000, 0.0012),
        (0.2646, 0.9643, 0.0000, 0.0341),
        (0.0932, 0.9957, 0.0000, 0.1080),
        (-0.0812, 0.9967, 0.0000, 0.2140),
        (-0.2530, 0.9675, 0.0000, 0.3393),
    ];
    let caler = bondprop::order::Calculator::new(&mole);
    println!(
        "{:>3} {:>8} {:>8} {:>9} {:>9} {:>9} {:>9} {:>9}",
        "#", "dx", "dy", "wm-val", "raw(0,2)", "raw(1,2)", "sqrt(0,2)", "sqrt(1,2)"
    );
    for (i, (dx, dy, _dz, wval)) in wm.iter().enumerate() {
        let dir = Vec3::new(*dx, *dy, 0.0).normalize();
        let dirs02: HashMap<usize, Vec3> = HashMap::from([(0usize, dir.clone()), (2usize, dir.clone())]);
        let omat02 = caler.pocv(&dirs02, false, false, "mayer").expect("pocv 失败");
        let dirs12: HashMap<usize, Vec3> = HashMap::from([(1usize, dir.clone()), (2usize, dir.clone())]);
        let omat12 = caler.pocv(&dirs12, false, false, "mayer").expect("pocv 失败");
        let v02 = omat02[(0, 2)];
        let v12 = omat12[(1, 2)];
        println!(
            "{:>3} {:>8.4} {:>8.4} {:>9.4} {:>9.4} {:>9.4} {:>9.4} {:>9.4}",
            i,
            dx,
            dy,
            wval,
            v02,
            v12,
            v02.max(0.0).sqrt(),
            v12.max(0.0).sqrt()
        );
    }
}

fn test8() {
    let lmat = Array2::eye(3);
    let nums = [3, 6, 10, 15];
    for num in nums {
        let tmat = mocv::get_tmat_car(num, &lmat);
        println!("num:{} tmat:\n{}", num, tmat);
    }
    let numt = [3, 5, 7, 9];
    for num in numt {
        let tmat = mocv::get_tmat_sph(num, &lmat);
        println!("num:{} tmat:\n{}", num, tmat);
    }
}

/// 运行所有测试的入口函数
pub fn test() {
    loger::init_log();
    let t0 = std::time::Instant::now();
    test_reader();
    println!("消耗时间 time: {:?}", t0.elapsed());
}
