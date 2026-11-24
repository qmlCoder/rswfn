// DHMO方法的第二个版本，计算任意朝向P轨道的重叠积分，其实就是在EHMO的基础上只考虑pz轨道
// 要能够计算两个p轨道的重叠积分，这两个原子轨道都在x轴上，积分数值只和距离有关

use crate::bondprop;
use crate::datas::dhmo2::get_paras;
use crate::{
    atomprop,
    base::{Mole, Stm, Vec3},
    datas::consts::BOHR,
    maths, utils,
};
use ndarray::{Array2, s};

#[rustfmt::skip]
const XS: [f64; 60] = [
    0.00000000, 0.05084746, 0.10169492, 0.15254237, 0.20338983, 0.25423729, 0.30508475, 0.35593220, 0.40677966, 0.45762712,
    0.50847458, 0.55932203, 0.61016949, 0.66101695, 0.71186441, 0.76271186, 0.81355932, 0.86440678, 0.91525424, 0.96610169,
    1.01694915, 1.06779661, 1.11864407, 1.16949153, 1.22033898, 1.27118644, 1.32203390, 1.37288136, 1.42372881, 1.47457627,
    1.52542373, 1.57627119, 1.62711864, 1.67796610, 1.72881356, 1.77966102, 1.83050847, 1.88135593, 1.93220339, 1.98305085,
    2.03389831, 2.08474576, 2.13559322, 2.18644068, 2.23728814, 2.28813559, 2.33898305, 2.38983051, 2.44067797, 2.49152542,
    2.54237288, 2.59322034, 2.64406780, 2.69491525, 2.74576271, 2.79661017, 2.84745763, 2.89830508, 2.94915254, 3.00000000,
];

/// 计算两个位于x轴上的基函数的重叠积分
///
/// 使用三次样条插值方法计算重叠积分。当距离大于3.0时返回0.0。
///
/// # 参数
/// - `zi`, `zj`: 两个原子的核电荷数
/// - `pi`, `pj`: 轨道角动量类型
/// - `dist`: 原子间距离
pub fn get_integ_x(zi: usize, zj: usize, pi: usize, pj: usize, dist: f64) -> f64 {
    if dist > 3.0 {
        return 0.0;
    }
    if dist == 0.0 {
        return 1.0;
    }
    let paras = get_paras(zi, zj, pi, pj); // 获取三次样条插值的参数
    let mut idx = 0;
    for i in 0..59 {
        if dist >= XS[i] && dist < XS[i + 1] {
            idx = i;
        }
    }
    let x = dist - XS[idx];
    let [a, b, c, d] = paras[idx];
    ((a * x + b) * x + c) * x + d
}

/// 计算两个任意朝向的Pz轨道的重叠积分
///
/// 通过将轨道投影到局部坐标系的三个轴上,分别计算各分量上的重叠积分并求和。
pub fn get_integ_d(zi: usize, zj: usize, di: &Vec3, dj: &Vec3, stm: Stm, dist: f64) -> f64 {
    // println!("Di:{}", Di.norm());
    // println!("Dj:{}", Dj.norm());
    // println!("Dj:{:?}", Dj);
    // println!("stm:{:?}", stm);
    if dist > 3.0 {
        return 0.0;
    }
    let ps = [0, 2, 2]; // px pz pz
    let mut val = 0.0;
    let axs = vec![Vec3::from_arr(stm.ex), Vec3::from_arr(stm.ey), Vec3::from_arr(stm.ez)];

    for i in 0..3 {
        let axi = &axs[i];
        // println!("axi:{}", axi.norm());
        let ci = di.dot(&axi);
        let cj = dj.dot(&axi);
        let sv = get_integ_x(zi, zj, ps[i], ps[i], dist);
        // println!("{} {:>10.4} {:>10.4} {:>10.4}", i, sv, ci, cj);
        val += sv * ci * cj;
    }
    // println!("val:{:>10.4}", val);
    val
}

/// dHMO方法版本2,计算分子中的所有原子轨道
///
/// 基于扩展Hückel方法(Extended Hückel),构建重叠矩阵和哈密顿矩阵,
/// 求解广义本征值问题,计算Mayer键级和电子布局。
pub fn dhmo_calc_v2(mole: &Mole) {
    let dir_caler = atomprop::direction::Calculator::new(mole);
    let atoms = mole.atoms();
    let natm = atoms.len();
    let mut ato_atms = vec![];
    let mut ato_dirs = vec![];
    // let ato_stms = vec![];
    // 原子轨道朝向
    for atm in 0..natm {
        if atoms.get(atm).sym == "H" {
            continue;
        }
        let dir = dir_caler.normal_vector(atm);
        if let Some(dir) = dir {
            let [dx, dy, dz] = dir.into_arr();
            println!("{}:{:>10.4}{:>10.4}{:>10.4}", atm + 1, dx, dy, dz);
            ato_atms.push(atm);
            ato_dirs.push(dir);
        }
    }
    // 计算重叠矩阵
    let nato = ato_atms.len();
    let nucs = atoms.nucs();
    let mut smat = Array2::eye(nato);
    for i in 0..nato - 1 {
        let ai = ato_atms[i];
        for j in i + 1..nato {
            let aj = ato_atms[j];
            // 创建一个局部坐标系，x轴在键轴上剩下两个轴应该对结果没影响，因为两个原子轨道使用同一个朝向的局部在坐标系，只是在x轴上有一定距离
            let vx = atoms.vec(ai, aj).normalize();
            let vr = Vec3::rand();
            let vz = vx.cross(&vr).normalize();
            let vy = vz.cross(&vx).normalize();
            let stm = Stm::new(vx.into_arr(), vy.into_arr(), vz.into_arr());
            let Zi = nucs[ai];
            let Zj = nucs[aj];
            let Di = &ato_dirs[i];
            let Dj = &ato_dirs[j];
            let dist = atoms.dist(ai, aj) / BOHR;

            let sval = get_integ_d(Zi, Zj, &Di, &Dj, stm, dist);
            println!("{}-{}:{:>10.4},{:>10.4}", ai + 1, aj + 1, dist, sval);
            smat[(i, j)] = sval;
            smat[(j, i)] = sval;
        }
    }
    println!("smat:");
    utils::show_matrix(&smat);
    utils::show_matrix(&(&smat * &smat.t()));
    let mut hmat = Array2::zeros((nato, nato));
    let K = 1.75;
    let Is = vec![-19.4; nato];
    for i in 0..nato {
        for j in 0..nato {
            if i == j {
                hmat[(i, j)] = Is[i];
            } else {
                hmat[(i, j)] = K * smat[(i, j)] / 2.0 * (Is[i] + Is[j]);
            }
        }
    }
    println!("hmat:");
    utils::show_matrix(&hmat);
    let xmat = maths::matrix_nhalf(&smat);
    let hmat_ = xmat.t().dot(&hmat).dot(&xmat);
    let (engs, cmat_) = maths::eigh(&hmat);
    let cmat = xmat.dot(&cmat_);
    println!("cmat:");
    utils::show_matrix(&cmat);
    println!("engs:{}", engs);
    let cscm = cmat.t().dot(&smat).dot(&cmat);
    println!("cscm:");
    utils::show_matrix(&cscm);
    let cmat = cmat.slice(s![.., ..nato / 2]).to_owned();
    let dmat = cmat.dot(&cmat.t()) * 2.0;
    let qmat = &dmat * &smat;
    println!("dmat:");
    utils::show_matrix(&dmat);
    println!("qmat:");
    utils::show_matrix(&qmat);
    println!("电子布局");
    for i in 0..nato {
        let mut sum = 0.0;
        for j in 0..nato {
            sum += qmat[(i, j)];
        }
        println!("{} {:>10.4}", i + 1, sum);
    }
    println!("mayer 键级");
    let mut uplo = vec![];
    for i in 0..nato {
        uplo.push([i, i + 1]);
    }
    let omat_ato = bondprop::order::mayer(&dmat, &smat);
    let omat_atm = bondprop::order::ato2atm(&omat_ato, &uplo);
    println!("omat_ato:");
    utils::show_matrix(&omat_ato);
    println!("omat_atm:");
    utils::show_matrix(&omat_atm);
    for i in 0..natm - 1 {
        for j in i + 1..natm {
            let val = omat_atm[(i, j)].abs().sqrt();
            if val.abs() < 1e-5 {
                continue;
            }
            println!("{} {} {:>10.4}", ato_atms[i] + 1, ato_atms[j] + 1, val);
        }
    }
}
