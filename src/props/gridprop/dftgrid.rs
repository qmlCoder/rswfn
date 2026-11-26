#![allow(non_snake_case)]
use ndarray::Array2;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::f64::consts::PI;

use crate::base::Mole;
use crate::datas::{elems, lebedev};
use crate::maths;
// use crate::utils;

/// DFT 网格计算器，用于生成球面、径向和分子格点
pub struct Calculator<'m> {
    pub mole: &'m Mole,
}

impl<'m> Calculator<'m> {
    /// 创建一个新的 Calculator 实例
    pub fn new(mole: &'m Mole) -> Self {
        Calculator { mole }
    }

    // 计算原子格点，径向格点与球面格点的结合
    #[rustfmt::skip]
    /// 计算单个原子的三维网格（径向×球面）
    pub fn atm_grids(
        &self,
        iatm: usize,
        nrad: usize, // 径向格点数量
        nsph: usize, // 球面格点函数
    ) -> (Vec<[f64; 3]>, Vec<f64>) {
        let (sph_grids, sph_weits) = sph_grids(nsph);
        let nsph = sph_grids.len(); //球面格点数量
        let atom = self.mole.atoms().get(iatm);
        let [ax, ay, az] = atom.xyz;
        let mut grids = vec![[0.0; 3]; nrad * nsph];
        let mut weits = vec![0.0; nrad * nsph];
        let R = atom.rad;
        let (rs, ws) = rad_grids(R, nrad);
        for i in 0..nrad {
            let r = rs[i]; // 径向距离
            for j in 0..nsph {
                let [x, y, z] = sph_grids[j];
                grids[i * nsph + j] = [x * r + ax, y * r + ay, z * r + az];
                weits[i * nsph + j] = ws[i] * sph_weits[j];
            }
        }
        return (grids, weits);
    }

    #[rustfmt::skip]
    /// 计算指定元素的原子网格（根据核电荷数）
    pub fn elm_grids(
        &self,
        nuc: usize,
        nrad: usize, // 径向格点数量
        nsph: usize, // 球面格点函数
    ) -> (Vec<[f64; 3]>, Vec<f64>) {
        let (sph_grids, sph_weits) = sph_grids(nsph);
        let nsph = sph_grids.len(); //球面格点数量
        let mut grids = vec![[0.0; 3]; nrad * nsph];
        let mut weits = vec![0.0; nrad * nsph];
        let R = elems::atm2rad(nuc);
        let (rs, ws) = rad_grids(R, nrad);
        for i in 0..nrad {
            let r = rs[i]; // 径向距离
            for j in 0..nsph {
                let [x, y, z] = sph_grids[j];
                grids[i * nsph + j] = [x * r, y * r, z * r ];
                weits[i * nsph + j] = ws[i] * sph_weits[j];
            }
        }
        return (grids, weits);
    }

    // 将原子的权重映射到分子中 Becke分区方法
    /// 使用 Becke 分区方法将原子权重映射到分子中
    pub fn a2m_weits(
        &self,
        iatm: usize,               // 第多少个原子
        atm_grids: &Vec<[f64; 3]>, // 原子格点
        atm_weits: &Vec<f64>,      // 原子格点权重
    ) -> Vec<f64> {
        let atoms = self.mole.atoms();
        let ngrid = atm_grids.len();
        let mut a2m_weits = vec![0.0; ngrid];
        let atm_xyzs = atoms.xyzs();
        let atm_rads = atoms.rads();
        let natm = atm_xyzs.len(); // 原子数量
        let mut S_u = vec![vec![1.0; natm]; natm];
        for g in 0..ngrid {
            let gp = atm_grids[g];
            for i in 0..natm {
                let ip = atm_xyzs[i];
                let ri = maths::calc_dist(&ip, &gp); // 格点到原子i的距离
                for j in 0..natm {
                    if i == j {
                        continue;
                    }
                    let jp = atm_xyzs[j];
                    let rj = maths::calc_dist(&jp, &gp); // 格点到原子j的距离
                    let dist_ij = maths::calc_dist(&ip, &jp); // 两个原子之间的距离
                    let miu_ij = (ri - rj) / dist_ij;
                    let chi = atm_rads[i] / atm_rads[j]; // 两个原子的半径之比
                    let mut nu_ij;
                    if (chi - 1.0).abs() < 1e-6 {
                        nu_ij = miu_ij;
                    } else {
                        let u_ij = (chi - 1.0) / (chi + 1.0);
                        let mut a_ij = u_ij / (u_ij.powf(2.0) - 1.0);
                        if a_ij > 0.5 {
                            a_ij = 0.5
                        };
                        if a_ij < -0.5 {
                            a_ij = -0.5
                        };
                        nu_ij = miu_ij + a_ij * (1.0 - miu_ij.powf(2.0));
                    }
                    nu_ij = 1.5 * nu_ij - 0.5 * nu_ij.powi(3);
                    nu_ij = 1.5 * nu_ij - 0.5 * nu_ij.powi(3);
                    nu_ij = 1.5 * nu_ij - 0.5 * nu_ij.powi(3);
                    S_u[j][i] = 0.5 * (1.0 - nu_ij);
                }
            }
            // 计算权重乘积
            let mut wt = vec![1.0; natm];
            for j in 0..natm {
                for i in 0..natm {
                    wt[j] *= S_u[i][j];
                }
            }

            let sum_wt: f64 = wt.iter().sum();
            let rat = if sum_wt.abs() > 1e-12 {
                wt[iatm as usize] / sum_wt
            } else {
                0.0
            };
            a2m_weits[g] = atm_weits[g] * rat;
        }
        return a2m_weits;
    }

    // 计算分子格点 nrad: 径向格点数量； fsph: 球面格点函数
    /// 计算全分子格点及其权重
    pub fn mol_grids(&self, nrad: usize, fsph: usize) -> (Vec<[f64; 3]>, Vec<f64>) {
        let natm = self.mole.atoms().len();
        let mut mol_grids: Vec<[f64; 3]> = Vec::new();
        let mut mol_weits: Vec<f64> = Vec::new();
        let datas: Vec<(Vec<[f64; 3]>, Vec<f64>)> = (0..natm)
            .into_par_iter()
            .map(|i| {
                let (atm_grids, atm_weits) = self.atm_grids(i as usize, nrad, fsph);
                let a2m_weits = self.a2m_weits(i as usize, &atm_grids, &atm_weits);
                (atm_grids, a2m_weits)
            })
            .collect();
        for (atm_grids, a2m_weits) in datas {
            mol_grids.extend(&atm_grids);
            mol_weits.extend(&a2m_weits);
        }
        return (mol_grids, mol_weits);
    }

    /// 计算指定原子片段的网格和权重
    pub fn frag_grids(&self, frag: &Vec<usize>) -> (Vec<[f64; 3]>, Vec<f64>) {
        let mut grids = Vec::new();
        let mut weits = Vec::new();
        for atm in frag {
            let (atm_grids, atm_weits) = self.atm_grids(*atm, 10, 74);
            let a2m_weits = self.a2m_weits(*atm, &atm_grids, &atm_weits);
            grids.extend(atm_grids);
            weits.extend(a2m_weits);
        }
        (grids, weits)
    }
}

// 球面格点
/// 获取 Lebedev 球面格点及其权重
/// nsph: 球面格点数 可选值为: 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590
pub fn sph_grids(nsph: usize) -> (Vec<[f64; 3]>, Vec<f64>) {
    match nsph {
        006 => lebedev::ld0006(),
        014 => lebedev::ld0014(),
        026 => lebedev::ld0026(),
        038 => lebedev::ld0038(),
        050 => lebedev::ld0050(),
        074 => lebedev::ld0074(),
        086 => lebedev::ld0086(),
        110 => lebedev::ld0110(),
        146 => lebedev::ld0146(),
        170 => lebedev::ld0170(),
        194 => lebedev::ld0194(),
        230 => lebedev::ld0230(),
        266 => lebedev::ld0266(),
        302 => lebedev::ld0302(),
        350 => lebedev::ld0350(),
        434 => lebedev::ld0434(),
        590 => lebedev::ld0590(),
        _ => {
            log::error!(
                "非法nsph:{}\n可选sph为 6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590",
                nsph
            );
            lebedev::ld0074()
        }
    }
}

// 计算dft原子径向格点
/// 计算 DFT 原子径向格点坐标和权重
pub fn rad_grids(R: f64, nrad: usize) -> (Vec<f64>, Vec<f64>) {
    // let R = 1.0;
    let mut rs = vec![0.0; nrad];
    let mut ws = vec![0.0; nrad];
    for i in 1..=nrad {
        let xi = f64::cos(i as f64 * PI / (nrad as f64 + 1.0));
        let ri = R * (1.0 + xi) / (1.0 - xi);
        let wi = 2.0 * PI / (nrad as f64 + 1.0) * R.powi(3) * (1.0 + xi).powf(2.5) / (1.0 - xi).powf(3.5);
        rs[i - 1] = ri;
        ws[i - 1] = wi * 4.0 * PI;
    }
    return (rs, ws);
}
