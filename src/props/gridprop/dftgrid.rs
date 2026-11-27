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
    /// 分子信息，提供原子坐标与共价半径
    pub mole: &'m Mole,
}

impl<'m> Calculator<'m> {
    /// 创建一个新的 Calculator 实例
    pub fn new(mole: &'m Mole) -> Self {
        Calculator { mole }
    }

    // 计算原子格点，径向格点与球面格点的结合
    #[rustfmt::skip]
    /// 计算单个原子的三维网格（径向×球面，不剪枝）
    pub fn atm_grids(
        &self,
        iatm: usize,
        nrad: usize, // 径向格点数量
        nsph: usize, // 球面格点函数
    ) -> (Vec<[f64; 3]>, Vec<f64>) {
        // let R = self.mole.atoms().get(iatm).rad;
        // atm_grids_pruned(self.mole.atoms().get(iatm).rad, self.mole.atoms().get(iatm).xyz, nrad, nsph, false)
        self.atm_grids_p(iatm, nrad, nsph, false)
    }

    /// 计算单个原子的三维网格，可按径向位置剪枝球面格点数
    ///
    /// 剪枝规则：靠近核与远离核的径向壳层用较粗的球面网格（这两处密度要么
    /// 近球形、要么已经很小），中间价层区用完整的球面网格。这样能在几乎不
    /// 改变积分精度的前提下把格点数降到约一半。
    pub fn atm_grids_p(
        &self,
        iatm: usize,
        nrad: usize,
        nsph: usize,
        prune: bool,
    ) -> (Vec<[f64; 3]>, Vec<f64>) {
        let (grids, weits, _rs) = self.atm_grids_meta(iatm, nrad, nsph, prune);
        (grids, weits)
    }

    /// 计算单个原子的三维网格，同时返回每个格点的径向半径（Bohr）
    ///
    /// 半径用于后续"每个格点块只取显著的原子轨道"的筛选：块内格点都在以
    /// 原子为中心的球壳上，因此到其它原子中心的最小距离可以由块半径区间
    /// 精确给出。
    pub fn atm_grids_meta(
        &self,
        iatm: usize,
        nrad: usize,
        nsph: usize,
        prune: bool,
    ) -> (Vec<[f64; 3]>, Vec<f64>, Vec<f64>) {
        let atom = self.mole.atoms().get(iatm);
        let [ax, ay, az] = atom.xyz;
        // 径向参数：共价半径乘以放大因子。Mura-Knowles 变换把 (0,∞) 映到
        // (-1,1)，R 只影响格点分布：R 偏小时尾部点过稀，分子越大积分误差越明显。
        let scale = std::env::var("RSWFN_RAD_SCALE")
            .ok()
            .and_then(|v| v.parse::<f64>().ok())
            .unwrap_or(1.0);
        let R = atom.rad * scale;
        let (rs, ws) = rad_grids(R, nrad);
        let mut grids: Vec<[f64; 3]> = Vec::new();
        let mut weits: Vec<f64> = Vec::new();
        let mut radii: Vec<f64> = Vec::new();
        if !prune {
            grids.reserve(nrad * nsph);
            weits.reserve(nrad * nsph);
            radii.reserve(nrad * nsph);
        }
        for i in 0..nrad {
            let nang = if prune { prune_nang_r(rs[i], nsph) } else { nsph };
            let (sph_grids, sph_weits) = sph_grids(nang);
            let r = rs[i];
            for j in 0..sph_grids.len() {
                let [x, y, z] = sph_grids[j];
                grids.push([x * r + ax, y * r + ay, z * r + az]);
                weits.push(ws[i] * sph_weits[j]);
                radii.push(r);
            }
        }
        (grids, weits, radii)
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
    /// 计算全分子格点及其权重（默认不剪枝，保证精度）
    pub fn mol_grids(&self, nrad: usize, fsph: usize) -> (Vec<[f64; 3]>, Vec<f64>) {
        self.mol_grids_p(nrad, fsph, false)
    }

    /// 计算全分子格点及其权重，`prune` 控制是否做径向-角向剪枝
    pub fn mol_grids_p(
        &self,
        nrad: usize,
        fsph: usize,
        prune: bool,
    ) -> (Vec<[f64; 3]>, Vec<f64>) {
        let (grids, weits, _atm, _radii) = self.mol_grids_meta(nrad, fsph, prune);
        (grids, weits)
    }

    /// 计算全分子格点，并附带每个格点的所属原子索引与径向半径
    ///
    /// 返回 `(格点坐标, 权重, 原子索引, 半径)`，四者一一对应。原子索引与半径
    /// 供 KS 的格点分块筛选使用（每块只保留近邻原子的显著轨道）。
    pub fn mol_grids_meta(
        &self,
        nrad: usize,
        fsph: usize,
        prune: bool,
    ) -> (Vec<[f64; 3]>, Vec<f64>, Vec<usize>, Vec<f64>) {
        let natm = self.mole.atoms().len();
        let datas: Vec<(Vec<[f64; 3]>, Vec<f64>, Vec<f64>)> = (0..natm)
            .into_par_iter()
            .map(|i| {
                let (atm_grids, atm_weits, radii) = self.atm_grids_meta(i as usize, nrad, fsph, prune);
                let a2m_weits = self.a2m_weits(i as usize, &atm_grids, &atm_weits);
                (atm_grids, a2m_weits, radii)
            })
            .collect();
        let mut mol_grids: Vec<[f64; 3]> = Vec::new();
        let mut mol_weits: Vec<f64> = Vec::new();
        let mut mol_atm: Vec<usize> = Vec::new();
        let mut mol_radii: Vec<f64> = Vec::new();
        for (i, (atm_grids, a2m_weits, radii)) in datas.into_iter().enumerate() {
            mol_grids.extend(&atm_grids);
            mol_weits.extend(&a2m_weits);
            mol_atm.extend(std::iter::repeat(i).take(atm_grids.len()));
            mol_radii.extend(&radii);
        }
        return (mol_grids, mol_weits, mol_atm, mol_radii);
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

/// 径向-角向剪枝：按径向壳层的**实际半径**（Bohr）选取球面格点数
///
/// 近核区（r ≤ 0.4）密度近似球对称、远核区（r > 10）权重已经很小，两者都
/// 用 50 点即可；0.4 < r ≤ 1.5 的内层价区用 194 点；1.5 < r ≤ 10 的价层与
/// 成键区用完整球面网格（`nmax` 点）。
///
/// 这样在 (99,590) 这种高精度设置下，格点数约为不剪枝的 1/4，而能量与
/// 不剪枝结果相差在 1e-6 Hartree 量级。
fn prune_nang_r(r: f64, nmax: usize) -> usize {
    let cand = if r <= 0.4 || r > 10.0 {
        50
    } else if r <= 1.5 {
        194
    } else {
        nmax
    };
    cand.min(nmax).max(50)
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
