#![allow(non_snake_case)]
// 此模块定义积分函数
//
// 单/双电子积分统一由 rsint crate(libcint 的 Rust 直译版)提供，
// 壳层块引擎 [`Engine`] 负责 SCF 相关的双电子部分。
use crate::base::{Mole, Stm};
use crate::datas::btrans;
use crate::props::orbtprop::mocv;
use ndarray::{Array2, s};
use rayon::prelude::*;
use std::collections::HashMap;

// 电子静电势的解析计算，与后端无关（不依赖 libcint 或手写积分的选择）
mod esp;
pub use esp::EspCalculator;

mod rsint;

// 壳层块积分引擎：SCF 的双电子部分走这里（依赖 rsint crate 的 C ABI）
mod engine;
pub use engine::Engine;

/// 双电子积分压缩存储结构体
///
/// 利用对称性将四维双电子积分矩阵压缩为一维存储。
pub struct Mat2e {
    /// 原子轨道数量
    pub nato: usize,
    /// 积分值列表
    pub vals: Vec<f64>,
}

impl Mat2e {
    /// 创建新的Mat2e
    pub fn new(nato: usize) -> Self {
        let nval = nato * (nato + 1) * (nato * (nato + 1) + 2) / 8; // 有效电子积分数量
        let vals = vec![0.0; nval];
        Self { nato, vals }
    }

    /// 获取积分值
    pub fn get(&self, i: usize, j: usize, k: usize, l: usize) -> f64 {
        let idx = get_2eid(i, j, k, l);
        self.vals[idx]
    }

    /// 设置积分值
    pub fn set(&mut self, i: usize, j: usize, k: usize, l: usize, val: f64) {
        let idx = get_2eid(i, j, k, l);
        self.vals[idx] = val;
    }

    /// 按两个原子轨道对的压缩索引取值 (p = (ij), q = (kl))
    ///
    /// 与连续调用 [`Mat2e::get`] 等价，但省去了四维索引的重复折算，
    /// 适合自洽场中反复构造库伦/交换矩阵的场景。
    pub fn get_by_pair(&self, p: usize, q: usize) -> f64 {
        let idx = if p >= q { p * (p + 1) / 2 + q } else { q * (q + 1) / 2 + p };
        self.vals[idx]
    }
}

/// Rust实现的分子积分计算器
///
/// 提供单电子积分(重叠、动能、核引)和双电子积分的计算功能,
/// 以及球谐与笛卡尔基函数的转换。
pub struct RsInt<'m> {
    /// 分子信息引用
    mole: &'m Mole,
}

impl<'m> RsInt<'m> {
    /// 创建新的RsInt实例
    ///
    /// # Parameters
    /// - `mole`: 分子信息，提供原子结构与基组用于积分计算
    pub fn new(mole: &'m Mole) -> Self {
        Self { mole }
    }

    /// 计算重叠矩阵
    pub fn smat(&self, form: &str) -> Array2<f64> {
        self.mat1e([1, 0, 0], form)[0].clone()
    }

    /// 计算动能矩阵
    pub fn tmat(&self, form: &str) -> Array2<f64> {
        self.mat1e([0, 1, 0], form)[1].clone()
    }

    /// 计算核吸引势矩阵
    pub fn vmat(&self, form: &str) -> Array2<f64> {
        self.mat1e([0, 0, 1], form)[2].clone()
    }

    /// 计算单电子积分(使用 rsint crate)
    /// cals 积分类型 [s, t, v]，1 表示计算该类型积分，0 表示不计算
    /// form 积分形式，"raw" 表示笛卡尔坐标，"sph" 表示球谐坐标
    pub fn mat1e(&self, cals: [usize; 3], form: &str) -> [Array2<f64>; 3] {
        let rsint = rsint::Rsint::new(self.mole);
        rsint.calc_int1e(form, cals)
    }

    /// 将笛卡尔单电子积分转换为球谐形式
    pub fn mat1e_sph(&self, car_mat: &Array2<f64>) -> Array2<f64> {
        let shls = &self.mole.basis().shls;
        let nshl = shls.len();
        let uplo_sph = shls.uplo("sph");
        let uplo_car = shls.uplo("car");
        let nato_sph = shls.nato("sph");
        let mut sph_mat = Array2::zeros((nato_sph, nato_sph));
        for i in 0..nshl {
            let angi = shls.get(i).ang;
            let [cui, cli] = uplo_car[i];
            let [sui, sli] = uplo_sph[i];
            let tmati = btrans::car2sph_tmat(angi);
            for j in 0..nshl {
                let angj = shls.get(j).ang;
                let [cuj, clj] = uplo_car[j];
                let [suj, slj] = uplo_sph[j];
                let tmatj = btrans::sph2car_tmat(angj);
                let car_part = car_mat.slice(s![cui..cli, cuj..clj]).to_owned();
                let sph_part: Array2<f64> = tmati.dot(&car_part);
                let sph_part = sph_part.dot(&tmatj);

                sph_mat.slice_mut(s![sui..sli, suj..slj]).assign(&sph_part);
            }
        }
        sph_mat
    }

    /// 将笛卡尔单电子积分转换为原始(非收缩)形式
    pub fn mat1e_raw(&self, car_mat: &Array2<f64>) -> Array2<f64> {
        let shls = &self.mole.basis().shls;
        let nshl = shls.len();
        let uplo_raw = shls.uplo("raw");
        let uplo_car = shls.uplo("car");
        let nato_raw = shls.nato("raw");
        let mut smat_raw = Array2::zeros((nato_raw, nato_raw));
        for i in 0..nshl {
            let ang_i = shls.get(i).ang;
            let [up_car_i, lo_car_i] = uplo_car[i];
            let [up_raw_i, lo_raw_i] = uplo_raw[i];
            let n_car_i = lo_car_i - up_car_i;
            let n_raw_i = lo_raw_i - up_raw_i;
            let tmati = btrans::sph2car_tmat(ang_i);
            for j in 0..nshl {
                let angj = shls.get(j).ang;
                let [up_car_j, lo_car_j] = uplo_car[j];
                let [up_raw_j, lo_raw_j] = uplo_raw[j];
                let n_car_j = lo_car_j - up_car_j;
                let n_raw_j = lo_raw_j - up_raw_j;
                let tmatj = btrans::sph2car_tmat(angj);
                let mut raw_part = car_mat.slice(s![up_car_i..lo_car_i, up_car_j..lo_car_j]).to_owned();
                if n_car_i != n_raw_i {
                    raw_part = tmati.t().dot(&raw_part)
                }
                if n_car_j != n_raw_j {
                    raw_part = raw_part.dot(&tmatj)
                }
                smat_raw.slice_mut(s![up_raw_i..lo_raw_i, up_raw_j..lo_raw_j]).assign(&raw_part);
            }
        }
        smat_raw
    }

    /// 获取局部坐标系的单电子积分
    pub fn mat1e_lcs(&self, stms: &HashMap<usize, Stm>, cals: [usize; 3]) -> [Array2<f64>; 3] {
        let shls = &self.mole.basis().shls;
        let uplo = shls.uplo("car");
        let nshl = shls.len();
        // let atms = &map_shl.atms;

        let mut Tmats = HashMap::new();
        for i in 0..nshl {
            let num = uplo[i][1] - uplo[i][0];
            let atm = shls.get(i).atm;
            let tmat = mocv::get_tmat_car(num, &stms[&atm].mat);
            Tmats.insert(i, tmat); // 局部坐标系系数到全局坐标系系数的转换矩阵
        }

        // 全局坐标系的重叠积分
        let [gsmat, gtmat, gvmat] = self.mat1e(cals, "car");
        let mut lsmat = gsmat.clone();
        let mut ltmat = gtmat.clone();
        let mut lvmat = gvmat.clone();

        for i in 0..nshl {
            let Mt_i = Tmats.get(&i).unwrap();
            let [ui, li] = uplo[i];

            for j in 0..nshl {
                let Mt_j = Tmats.get(&j).unwrap();
                let [uj, lj] = uplo[j];

                let info = s![ui..li, uj..lj];
                let gs = gsmat.slice(info);
                let gt = gtmat.slice(info);
                let gv = gvmat.slice(info);
                let gs: Array2<f64> = Mt_i.dot(&gs).dot(&Mt_j.t());
                let gt: Array2<f64> = Mt_i.dot(&gt).dot(&Mt_j.t());
                let gv: Array2<f64> = Mt_i.dot(&gv).dot(&Mt_j.t());
                lsmat.slice_mut(info).assign(&gs);
                ltmat.slice_mut(info).assign(&gt);
                lvmat.slice_mut(info).assign(&gv);
            }
        }
        [lsmat, ltmat, lvmat]
    }
}

/// 由压缩的对称轨道对索引 `p`（`0 ≤ j ≤ i < nato`）解码出 `(i, j)`
///
/// 其中 `p = i(i+1)/2 + j`。
fn pair_index_to_ij(p: usize) -> (usize, usize) {
    let mut i = (2.0 * p as f64).sqrt() as usize;
    while i * (i + 1) / 2 > p {
        i -= 1;
    }
    while (i + 1) * (i + 2) / 2 <= p {
        i += 1;
    }
    (i, p - i * (i + 1) / 2)
}

/// 由密度矩阵与双电子积分构造库伦矩阵 J 与交换矩阵 K
///
/// - `J_ij = Σ_kl (ij|kl) D_kl`
/// - `K_ij = Σ_kl (ik|jl) D_kl`
///
/// 按原子轨道对 `(i, j)` 并行；预先算好轨道对的压缩索引，
/// 避免在内层循环中重复折算四维下标。J、K 均为对称矩阵，只计算 `i ≥ j` 的一半。
pub fn jk_mat(dmat: &Array2<f64>, emat: &Mat2e) -> (Array2<f64>, Array2<f64>) {
    let nato = dmat.nrows();
    let mut pair = vec![0usize; nato * nato];
    for i in 0..nato {
        for j in 0..nato {
            pair[i * nato + j] = if i >= j { i * (i + 1) / 2 + j } else { j * (j + 1) / 2 + i };
        }
    }
    let dvals: Vec<f64> = dmat.iter().copied().collect();
    // J、K 均为对称矩阵，只需计算 i ≥ j 的一半，再镜像填充
    let npair = nato * (nato + 1) / 2;
    let jk: Vec<[f64; 2]> = (0..npair)
        .into_par_iter()
        .map(|p| {
            let (i, j) = pair_index_to_ij(p);
            let pij = pair[i * nato + j];
            let jrow = &pair[j * nato..j * nato + nato];
            let (mut jsum, mut ksum) = (0.0, 0.0);
            for k in 0..nato {
                let pik = pair[i * nato + k];
                let krow = k * nato;
                for l in 0..nato {
                    let d = dvals[krow + l];
                    if d == 0.0 {
                        continue;
                    }
                    jsum += emat.get_by_pair(pij, pair[krow + l]) * d;
                    ksum += emat.get_by_pair(pik, jrow[l]) * d;
                }
            }
            [jsum, ksum]
        })
        .collect();
    let mut jmat = Array2::zeros((nato, nato));
    let mut kmat = Array2::zeros((nato, nato));
    for (p, [jval, kval]) in jk.into_iter().enumerate() {
        let (i, j) = pair_index_to_ij(p);
        jmat[(i, j)] = jval;
        jmat[(j, i)] = jval;
        kmat[(i, j)] = kval;
        kmat[(j, i)] = kval;
    }
    (jmat, kmat)
}

/// 将双电子积分的四维索引压缩为一维索引
///
/// 利用(ij|kl)积分的对称性,将(i,j,k,l)四维索引映射到一维数组。
pub fn get_2eid(i: usize, j: usize, k: usize, l: usize) -> usize {
    let p = if i >= j { i * (i + 1) / 2 + j } else { j * (j + 1) / 2 + i };
    let q = if k >= l { k * (k + 1) / 2 + l } else { l * (l + 1) / 2 + k };
    if p >= q { p * (p + 1) / 2 + q } else { q * (q + 1) / 2 + p }
}
