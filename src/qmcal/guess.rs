//! SCF 初猜：叠加原子密度（SAD-style）
//!
//! 经典做法是把"整分子核心哈密顿量 + 上标"作为初猜，但对长链/大分子，这个
//! 初猜几乎没有原子壳层结构，前几轮 SCF 会大幅摆动（C18H38 的 B3LYP 实测
//! 需要 28 轮，而 g16 用叠加原子密度只需 11 轮）。
//!
//! 这里按原子切开：每个原子在自己的 AO 子空间里对 `H_a` 做广义本征值分解
//! （用 `S_a^{-1/2}` 正交化），再按该原子的电子数填充，最后把各原子块直和
//! 成整分子密度。它保留了孤对/价层的原子结构，是 Gaussian 默认 SAD 初猜的
//! 廉价近似（不含原子密度表，直接由基组现算）。

use ndarray::Array2;

use crate::base::Mole;
use crate::maths;

/// 每个原子轨道所属的原子索引（顺序与 `nato` 一致）
pub fn ao_owners(mole: &Mole, form: &str) -> Vec<usize> {
    let atos = mole.basis().atos();
    let list = if form == "car" { &atos.cars } else { &atos.sphs };
    list.iter().map(|a| a.atm).collect()
}

/// 是否使用 SAD 初猜（环境变量 `RSWFN_GUESS=core` 可退回核心哈密顿量）
pub fn use_sad_enabled() -> bool {
    std::env::var("RSWFN_GUESS").map(|v| v != "core").unwrap_or(true)
}

/// 由"逐原子求解核哈密顿量块"构造叠加原子密度矩阵
///
/// - `smat`/`hmat`：整分子的重叠矩阵与核心哈密顿量
/// - `ao_atm`：每个原子轨道所属的原子索引
/// - `nucs`：各原子核电荷（用于按比例分配总电子数）
/// - `nele`：分子总电子数
///
/// 返回 `None` 表示某个原子块无法正交化（基组近线性相关），调用方应退回
/// 核心哈密顿量初猜。
pub fn sad_dmat(
    smat: &Array2<f64>,
    hmat: &Array2<f64>,
    ao_atm: &[usize],
    nucs: &[usize],
    nele: usize,
) -> Option<Array2<f64>> {
    let nato = smat.shape()[0];
    let natm = nucs.len();
    let zsum: usize = nucs.iter().sum();
    if nato == 0 || natm == 0 || zsum == 0 || nele == 0 {
        return None;
    }
    let mut dmat = Array2::<f64>::zeros((nato, nato));
    for atm in 0..natm {
        let idx: Vec<usize> = (0..nato).filter(|&i| ao_atm[i] == atm).collect();
        let m = idx.len();
        if m == 0 {
            continue;
        }
        let mut s = Array2::<f64>::zeros((m, m));
        let mut h = Array2::<f64>::zeros((m, m));
        for (a, &i) in idx.iter().enumerate() {
            for (b, &j) in idx.iter().enumerate() {
                s[(a, b)] = smat[(i, j)];
                h[(a, b)] = hmat[(i, j)];
            }
        }
        let x = maths::matrix_nhalf(&s);
        for v in x.iter() {
            if !v.is_finite() {
                return None;
            }
        }
        let f = x.t().dot(&h).dot(&x);
        let (_engs, cmat) = maths::eigh(&f);
        let c = x.dot(&cmat);
        // 该原子的电子数：按核电荷比例分摊总电子数（中性分子即 Z_a）
        let nel = (nucs[atm] as f64 * nele as f64 / zsum as f64).round() as usize;
        let mut left = nel;
        for k in 0..m {
            if left == 0 {
                break;
            }
            let occ = left.min(2) as f64;
            left -= occ as usize;
            for a in 0..m {
                let ca = c[(a, k)] * occ;
                if ca == 0.0 {
                    continue;
                }
                for b in 0..m {
                    dmat[(idx[a], idx[b])] += ca * c[(b, k)];
                }
            }
        }
    }
    Some(dmat)
}
