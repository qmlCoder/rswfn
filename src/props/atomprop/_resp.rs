#![allow(non_snake_case)]
//! RESP（Restrained ElectroStatic Potential）原子电荷。
//!
//! 实现参考 Multiwfn 3.8 的 `population.f90` 中 `RESP` / `RESPiter` / `setMKpt`：
//!
//! 1. 用 MK（Merz-Kollman）方案在分子周围生成拟合点：每个原子按 van der Waals
//!    半径的 1.4/1.6/1.8/2.0 倍生成 4 层球面点，密度为 6 点/Å²；落在其它原子
//!    1.4 倍 vdW 球内的点被剔除。
//! 2. 在拟合点上用解析高斯积分计算静电势（核项 + 电子项，原子单位）。
//! 3. 以 `χ² = Σ(V_pt - Σ q_i/r_i)² + a Σ f(q_i)` 为目标做带约束的最小二乘，
//!    其中双曲限制函数 `f(q) = sqrt(q² + b²) - b`，因此法方程对角元额外加上
//!    `a / sqrt(q² + b²)`。法方程同时包含总电荷约束。
//! 4. 标准两阶段拟合：第一阶段对所有重原子用 `a = 0.0005`；第二阶段把非
//!    sp³/亚甲基碳及其所连氢的原子固定为第一阶段电荷，对剩余的 sp³ 碳和
//!    CH₃/CH₂/=CH₂ 上的氢用 `a = 0.001` 继续拟合，并约束同一碳上的氢等价。
//!
//! 与 Multiwfn 一致，双曲限制只施加在非氢原子上，氢原子的约束通过等价约束
//! （同一 CH₃/CH₂/=CH₂ 上的 H 电荷相同）体现。
//!
//! # Example
//! ```no_run
//! use rswfn::props::atomprop::resp::Calculator;
//! use rswfn::base::Mole;
//!
//! let mole = Mole::from_file("moles/C2H4.out");
//! let caler = Calculator::new(&mole);
//! let res = caler.resp().unwrap();
//! println!("{:?}", res.charges);
//! ```

use std::f64::consts::PI;

use ndarray::{Array1, Array2};
use rayon::prelude::*;

use crate::base::Mole;
use crate::datas::consts::BOHR;
use crate::maths;
use crate::props::gridprop::density;

/// 静电势类型。
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EspType {
    /// 核项 + 电子项（默认，得到的量是原子电荷）
    Total,
    /// 仅电子项（得到的量是电子布居的负值）
    Electronic,
}

/// RESP 拟合参数，默认值与 Multiwfn 3.8 一致。
#[derive(Clone, Debug)]
pub struct RespOptions {
    /// MK 各层半径相对 vdW 半径的倍数，默认 `[1.4, 1.6, 1.8, 2.0]`
    pub layer_scales: Vec<f64>,
    /// MK 拟合点密度（点/Å²），默认 6.0
    pub pt_density: f64,
    /// 第一阶段双曲限制强度，默认 0.0005
    pub hyper_a1: f64,
    /// 第二阶段双曲限制强度，默认 0.001
    pub hyper_a2: f64,
    /// 双曲限制的 tightness 参数 b，默认 0.1
    pub hyper_b: f64,
    /// 最大迭代次数，默认 300
    pub max_iter: usize,
    /// 电荷收敛阈值，默认 1e-6
    pub conv: f64,
    /// 是否执行标准两阶段拟合。`false` 时只执行第一阶段
    pub two_stage: bool,
    /// 静电势类型
    pub esp_type: EspType,
}

impl Default for RespOptions {
    fn default() -> Self {
        Self {
            layer_scales: vec![1.4, 1.6, 1.8, 2.0],
            pt_density: 6.0,
            hyper_a1: 0.0005,
            hyper_a2: 0.001,
            hyper_b: 0.1,
            max_iter: 300,
            conv: 1.0e-6,
            two_stage: true,
            esp_type: EspType::Total,
        }
    }
}

/// RESP 计算结果。
#[derive(Clone, Debug)]
pub struct RespResult {
    /// 最终 RESP 电荷（单位 e）
    pub charges: Vec<f64>,
    /// 第一阶段电荷（未做两阶段时为最终电荷）
    pub stage1: Vec<f64>,
    /// MK 拟合点数量
    pub n_points: usize,
    /// 拟合静电势的 RMSE（原子单位）
    pub rmse: f64,
    /// 拟合静电势的 RRMSE
    pub rrmse: f64,
    /// 第一阶段迭代次数
    pub stage1_iter: usize,
    /// 第二阶段迭代次数（未执行时为 0）
    pub stage2_iter: usize,
    /// 是否执行了第二阶段
    pub stage2: bool,
    /// 两次拟合是否都收敛
    pub converged: bool,
}

/// RESP 电荷计算器。
pub struct Calculator<'m> {
    /// 分子引用
    pub mole: &'m Mole,
    /// 参与拟合的原子索引列表（默认为全部分子原子）
    pub atms: &'m Vec<usize>,
}

impl<'m> Calculator<'m> {
    /// 创建 RESP 计算器。
    pub fn new(mole: &'m Mole) -> Self {
        let atms = mole.atoms().idxs();
        Self { mole, atms }
    }

    /// 使用默认参数计算标准两阶段 RESP 电荷。
    pub fn resp(&self) -> Result<RespResult, String> {
        self.resp_with(&RespOptions::default())
    }

    /// 使用指定参数计算 RESP 电荷。
    pub fn resp_with(&self, opts: &RespOptions) -> Result<RespResult, String> {
        if opts.layer_scales.is_empty() {
            return Err("RESP: layer_scales 不能为空".to_string());
        }
        if opts.layer_scales.iter().any(|&s| !(s.is_finite() && s > 0.0)) {
            return Err("RESP: layer_scales 必须全部为正数".to_string());
        }
        if !(opts.pt_density.is_finite() && opts.pt_density > 0.0) {
            return Err("RESP: pt_density 必须为正数".to_string());
        }

        let atoms = self.mole.atoms();
        let natm = atoms.len();
        if natm == 0 {
            return Err("RESP: 分子中没有原子".to_string());
        }
        if self.mole.coefs.is_none() || self.mole.basis.is_none() {
            return Err("RESP: 输入文件不包含波函数信息，无法计算静电势".to_string());
        }

        // 1. 生成 MK 拟合点并计算静电势
        let pts = self.mk_points(opts);
        if pts.is_empty() {
            return Err("RESP: 未能生成任何 MK 拟合点".to_string());
        }
        log::info!("RESP: MK 拟合点数量 {}", pts.len());
        let vals = self.esp_values(&pts, opts.esp_type);

        // 2. 构造法方程（与迭代无关，两个阶段共用）
        let centers: Vec<[f64; 3]> = atoms.xyzs().clone();
        let (a_base, b_base) = normal_equations(&centers, &pts, &vals);

        let is_h: Vec<bool> = (0..natm).map(|i| atoms.get(i).nuc == 1).collect();
        let total_rhs = match opts.esp_type {
            EspType::Total => self.mole.charge() as f64,
            EspType::Electronic => {
                let [nela, nelb] = atoms.neles;
                -((nela + nelb) as f64)
            }
        };

        // 3. 第一阶段：弱限制，所有原子均为独立拟合中心
        let groups1: Vec<Vec<usize>> = (0..natm).map(|i| vec![i]).collect();
        let (stage1, iter1, delta1, ok1) = fit_stage(&a_base, &b_base, total_rhs, &groups1, &[], opts.hyper_a1, opts.hyper_b, &is_h, opts.max_iter, opts.conv)?;
        log::info!("RESP 第一阶段：迭代 {} 次，最大电荷变化 {:.3e}{}", iter1, delta1, if ok1 { "" } else { "（未收敛）" });

        // 4. 第二阶段：固定非 sp³/亚甲基碳及其氢，拟合 CH 部分
        let (mut charges, mut stage2, mut iter2, mut converged) = (stage1.clone(), false, 0, ok1);
        if opts.two_stage {
            let (ch_list, groups2) = identify_ch_groups(self.mole);
            if !ch_list.is_empty() {
                let active: Vec<bool> = {
                    let mut m = vec![false; natm];
                    for &i in &ch_list {
                        m[i] = true;
                    }
                    m
                };
                let constraints: Vec<(Vec<usize>, f64)> = (0..natm).filter(|&i| !active[i]).map(|i| (vec![i], stage1[i])).collect();

                let (stage2_chg, iter2_, delta2, ok2) = fit_stage(&a_base, &b_base, total_rhs, &groups2, &constraints, opts.hyper_a2, opts.hyper_b, &is_h, opts.max_iter, opts.conv)?;
                log::info!("RESP 第二阶段：迭代 {} 次，最大电荷变化 {:.3e}{}", iter2_, delta2, if ok2 { "" } else { "（未收敛）" });
                charges = stage2_chg;
                stage2 = true;
                iter2 = iter2_;
                converged = ok1 && ok2;
            }
        }

        // 5. 拟合误差
        let (rmse, rrmse) = fit_errors(&pts, &vals, &centers, &charges);

        Ok(RespResult {
            charges,
            stage1,
            n_points: pts.len(),
            rmse,
            rrmse,
            stage1_iter: iter1,
            stage2_iter: iter2,
            stage2,
            converged,
        })
    }

    /// 生成 Multiwfn 风格的 MK 拟合点。
    pub fn mk_points(&self, opts: &RespOptions) -> Vec<[f64; 3]> {
        let atoms = self.mole.atoms();
        let natm = atoms.len();
        let radii: Vec<f64> = (0..natm).map(|i| resp_vdw_radius(atoms.get(i).nuc)).collect();
        let cut_scale = opts.layer_scales.iter().cloned().fold(f64::INFINITY, f64::min);

        let mut pts = Vec::new();
        for i in 0..natm {
            let center = atoms.get(i).xyz;
            for &scale in &opts.layer_scales {
                let radius = radii[i] * scale;
                // 球面面积（Å²）× 点密度（点/Å²）
                let radius_ang = radius / BOHR;
                let n_expect = (4.0 * PI * radius_ang * radius_ang * opts.pt_density).round() as usize;
                if n_expect == 0 {
                    continue;
                }
                for dir in unit_sphere_points(n_expect) {
                    let p = [center[0] + dir[0] * radius, center[1] + dir[1] * radius, center[2] + dir[2] * radius];
                    // 剔除落在其它原子最内层 vdW 球内的点
                    let mut keep = true;
                    for j in 0..natm {
                        if j == i {
                            continue;
                        }
                        let cut = radii[j] * cut_scale;
                        if dist2(&p, &atoms.get(j).xyz) < cut * cut {
                            keep = false;
                            break;
                        }
                    }
                    if keep {
                        pts.push(p);
                    }
                }
            }
        }
        pts
    }

    /// 计算一组点上的静电势。
    pub fn esp_values(&self, pts: &[[f64; 3]], esp_type: EspType) -> Vec<f64> {
        let caler = density::Calculator::new(self.mole);
        // 预热解析积分器，避免并行线程同时初始化
        if let Some(p0) = pts.first() {
            let _ = match esp_type {
                EspType::Total => caler.tot_esp(p0),
                EspType::Electronic => caler.ele_esp(p0),
            };
        }
        pts.par_iter()
            .map(|p| match esp_type {
                EspType::Total => caler.tot_esp(p),
                EspType::Electronic => caler.ele_esp(p),
            })
            .collect()
    }
}

/// 构造法方程 `A q = B`，其中
/// `A_ij = Σ_pt 1/(r_i r_j)`，`B_i = Σ_pt V_pt / r_i`。
fn normal_equations(centers: &[[f64; 3]], pts: &[[f64; 3]], vals: &[f64]) -> (Array2<f64>, Vec<f64>) {
    let n = centers.len();
    // 并行计算每个拟合点到所有中心的倒数距离
    let inv_r: Vec<Vec<f64>> = pts.par_iter().map(|p| centers.iter().map(|c| 1.0 / dist(p, c)).collect()).collect();

    let mut a = Array2::<f64>::zeros((n, n));
    let mut b = vec![0.0; n];
    for (ipt, inv) in inv_r.iter().enumerate() {
        let v = vals[ipt];
        for i in 0..n {
            b[i] += v * inv[i];
            for j in i..n {
                a[(i, j)] += inv[i] * inv[j];
            }
        }
    }
    for i in 0..n {
        for j in 0..i {
            a[(i, j)] = a[(j, i)];
        }
    }
    (a, b)
}

/// 执行一个阶段的带约束 RESP 拟合。
///
/// - `groups`：等价原子分组，同一组内的原子共享电荷，且必须是所有原子的一个划分。
/// - `constraints`：额外电荷约束 `(原子集合, 目标电荷总和)`。
/// - `a_penalty`：双曲限制强度；为 0 时只求解一次。
///
/// 返回 `(原子电荷, 迭代次数, 最后一次最大电荷变化, 是否收敛)`。
#[allow(clippy::too_many_arguments)]
fn fit_stage(
    a_base: &Array2<f64>,
    b_base: &[f64],
    total_charge: f64,
    groups: &[Vec<usize>],
    constraints: &[(Vec<usize>, f64)],
    a_penalty: f64,
    b_penalty: f64,
    is_h: &[bool],
    max_iter: usize,
    conv: f64,
) -> Result<(Vec<f64>, usize, f64, bool), String> {
    let natm = a_base.nrows();
    let ng = groups.len();
    let nc = constraints.len();
    let dim = ng + nc + 1;
    let total_idx = ng;
    let cons_start = ng + 1;

    // 原子 -> 等价组映射，并检查分组是否是一个完整划分
    let mut group_of = vec![0usize; natm];
    let mut group_sizes = vec![0usize; ng];
    let mut seen = vec![0usize; natm];
    for (g, members) in groups.iter().enumerate() {
        for &i in members {
            if i >= natm {
                return Err(format!("RESP: 等价约束中的原子索引越界 ({i})"));
            }
            group_of[i] = g;
            group_sizes[g] += 1;
            seen[i] += 1;
        }
    }
    if seen.iter().any(|&c| c != 1) {
        return Err("RESP: 等价原子分组未覆盖每个原子恰好一次".to_string());
    }

    let mut q = vec![0.0; natm];
    let mut iter_used = 0;
    let mut last_delta = 0.0;
    let mut converged = false;

    for iter in 1..=max_iter.max(1) {
        let q_old = q.clone();
        let mut amat = Array2::<f64>::zeros((dim, dim));
        let mut bvec = vec![0.0; dim];

        // 原子-原子块，同时加入双曲限制的对角项
        for i in 0..natm {
            let gi = group_of[i];
            bvec[gi] += b_base[i];
            for j in 0..natm {
                let gj = group_of[j];
                let mut v = a_base[(i, j)];
                if i == j && !is_h[i] {
                    v += a_penalty / (q_old[i] * q_old[i] + b_penalty * b_penalty).sqrt();
                }
                amat[(gi, gj)] += v;
            }
        }

        // 总电荷约束
        for g in 0..ng {
            amat[(g, total_idx)] = group_sizes[g] as f64;
            amat[(total_idx, g)] = group_sizes[g] as f64;
        }
        bvec[total_idx] = total_charge;

        // 额外电荷约束
        for (k, (members, val)) in constraints.iter().enumerate() {
            let ci = cons_start + k;
            for &i in members {
                let g = group_of[i];
                amat[(g, ci)] += 1.0;
                amat[(ci, g)] += 1.0;
            }
            bvec[ci] = *val;
        }

        let q_group = solve_dense(&amat, &bvec)?;
        for i in 0..natm {
            q[i] = q_group[group_of[i]];
        }

        last_delta = q.iter().zip(q_old.iter()).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        iter_used = iter;

        if a_penalty == 0.0 || last_delta < conv {
            converged = true;
            break;
        }
    }

    Ok((q, iter_used, last_delta, converged))
}

/// 用 LU 分解求解线性方程组。
fn solve_dense(a: &Array2<f64>, b: &[f64]) -> Result<Vec<f64>, String> {
    if a.iter().any(|v| !v.is_finite()) || b.iter().any(|v| !v.is_finite()) {
        return Err("RESP: 法方程包含非有限值".to_string());
    }
    let inv = maths::inv(a);
    if inv.iter().any(|v| !v.is_finite()) {
        return Err("RESP: 法方程奇异，矩阵求逆失败".to_string());
    }
    let rhs = Array1::from_vec(b.to_vec());
    let q = inv.dot(&rhs).to_vec();
    if q.iter().any(|v| !v.is_finite()) {
        return Err("RESP: 线性方程组求解失败".to_string());
    }
    Ok(q)
}

/// 识别标准 RESP 第二阶段需要继续拟合的 sp³ 碳/亚甲基碳及其所连氢。
///
/// 返回 `(活跃原子列表, 第二阶段等价分组)`。等价分组为：
/// 每个活跃碳上的一组等价氢 + 其余原子的单原子组。
fn identify_ch_groups(mole: &Mole) -> (Vec<usize>, Vec<Vec<usize>>) {
    let atoms = mole.atoms();
    let natm = atoms.len();

    let mut ch_list = Vec::new();
    let mut h_groups: Vec<Vec<usize>> = Vec::new();

    for i in 0..natm {
        if atoms.get(i).nuc != 6 {
            continue;
        }
        let nebs = atoms.nebs(i);
        let nbond = nebs.len();
        let mut h_idx = Vec::new();
        let mut n_c = 0;
        for &j in &nebs {
            let nuc = atoms.get(j).nuc;
            if nuc == 1 {
                h_idx.push(j);
            } else if nuc == 6 {
                n_c += 1;
            }
        }
        // sp³ 碳（4 配位）或 =CH₂ / -CH₂-（2 个 H 且只连 1 个 C）
        if nbond == 4 || (h_idx.len() == 2 && n_c == 1) {
            ch_list.push(i);
            if !h_idx.is_empty() {
                for &h in &h_idx {
                    ch_list.push(h);
                }
                h_groups.push(h_idx);
            }
        }
    }

    // 构建完整划分：活跃碳上的等价氢组 + 其余原子的单原子组
    let mut covered = vec![false; natm];
    for g in &h_groups {
        for &i in g {
            covered[i] = true;
        }
    }
    let mut groups = h_groups;
    for i in 0..natm {
        if !covered[i] {
            groups.push(vec![i]);
        }
    }

    (ch_list, groups)
}

/// 计算 RESP 拟合的 RMSE / RRMSE。
fn fit_errors(pts: &[[f64; 3]], vals: &[f64], centers: &[[f64; 3]], charges: &[f64]) -> (f64, f64) {
    let ncen = centers.len();
    let mut sse = 0.0;
    let mut ssq = 0.0;
    for (p, &v) in pts.iter().zip(vals.iter()) {
        let mut v_fit = 0.0;
        for i in 0..ncen {
            v_fit += charges[i] / dist(p, &centers[i]);
        }
        sse += (v - v_fit).powi(2);
        ssq += v * v;
    }
    let n = pts.len() as f64;
    let rmse = (sse / n).sqrt();
    let rrmse = if ssq > 0.0 { (sse / ssq).sqrt() } else { f64::NAN };
    (rmse, rrmse)
}

/// Multiwfn `unitspherept`：在单位球面上按经纬度均匀撒点。
///
/// 输入期望点数，返回实际生成的（不超过期望值的）点。
fn unit_sphere_points(numpt: usize) -> Vec<[f64; 3]> {
    if numpt == 0 {
        return vec![];
    }
    let pi = PI;
    let numequ = (numpt as f64 * pi).sqrt() as usize;
    let numvert = numequ / 2;
    if numvert == 0 {
        return vec![[0.0, 0.0, 1.0]];
    }

    let mut pts = Vec::with_capacity(numpt);
    for ivert in 0..=numvert {
        let angz = ivert as f64 / numvert as f64 * pi;
        let scalexy = angz.sin();
        let z = angz.cos();
        let numxy = ((numequ as f64) * scalexy) as usize;
        let numxy = numxy.max(1);
        for ihori in 1..=numxy {
            if pts.len() >= numpt {
                return pts;
            }
            let angxy = 2.0 * pi * ihori as f64 / numxy as f64;
            pts.push([angxy.cos() * scalexy, angxy.sin() * scalexy, z]);
        }
    }
    pts
}

/// RESP/MK 拟合使用的 vdW 半径（Bohr）。
///
/// 对 H–Cl 使用 Multiwfn 的 MK 专用半径；其它元素使用 UFF 半径除以 1.2。
fn resp_vdw_radius(nuc: usize) -> f64 {
    // H–Cl，单位 Å
    const MK_RADII: [f64; 17] = [1.20, 1.20, 1.37, 1.45, 1.45, 1.50, 1.50, 1.40, 1.35, 1.30, 1.57, 1.65, 1.65, 1.80, 1.80, 1.75, 1.70];
    let ang = if (1..=17).contains(&nuc) { MK_RADII[nuc - 1] } else { uff_radius(nuc) / 1.2 };
    ang * BOHR
}

/// UFF vdW 半径（Å），取自 JACS 114, 10024 (1992) Table 1（即 Multiwfn 的 `vdwr_UFF`）。
fn uff_radius(nuc: usize) -> f64 {
    const UFF: [f64; 104] = [
        0.4, // 0 Ghost
        1.443, 1.181, // H, He
        1.2255, 1.3725, 2.0415, 1.9255, 1.83, 1.75, 1.682, 1.6215, // Li–Ne
        1.4915, 1.5105, 2.2495, 2.1475, 2.0735, 2.0175, 1.9735, 1.934, // Na–Ar
        1.906, 1.6995, 1.6475, 1.5875, 1.572, 1.5115, 1.4805, 1.456, 1.436, 1.417, 1.7475, 1.3815, 2.1915, 2.14, 2.115, 2.1025, 2.0945, 2.0705, // K–Kr
        2.057, 1.8205, 1.6725, 1.562, 1.5825, 1.526, 1.499, 1.4815, 1.4645, 1.4495, 1.574, 1.424, 2.2315, 2.196, 2.21, 2.235, 2.25, 2.202, // Rb–Xe
        2.2585, 1.8515, 1.761, 1.778, 1.803, 1.7875, 1.7735, 1.76, 1.7465, 1.684, 1.7255, 1.714, 1.7045, 1.6955, 1.687, 1.6775, 1.82, // Cs–Lu
        1.5705, 1.585, 1.5345, 1.477, 1.56, 1.42, 1.377, 1.6465, 1.3525, 2.1735, 2.1485, 2.185, 2.3545, 2.375, 2.3825, // Hf–Rn
        2.45, 1.8385, 1.739, 1.698, 1.712, 1.6975, 1.712, 1.712, 1.6905, 1.663, 1.6695, 1.6565, 1.6495, 1.643, 1.637, 1.624, 1.618, // Fr–Lr
    ];
    if nuc < UFF.len() { UFF[nuc] } else { 2.0 }
}

#[inline]
fn dist2(p1: &[f64; 3], p2: &[f64; 3]) -> f64 {
    (p1[0] - p2[0]).powi(2) + (p1[1] - p2[1]).powi(2) + (p1[2] - p2[2]).powi(2)
}

#[inline]
fn dist(p1: &[f64; 3], p2: &[f64; 3]) -> f64 {
    dist2(p1, p2).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_unit_sphere_points() {
        for &n in &[6usize, 50, 110, 230] {
            let pts = unit_sphere_points(n);
            assert!(!pts.is_empty());
            assert!(pts.len() <= n);
            for p in pts {
                let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
                assert!((r - 1.0).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_uff_table() {
        assert!((uff_radius(6) - 1.9255).abs() < 1e-12);
        assert!((uff_radius(1) - 1.443).abs() < 1e-12);
        assert!((uff_radius(200) - 2.0).abs() < 1e-12);
    }
}
