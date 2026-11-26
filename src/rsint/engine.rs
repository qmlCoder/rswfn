#![allow(non_snake_case)]
//! 基于壳层块的电子积分引擎（SCF 的双电子部分）
//!
//! 旧的 `Mat2e` 把 (ij|kl) 压缩成一个 nato⁴/8 的一维数组，构造 J/K 时按
//! 原子轨道对随机访存。这个设计有三个问题：
//!
//! 1. 内存是 O(N⁴)。310 个基函数就要 9.3 GB，460 个基函数要 45 GB，
//!    后者在 64 GB 机器上直接把进程拖死；
//! 2. 构造 J/K 的内层循环每次都要做两次整数除法和一次跨上百 MB 的随机
//!    读取，实测占 SCF 的主要时间；
//! 3. 积分按 8 重对称压缩存储，写入时又是随机散射，构造积分本身也慢。
//!
//! 本模块改为按**壳层四元组**存储：每个块是 `(di,dj,dk,dl)` 的一段连续
//! 内存，内容是 libcint 原始顺序、未归一化的积分值。由此得到：
//!
//! - 内存 = 被 Schwarz 筛选保留下来的块，通常比 N⁴/8 小一个数量级以上；
//! - 构造 J/K 时块内是连续访存，J/K/D 三个矩阵都能驻留 L2，做分块稠密
//!   累加，不再有随机访存；
//! - 归一化因子被吸收进密度矩阵（`D' = f_a f_b D_ab`），逐元乘法消失。
//!
//! 当筛选后的块仍超过内存上限时，自动退化为**直接 SCF**：不缓存积分，
//! 每轮迭代重新计算通过 Schwarz + 密度筛选的块，内存回到 O(N²)。

use std::os::raw::{c_double, c_int};
use std::sync::Arc;
use std::time::{Duration, Instant};

use ::rsint::src::cint2e::{CINTOpt, cint2e_cart, cint2e_cart_optimizer, cint2e_sph, cint2e_sph_optimizer};
use ::rsint::src::optimizer::CINTdel_optimizer;
use ndarray::Array2;
use rayon::prelude::*;

use crate::base::Mole;

use super::Mat2e;

/// libcint 优化器指针的线程安全包装
///
/// `CINTOpt` 只在积分计算时被只读访问（优化器是纯查表缓存），因此多线程
/// 共享是安全的；libcint 的 C ABI 里没有 `const`，只能这样包一层。
struct OptPtr(*mut CINTOpt);
unsafe impl Send for OptPtr {}
unsafe impl Sync for OptPtr {}

/// 裸指针的线程安全包装，用于并行填充互不重叠的积分块
struct SendPtr(*mut f64);
impl SendPtr {
    #[inline]
    fn get(&self) -> *mut f64 {
        self.0
    }
}
unsafe impl Send for SendPtr {}
unsafe impl Sync for SendPtr {}

/// 壳层对（i ≥ j）及其 Schwarz 上界
#[derive(Clone, Copy, Debug)]
struct Pair {
    i: u16,
    j: u16,
    /// sqrt(max|(ij|ij)|)，用于 |(ij|kl)| ≤ q_ij · q_kl 的筛选
    q: f64,
}

/// 筛选后保留的积分块（连续存储）
struct Store {
    /// 左侧壳层对索引
    ps: Vec<u32>,
    /// 右侧壳层对索引
    qs: Vec<u32>,
    /// 每块在一维数据中的起点
    offs: Vec<u64>,
    /// 全部积分值（libcint 原始顺序）
    data: Vec<f64>,
}

impl Store {
    fn len(&self) -> usize {
        self.ps.len()
    }
    fn bytes(&self) -> usize {
        self.data.len() * 8 + self.offs.len() * 8 + self.ps.len() * 8
    }
}

/// 局部旋转/裁剪时复用的暂存缓冲（避免每个壳层块都重新分配）
#[derive(Default)]
struct PiScratch {
    g: Vec<f64>,
    a1: Vec<f64>,
    a2: Vec<f64>,
    a3: Vec<f64>,
    l: Vec<f64>,
}

impl PiScratch {
    fn reserve(&mut self, n: usize) {
        if self.g.len() < n {
            self.g.resize(n, 0.0);
            self.a1.resize(n, 0.0);
            self.a2.resize(n, 0.0);
            self.a3.resize(n, 0.0);
            self.l.resize(n, 0.0);
        }
    }
}

/// 壳层块积分引擎
pub struct Engine<'m> {
    mole: &'m Mole,
    form: String,
    /// 原子轨道数
    nato: usize,
    /// 壳层数
    nshl: usize,
    /// 每个壳层的 AO 起点（两种顺序下相同）
    offs: Vec<usize>,
    /// 每个壳层的 AO 个数
    dims: Vec<usize>,
    /// libcint 序号 -> 高斯序号
    gidx: Vec<usize>,
    /// libcint 序号下的归一化因子
    fix: Vec<f64>,
    natm: i32,
    nbas: i32,
    atm: Arc<Vec<c_int>>,
    bas: Arc<Vec<c_int>>,
    env: Arc<Vec<c_double>>,
    opt: OptPtr,
    /// 壳层对表
    pairs: Vec<Pair>,
    /// 保留的积分块；None 表示直接 SCF（不缓存积分）
    store: Option<Store>,
    /// Schwarz 筛选阈值
    pub screen: f64,
    /// 单块最大元素数（直接模式下的线程局部缓冲）
    max_block: usize,
    /// 建库耗时
    pub time_build: Duration,
    /// 保留的块数
    pub nblock: usize,
    /// 积分数据占用字节数
    pub mem_bytes: usize,
}

impl Drop for Engine<'_> {
    fn drop(&mut self) {
        unsafe {
            let mut opt = self.opt.0;
            CINTdel_optimizer(&mut opt as *mut *mut CINTOpt as *mut *mut ::rsint::src::optimizer::CINTOpt);
        }
    }
}

unsafe impl Send for Engine<'_> {}
unsafe impl Sync for Engine<'_> {}

impl<'m> Engine<'m> {
    /// 创建引擎并预计算筛选后的积分块
    ///
    /// - `screen`: Schwarz 筛选阈值（|(ij|kl)| 的上界小于它时丢弃该块）
    /// - `max_mem`: 积分块允许占用的最大字节数；超过则退化为直接 SCF
    pub fn new(mole: &'m Mole, form: &str, screen: f64, max_mem: usize) -> Self {
        let t0 = Instant::now();
        let raw = super::rsint::Rsint::new(mole);
        let (nato, fix) = match form {
            "sph" => (raw.nato[0], raw.fixs[0].clone()),
            "car" => (raw.nato[1], raw.fixs[1].clone()),
            _ => panic!("未知积分形式: {form}"),
        };
        let nums: [usize; 7] = if form == "car" { [1, 3, 6, 10, 15, 21, 28] } else { [1, 3, 5, 7, 9, 11, 13] };
        let nshl = raw.bas.len() / 8;
        let mut offs = vec![0usize; nshl];
        let mut dims = vec![0usize; nshl];
        let mut loc = 0usize;
        for s in 0..nshl {
            let ang = raw.bas[s * 8 + 1] as usize;
            offs[s] = loc;
            dims[s] = nums[ang];
            loc += nums[ang];
        }
        assert_eq!(loc, nato, "壳层维度累加与 nato 不一致");

        // libcint AO 序号 -> 高斯 AO 序号
        let mut gidx = vec![0usize; nato];
        for s in 0..nshl {
            let perm = shell_perm(dims[s]); // 高斯位置 -> libcint 位置
            for (g, lc) in perm.iter().enumerate() {
                gidx[offs[s] + *lc] = offs[s] + g;
            }
        }

        let natm = raw.atm.len() / 6;
        let nbas = raw.bas.len() / 8;
        let atm = Arc::new(raw.atm.iter().map(|&v| v as c_int).collect::<Vec<_>>());
        let bas = Arc::new(raw.bas.iter().map(|&v| v as c_int).collect::<Vec<_>>());
        let env = Arc::new(raw.env.clone());

        let mut opt: *mut CINTOpt = std::ptr::null_mut();
        unsafe {
            match form {
                "sph" => cint2e_sph_optimizer(&mut opt, atm.as_ptr() as *mut c_int, natm as c_int, bas.as_ptr() as *mut c_int, nbas as c_int, env.as_ptr() as *mut c_double),
                _ => cint2e_cart_optimizer(&mut opt, atm.as_ptr() as *mut c_int, natm as c_int, bas.as_ptr() as *mut c_int, nbas as c_int, env.as_ptr() as *mut c_double),
            }
        }

        let mut this = Self {
            mole,
            form: form.to_string(),
            nato,
            nshl,
            offs,
            dims,
            gidx,
            fix,
            natm: natm as i32,
            nbas: nbas as i32,
            atm,
            bas,
            env,
            opt: OptPtr(opt),
            pairs: Vec::new(),
            store: None,
            screen,
            max_block: 1,
            time_build: Duration::ZERO,
            nblock: 0,
            mem_bytes: 0,
        };
        this.build_pairs();
        this.build_store(max_mem);
        this.time_build = t0.elapsed();
        this
    }

    /// 自检：把存储的积分块与重新计算的结果逐块比较（调试用）
    pub fn debug_check_blocks(&self) -> (usize, usize, f64) {
        let Some(store) = &self.store else {
            return (0, 0, 0.0);
        };
        let mut worst = 0.0f64;
        let mut worst_bi = 0usize;
        for bi in 0..store.len() {
            let p = store.ps[bi] as usize;
            let q = store.qs[bi] as usize;
            let (si, sj) = self.pair_shells(p);
            let (sk, sl) = self.pair_shells(q);
            let nels = self.dims[si] * self.dims[sj] * self.dims[sk] * self.dims[sl];
            let off = store.offs[bi] as usize;
            let mut buf = vec![0.0f64; nels];
            let mut shls = [si as c_int, sj as c_int, sk as c_int, sl as c_int];
            self.call_cint(&mut buf, &mut shls);
            for i in 0..nels {
                let d = (buf[i] - store.data[off + i]).abs();
                if d > worst {
                    worst = d;
                    worst_bi = bi;
                }
            }
        }
        (store.len(), worst_bi, worst)
    }

    /// 原子轨道数
    pub fn nato(&self) -> usize {
        self.nato
    }

    /// 关联的分子
    pub fn mole(&self) -> &'m Mole {
        self.mole
    }

    /// 是否使用直接 SCF（不缓存积分）
    pub fn is_direct(&self) -> bool {
        self.store.is_none()
    }

    /// 构建局部坐标系下的 π 子块双电子积分（供 π-only 的 LKS 使用）
    ///
    /// `tmats[s]` 是壳层 `s` 的局部坐标变换（局部 = `T ·` 全局），
    /// `pi_comp[s]` 是该壳层内部的 π 分量：`(局部分量下标, π 约化基下标)`。
    ///
    /// 直接复用已筛选的壳层块，不再做全维 4 指标变换：对每个壳层四元组块，
    /// 先用归一化因子和 libcint→高斯置换得到高斯顺序的块，再沿 4 个方向做小
    /// 矩阵乘 `(Ti ⊗ Tj) G (Tk ⊗ Tl)ᵀ`，最后只取 π 分量写入 `nπ` 维的 [`Mat2e`]。
    /// 任何一层壳层没有 π 分量时整块跳过（s 壳层占多数，因此大部分块不进变换）。
    pub fn mat2e_local_pi(&self, tmats: &[Array2<f64>], pi_comp: &[Vec<(usize, usize)>], npi: usize) -> Mat2e {
        let mut emat = Mat2e::new(npi);
        let mut scratch = PiScratch::default();
        if let Some(store) = &self.store {
            for bi in 0..store.len() {
                let p = store.ps[bi] as usize;
                let q = store.qs[bi] as usize;
                let (si, sj) = self.pair_shells(p);
                let (sk, sl) = self.pair_shells(q);
                let off = store.offs[bi] as usize;
                let nels = self.dims[si] * self.dims[sj] * self.dims[sk] * self.dims[sl];
                self.pour_pi_block(&store.data[off..off + nels], si, sj, sk, sl, tmats, pi_comp, &mut emat, &mut scratch);
            }
        } else {
            // 直接模式：没有缓存，按 Schwarz 筛选临时重算每个块
            let npair = self.pairs.len();
            let mut buf = vec![0.0f64; self.max_block.max(1)];
            for p in 0..npair {
                let (si, sj) = self.pair_shells(p);
                let qp = self.pairs[p].q;
                for q in 0..=p {
                    if qp * self.pairs[q].q < self.screen {
                        continue;
                    }
                    let (sk, sl) = self.pair_shells(q);
                    let nels = self.dims[si] * self.dims[sj] * self.dims[sk] * self.dims[sl];
                    let mut shls = [si as c_int, sj as c_int, sk as c_int, sl as c_int];
                    self.call_cint(&mut buf[..nels], &mut shls);
                    self.pour_pi_block(&buf[..nels], si, sj, sk, sl, tmats, pi_comp, &mut emat, &mut scratch);
                }
            }
        }
        emat
    }

    /// 把一个壳层块旋到局部坐标系并只保留 π 分量，累加进 π 约化基的双电子积分
    #[allow(clippy::too_many_arguments)]
    fn pour_pi_block(&self, buf: &[f64], si: usize, sj: usize, sk: usize, sl: usize, tmats: &[Array2<f64>], pi_comp: &[Vec<(usize, usize)>], emat: &mut Mat2e, scratch: &mut PiScratch) {
        let (ci, cj, ck, cl) = (&pi_comp[si], &pi_comp[sj], &pi_comp[sk], &pi_comp[sl]);
        if ci.is_empty() || cj.is_empty() || ck.is_empty() || cl.is_empty() {
            return;
        }
        let (di, dj, dk, dl) = (self.dims[si], self.dims[sj], self.dims[sk], self.dims[sl]);
        let (oi, oj, ok, ol) = (self.offs[si], self.offs[sj], self.offs[sk], self.offs[sl]);
        let nels = di * dj * dk * dl;
        scratch.reserve(nels);
        // 1) 高斯顺序 + 归一化因子
        for a in 0..di {
            let fa = self.fix[oi + a];
            let ga = self.gidx[oi + a] - oi;
            for b in 0..dj {
                let fb = fa * self.fix[oj + b];
                let gb = self.gidx[oj + b] - oj;
                for c in 0..dk {
                    let fc = fb * self.fix[ok + c];
                    let gc = self.gidx[ok + c] - ok;
                    for d in 0..dl {
                        let raw = buf[a + di * (b + dj * (c + dk * d))];
                        let gd = self.gidx[ol + d] - ol;
                        let idx = ga + di * (gb + dj * (gc + dk * gd));
                        scratch.g[idx] = raw * fc * self.fix[ol + d];
                    }
                }
            }
        }
        // 2) 沿 4 个方向依次做局部坐标变换：L = (Ti ⊗ Tj) G (Tk ⊗ Tl)ᵀ
        let t0 = &tmats[si];
        for i in 0..di {
            for b in 0..dj {
                for c in 0..dk {
                    for d in 0..dl {
                        let mut s = 0.0;
                        for aa in 0..di {
                            s += t0[(i, aa)] * scratch.g[aa + di * (b + dj * (c + dk * d))];
                        }
                        scratch.a1[i + di * (b + dj * (c + dk * d))] = s;
                    }
                }
            }
        }
        let t1 = &tmats[sj];
        for i in 0..di {
            for j in 0..dj {
                for c in 0..dk {
                    for d in 0..dl {
                        let mut s = 0.0;
                        for bb in 0..dj {
                            s += t1[(j, bb)] * scratch.a1[i + di * (bb + dj * (c + dk * d))];
                        }
                        scratch.a2[i + di * (j + dj * (c + dk * d))] = s;
                    }
                }
            }
        }
        let t2 = &tmats[sk];
        for i in 0..di {
            for j in 0..dj {
                for k in 0..dk {
                    for d in 0..dl {
                        let mut s = 0.0;
                        for cc in 0..dk {
                            s += t2[(k, cc)] * scratch.a2[i + di * (j + dj * (cc + dk * d))];
                        }
                        scratch.a3[i + di * (j + dj * (k + dk * d))] = s;
                    }
                }
            }
        }
        let t3 = &tmats[sl];
        for i in 0..di {
            for j in 0..dj {
                for k in 0..dk {
                    for ll in 0..dl {
                        let mut s = 0.0;
                        for dd in 0..dl {
                            s += t3[(ll, dd)] * scratch.a3[i + di * (j + dj * (k + dk * dd))];
                        }
                        scratch.l[i + di * (j + dj * (k + dk * ll))] = s;
                    }
                }
            }
        }
        // 3) 只取 π 分量写入 π 约化基的双电子积分
        for &(li, pi) in ci {
            for &(lj, pj) in cj {
                for &(lk, pk) in ck {
                    for &(ll, pl) in cl {
                        emat.set(pi, pj, pk, pl, scratch.l[li + di * (lj + dj * (lk + dk * ll))]);
                    }
                }
            }
        }
    }

    /// 计算两个壳层之间的 Schwarz 上界
    fn pair_bound(&self, i: usize, j: usize) -> f64 {
        let (di, dj) = (self.dims[i], self.dims[j]);
        let mut buf = vec![0.0f64; di * dj * di * dj];
        let mut shls = [i as c_int, j as c_int, i as c_int, j as c_int];
        self.call_cint(&mut buf, &mut shls);
        let mut mx = 0.0f64;
        for v in buf.iter() {
            let a = v.abs();
            if a > mx {
                mx = a;
            }
        }
        mx.sqrt()
    }

    /// 调用 libcint 计算一个壳层四元组（原始顺序、未归一化）
    #[inline]
    fn call_cint(&self, buf: &mut [f64], shls: &mut [c_int; 4]) {
        let ret = unsafe {
            match self.form.as_str() {
                "sph" => cint2e_sph(
                    buf.as_mut_ptr(),
                    shls.as_mut_ptr(),
                    self.atm.as_ptr() as *mut c_int,
                    self.natm,
                    self.bas.as_ptr() as *mut c_int,
                    self.nbas,
                    self.env.as_ptr() as *mut c_double,
                    self.opt.0,
                ),
                _ => cint2e_cart(
                    buf.as_mut_ptr(),
                    shls.as_mut_ptr(),
                    self.atm.as_ptr() as *mut c_int,
                    self.natm,
                    self.bas.as_ptr() as *mut c_int,
                    self.nbas,
                    self.env.as_ptr() as *mut c_double,
                    self.opt.0,
                ),
            }
        };
        // rsint 的 cint2e 返回 `(empty == 0)`：0 表示该块被优化器判为全零，
        // 1 表示正常算出。两种都不是错误，这里不据此报错。
        let _ = ret;
    }

    /// 建立壳层对表与 Schwarz 上界
    fn build_pairs(&mut self) {
        for i in 0..self.nshl {
            for j in 0..=i {
                let q = self.pair_bound(i, j);
                self.pairs.push(Pair { i: i as u16, j: j as u16, q });
            }
        }
    }

    /// 预计算筛选后的积分块；超过内存上限时放弃缓存，改为直接 SCF
    ///
    /// 第一遍只做筛选与容量统计（很快），随后在壳层四元组上并行调用
    /// libcint 填充数据；各块的写入区间互不重叠，因此可以安全地并行写。
    fn build_store(&mut self, max_mem: usize) {
        let npair = self.pairs.len();
        // 第一遍：并行统计每个壳层对 p 通过 Schwarz 筛选的 q 与块大小。
        // 旧实现是串行 O(npair²)（C18H38 上约 1 s），这里按 p 并行，
        // 块顺序仍按 (p, q) 升序，因此结果与串行完全一致。
        let per_p: Vec<Vec<(u32, u32)>> = (0..npair)
            .into_par_iter()
            .map(|p| {
                let (pi, pj) = self.pair_shells(p);
                let qp = self.pairs[p].q;
                let mut list: Vec<(u32, u32)> = Vec::new();
                for q in 0..=p {
                    if qp * self.pairs[q].q < self.screen {
                        continue;
                    }
                    let (qk, ql) = self.pair_shells(q);
                    let n = self.dims[pi] * self.dims[pj] * self.dims[qk] * self.dims[ql];
                    list.push((q as u32, n as u32));
                }
                list
            })
            .collect();
        let mut ps: Vec<u32> = Vec::new();
        let mut qs: Vec<u32> = Vec::new();
        let mut offs: Vec<u64> = Vec::new();
        let mut sizes: Vec<u32> = Vec::new();
        let mut total = 0usize;
        for (p, list) in per_p.into_iter().enumerate() {
            let (pi, pj) = self.pair_shells(p);
            let blk = self.dims[pi] * self.dims[pj] * self.dims[pi] * self.dims[pj];
            if blk > self.max_block {
                self.max_block = blk;
            }
            for (q, n) in list {
                if n as usize > self.max_block {
                    self.max_block = n as usize;
                }
                ps.push(p as u32);
                qs.push(q);
                offs.push(total as u64);
                sizes.push(n);
                total += n as usize;
            }
        }
        let bytes = total * 8 + offs.len() * 16;
        if bytes > max_mem {
            log::info!("双电子积分块 {:.2} GB 超过上限 {:.2} GB，改用直接 SCF（不缓存积分）", bytes as f64 / 1e9, max_mem as f64 / 1e9);
            self.store = None;
            self.mem_bytes = 0;
            self.nblock = 0;
            return;
        }
        let mut data = vec![0.0f64; total];
        let base = SendPtr(data.as_mut_ptr());
        ps.par_iter().zip(qs.par_iter()).zip(offs.par_iter().zip(sizes.par_iter())).for_each(|((p, q), (off, n))| {
            let (si, sj) = self.pair_shells(*p as usize);
            let (sk, sl) = self.pair_shells(*q as usize);
            let mut shls = [si as c_int, sj as c_int, sk as c_int, sl as c_int];
            let dst = unsafe { std::slice::from_raw_parts_mut(base.get().add(*off as usize), *n as usize) };
            self.call_cint(dst, &mut shls);
        });
        let store = Store { ps, qs, offs, data };
        self.nblock = store.len();
        self.mem_bytes = store.bytes();
        self.store = Some(store);
    }

    #[inline]
    fn pair_shells(&self, p: usize) -> (usize, usize) {
        (self.pairs[p].i as usize, self.pairs[p].j as usize)
    }

    /// 把密度矩阵转到 libcint 序号空间并吸收归一化因子
    fn scaled_density(&self, dmat: &Array2<f64>) -> Vec<f64> {
        let n = self.nato;
        let mut dl = vec![0.0f64; n * n];
        for a in 0..n {
            let ga = self.gidx[a];
            let fa = self.fix[a];
            for b in 0..n {
                dl[a * n + b] = dmat[(ga, self.gidx[b])] * fa * self.fix[b];
            }
        }
        dl
    }

    /// libcint 序号空间的 J/K 转回高斯序号空间并施加归一化因子
    fn unscale(&self, jl: &[f64], out: &mut Array2<f64>) {
        let n = self.nato;
        for a in 0..n {
            let ga = self.gidx[a];
            let mut row = out.row_mut(ga);
            for b in 0..n {
                row[self.gidx[b]] = jl[a * n + b] * self.fix[a] * self.fix[b];
            }
        }
    }

    /// 计算库伦矩阵 J 与交换矩阵 K（闭壳层，单密度）
    pub fn jk(&self, dmat: &Array2<f64>) -> (Array2<f64>, Array2<f64>) {
        let n = self.nato;
        let dl = self.scaled_density(dmat);
        let mut jl = vec![0.0f64; n * n];
        let mut kl = vec![0.0f64; n * n];
        self.accumulate(&dl, &mut jl, &mut kl);
        let mut jmat = Array2::zeros((n, n));
        let mut kmat = Array2::zeros((n, n));
        self.unscale(&jl, &mut jmat);
        self.unscale(&kl, &mut kmat);
        (jmat, kmat)
    }

    /// 计算 J、K_α、K_β（开壳层）
    pub fn jk_spin(&self, d_tot: &Array2<f64>, d_a: &Array2<f64>, d_b: &Array2<f64>) -> (Array2<f64>, Array2<f64>, Array2<f64>) {
        let n = self.nato;
        let dt = self.scaled_density(d_tot);
        let da = self.scaled_density(d_a);
        let db = self.scaled_density(d_b);
        let mut jl = vec![0.0f64; n * n];
        let mut kal = vec![0.0f64; n * n];
        let mut kbl = vec![0.0f64; n * n];
        self.accumulate_spin(&dt, &da, &db, &mut jl, &mut kal, &mut kbl);
        let mut jmat = Array2::zeros((n, n));
        let mut kamat = Array2::zeros((n, n));
        let mut kbmat = Array2::zeros((n, n));
        self.unscale(&jl, &mut jmat);
        self.unscale(&kal, &mut kamat);
        self.unscale(&kbl, &mut kbmat);
        (jmat, kamat, kbmat)
    }

    /// 遍历全部保留的块累加 J、K
    fn accumulate(&self, dl: &[f64], jl: &mut [f64], kl: &mut [f64]) {
        let n = self.nato;
        if let Some(store) = &self.store {
            let (aj, ak) = (0..store.len())
                .into_par_iter()
                .with_min_len(1024)
                .fold(
                    || (vec![0.0f64; n * n], vec![0.0f64; n * n]),
                    |(mut j, mut k), bi| {
                        let p = store.ps[bi] as usize;
                        let q = store.qs[bi] as usize;
                        let (si, sj) = self.pair_shells(p);
                        let (sk, sl) = self.pair_shells(q);
                        let off = store.offs[bi] as usize;
                        let nels = self.dims[si] * self.dims[sj] * self.dims[sk] * self.dims[sl];
                        self.contract(&store.data[off..off + nels], si, sj, sk, sl, dl, &mut j, &mut k);
                        (j, k)
                    },
                )
                .reduce(
                    || (vec![0.0f64; n * n], vec![0.0f64; n * n]),
                    |(mut ja, mut ka), (jb, kb)| {
                        for (x, y) in ja.iter_mut().zip(jb.iter()) {
                            *x += y;
                        }
                        for (x, y) in ka.iter_mut().zip(kb.iter()) {
                            *x += y;
                        }
                        (ja, ka)
                    },
                );
            jl.copy_from_slice(&aj);
            kl.copy_from_slice(&ak);
        } else {
            self.accumulate_direct(dl, jl, kl);
        }
    }

    /// 直接 SCF：每轮重新计算通过 Schwarz + 密度筛选的块
    fn accumulate_direct(&self, dl: &[f64], jl: &mut [f64], kl: &mut [f64]) {
        let n = self.nato;
        let npair = self.pairs.len();
        let dmax = dl.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let scr = (self.screen * 1.0e-3).max(1.0e-14);
        let parts: Vec<(Vec<f64>, Vec<f64>)> = (0..npair)
            .into_par_iter()
            .with_min_len(8)
            .fold(
                || (vec![0.0f64; n * n], vec![0.0f64; n * n], vec![0.0f64; self.max_block]),
                |(mut j, mut k, mut buf), p| {
                    let (i, jsh) = self.pair_shells(p);
                    let qp = self.pairs[p].q;
                    for q in 0..=p {
                        if qp * self.pairs[q].q * dmax < scr {
                            continue;
                        }
                        let (k_, l_) = self.pair_shells(q);
                        let nels = self.dims[i] * self.dims[jsh] * self.dims[k_] * self.dims[l_];
                        let mut shls = [i as c_int, jsh as c_int, k_ as c_int, l_ as c_int];
                        self.call_cint(&mut buf[..nels], &mut shls);
                        self.contract(&buf[..nels], i, jsh, k_, l_, dl, &mut j, &mut k);
                    }
                    (j, k, buf)
                },
            )
            .map(|(j, k, _)| (j, k))
            .collect();
        let mut ja = vec![0.0f64; n * n];
        let mut ka = vec![0.0f64; n * n];
        for (j, k) in parts {
            for (x, y) in ja.iter_mut().zip(j.iter()) {
                *x += y;
            }
            for (x, y) in ka.iter_mut().zip(k.iter()) {
                *x += y;
            }
        }
        jl.copy_from_slice(&ja);
        kl.copy_from_slice(&ka);
    }

    /// 开壳层：J 用总密度，K 用 α/β 密度
    fn accumulate_spin(&self, dt: &[f64], da: &[f64], db: &[f64], jl: &mut [f64], kal: &mut [f64], kbl: &mut [f64]) {
        let n = self.nato;
        let zero = || (vec![0.0f64; n * n], vec![0.0f64; n * n], vec![0.0f64; n * n]);
        let (pj, pka, pkb) = match &self.store {
            Some(store) => (0..store.len())
                .into_par_iter()
                .with_min_len(1024)
                .fold(zero, |(mut j, mut ka, mut kb), bi| {
                    let p = store.ps[bi] as usize;
                    let q = store.qs[bi] as usize;
                    let (si, sj) = self.pair_shells(p);
                    let (sk, sl) = self.pair_shells(q);
                    let off = store.offs[bi] as usize;
                    let nels = self.dims[si] * self.dims[sj] * self.dims[sk] * self.dims[sl];
                    self.contract_spin(&store.data[off..off + nels], si, sj, sk, sl, dt, da, db, &mut j, &mut ka, &mut kb);
                    (j, ka, kb)
                })
                .reduce(zero, |(mut ja, mut aa, mut ba), (jb, ab, bb)| {
                    for (x, y) in ja.iter_mut().zip(jb.iter()) {
                        *x += y;
                    }
                    for (x, y) in aa.iter_mut().zip(ab.iter()) {
                        *x += y;
                    }
                    for (x, y) in ba.iter_mut().zip(bb.iter()) {
                        *x += y;
                    }
                    (ja, aa, ba)
                }),
            None => {
                let dmax = dt.iter().fold(0.0f64, |m, v| m.max(v.abs()));
                let scr = (self.screen * 1.0e-3).max(1.0e-14);
                let npair = self.pairs.len();
                let parts: Vec<(Vec<f64>, Vec<f64>, Vec<f64>)> = (0..npair)
                    .into_par_iter()
                    .with_min_len(8)
                    .fold(
                        || (vec![0.0f64; n * n], vec![0.0f64; n * n], vec![0.0f64; n * n], vec![0.0f64; self.max_block]),
                        |(mut j, mut ka, mut kb, mut buf), p| {
                            let (i, jsh) = self.pair_shells(p);
                            let qp = self.pairs[p].q;
                            for q in 0..=p {
                                if qp * self.pairs[q].q * dmax < scr {
                                    continue;
                                }
                                let (k_, l_) = self.pair_shells(q);
                                let nels = self.dims[i] * self.dims[jsh] * self.dims[k_] * self.dims[l_];
                                let mut shls = [i as c_int, jsh as c_int, k_ as c_int, l_ as c_int];
                                self.call_cint(&mut buf[..nels], &mut shls);
                                self.contract_spin(&buf[..nels], i, jsh, k_, l_, dt, da, db, &mut j, &mut ka, &mut kb);
                            }
                            (j, ka, kb, buf)
                        },
                    )
                    .map(|(j, ka, kb, _)| (j, ka, kb))
                    .collect();
                let mut ja = vec![0.0f64; n * n];
                let mut aa = vec![0.0f64; n * n];
                let mut ba = vec![0.0f64; n * n];
                for (j, ka, kb) in parts {
                    for (x, y) in ja.iter_mut().zip(j.iter()) {
                        *x += y;
                    }
                    for (x, y) in aa.iter_mut().zip(ka.iter()) {
                        *x += y;
                    }
                    for (x, y) in ba.iter_mut().zip(kb.iter()) {
                        *x += y;
                    }
                }
                (ja, aa, ba)
            }
        };
        jl.copy_from_slice(&pj);
        kal.copy_from_slice(&pka);
        kbl.copy_from_slice(&pkb);
    }

    /// 把一个壳层块累加进 J/K
    ///
    /// libcint 的四维积分按**列主序**存放：`idx = a + di*(b + dj*(c + dk*d))`。
    ///
    /// 块的布局始终是规范顺序 (i ≥ j, k ≥ l, ij ≥ kl)。四种壳层互不相同时
    /// 走快速路径（块内 6 次点积）；存在简并（i=j 或 k=l 或同一个壳层对）时
    /// 走通用路径，把 8 个置换逐个枚举并按壳层四元组去重，保证与朴素的四重
    /// 循环求和逐项一致。
    #[inline]
    fn contract(&self, buf: &[f64], si: usize, sj: usize, sk: usize, sl: usize, dl: &[f64], jl: &mut [f64], kl: &mut [f64]) {
        let n = self.nato;
        let (di, dj, dk, dl_) = (self.dims[si], self.dims[sj], self.dims[sk], self.dims[sl]);
        let (oi, oj, ok, ol) = (self.offs[si], self.offs[sj], self.offs[sk], self.offs[sl]);
        // 列主序步长
        let s0 = 1usize;
        let s1 = di;
        let s2 = di * dj;
        let s3 = di * dj * dk;
        let pij = si * (si + 1) / 2 + sj;
        let pkl = sk * (sk + 1) / 2 + sl;
        if si != sj && sk != sl && pij != pkl {
            // ---- 快速路径 ----
            let bkl = dk * dl_;
            assert!(bkl <= 256, "壳层块过大: {bkl}");
            let mut t = [0.0f64; 256];
            t[..bkl].fill(0.0);
            // J[i,j] 与 J[k,l] 共用一次块遍历
            for a in 0..di {
                for b in 0..dj {
                    let base = a * s0 + b * s1;
                    let dav = dl[(oi + a) * n + oj + b];
                    let mut sum = 0.0;
                    for c in 0..dk {
                        let dkrow = (ok + c) * n + ol;
                        let brow = base + c * s2;
                        for d in 0..dl_ {
                            let w = buf[brow + d * s3];
                            sum += w * dl[dkrow + d];
                            t[c * dl_ + d] += w * dav;
                        }
                    }
                    let v = 2.0 * sum;
                    jl[(oi + a) * n + oj + b] += v;
                    jl[(oj + b) * n + oi + a] += v;
                }
            }
            for c in 0..dk {
                for d in 0..dl_ {
                    let v = 2.0 * t[c * dl_ + d];
                    jl[(ok + c) * n + ol + d] += v;
                    jl[(ol + d) * n + ok + c] += v;
                }
            }
            // K[i,k] = Σ (ij|kl) D[j,l]
            for a in 0..di {
                for c in 0..dk {
                    let mut sum = 0.0;
                    for b in 0..dj {
                        let drow = (oj + b) * n + ol;
                        let brow = a * s0 + b * s1 + c * s2;
                        for d in 0..dl_ {
                            sum += buf[brow + d * s3] * dl[drow + d];
                        }
                    }
                    kl[(oi + a) * n + ok + c] += sum;
                    kl[(ok + c) * n + oi + a] += sum;
                }
            }
            // K[i,l] = Σ (ij|kl) D[j,k]
            for a in 0..di {
                for d in 0..dl_ {
                    let mut sum = 0.0;
                    for b in 0..dj {
                        let brow = a * s0 + b * s1;
                        let drow = (oj + b) * n + ok;
                        for c in 0..dk {
                            sum += buf[brow + c * s2 + d * s3] * dl[drow + c];
                        }
                    }
                    kl[(oi + a) * n + ol + d] += sum;
                    kl[(ol + d) * n + oi + a] += sum;
                }
            }
            // K[j,k] = Σ (ij|kl) D[i,l]
            for b in 0..dj {
                for c in 0..dk {
                    let mut sum = 0.0;
                    for a in 0..di {
                        let brow = a * s0 + b * s1 + c * s2;
                        let drow = (oi + a) * n + ol;
                        for d in 0..dl_ {
                            sum += buf[brow + d * s3] * dl[drow + d];
                        }
                    }
                    kl[(oj + b) * n + ok + c] += sum;
                    kl[(ok + c) * n + oj + b] += sum;
                }
            }
            // K[j,l] = Σ (ij|kl) D[i,k]
            for b in 0..dj {
                for d in 0..dl_ {
                    let mut sum = 0.0;
                    for a in 0..di {
                        let brow = a * s0 + b * s1;
                        let drow = (oi + a) * n + ok;
                        for c in 0..dk {
                            sum += buf[brow + c * s2 + d * s3] * dl[drow + c];
                        }
                    }
                    kl[(oj + b) * n + ol + d] += sum;
                    kl[(ol + d) * n + oj + b] += sum;
                }
            }
            return;
        }
        // ---- 通用路径：枚举 8 个置换并按壳层四元组去重 ----
        let shs = [si, sj, sk, sl];
        let dims = [di, dj, dk, dl_];
        let offs = [oi, oj, ok, ol];
        let strides = [s0, s1, s2, s3];
        let mut seen: [[usize; 4]; 8] = [[usize::MAX; 4]; 8];
        let mut nseen = 0usize;
        for perm in perms().iter() {
            let key = [shs[perm[0]], shs[perm[1]], shs[perm[2]], shs[perm[3]]];
            if seen[..nseen].contains(&key) {
                continue;
            }
            seen[nseen] = key;
            nseen += 1;
            let dm = [dims[perm[0]], dims[perm[1]], dims[perm[2]], dims[perm[3]]];
            let of = [offs[perm[0]], offs[perm[1]], offs[perm[2]], offs[perm[3]]];
            let st = [strides[perm[0]], strides[perm[1]], strides[perm[2]], strides[perm[3]]];
            for a in 0..dm[0] {
                let ia = of[0] + a;
                for b in 0..dm[1] {
                    let ib = of[1] + b;
                    for c in 0..dm[2] {
                        let ic = of[2] + c;
                        let base = a * st[0] + b * st[1] + c * st[2];
                        for d in 0..dm[3] {
                            let w = buf[base + d * st[3]];
                            if w == 0.0 {
                                continue;
                            }
                            let id = of[3] + d;
                            jl[ia * n + ib] += w * dl[ic * n + id];
                            kl[ia * n + ic] += w * dl[ib * n + id];
                        }
                    }
                }
            }
        }
    }

    /// 开壳层版本：J 用总密度，K_α/K_β 分别用 α/β 密度
    #[inline]
    fn contract_spin(&self, buf: &[f64], si: usize, sj: usize, sk: usize, sl: usize, dt: &[f64], da: &[f64], db: &[f64], jl: &mut [f64], kal: &mut [f64], kbl: &mut [f64]) {
        let n = self.nato;
        let (di, dj, dk, dl_) = (self.dims[si], self.dims[sj], self.dims[sk], self.dims[sl]);
        let (oi, oj, ok, ol) = (self.offs[si], self.offs[sj], self.offs[sk], self.offs[sl]);
        let shs = [si, sj, sk, sl];
        let dims = [di, dj, dk, dl_];
        let offs = [oi, oj, ok, ol];
        let strides = [1usize, di, di * dj, di * dj * dk];
        let mut seen: [[usize; 4]; 8] = [[usize::MAX; 4]; 8];
        let mut nseen = 0usize;
        for perm in perms().iter() {
            let key = [shs[perm[0]], shs[perm[1]], shs[perm[2]], shs[perm[3]]];
            if seen[..nseen].contains(&key) {
                continue;
            }
            seen[nseen] = key;
            nseen += 1;
            let dm = [dims[perm[0]], dims[perm[1]], dims[perm[2]], dims[perm[3]]];
            let of = [offs[perm[0]], offs[perm[1]], offs[perm[2]], offs[perm[3]]];
            let st = [strides[perm[0]], strides[perm[1]], strides[perm[2]], strides[perm[3]]];
            for a in 0..dm[0] {
                let ia = of[0] + a;
                for b in 0..dm[1] {
                    let ib = of[1] + b;
                    for c in 0..dm[2] {
                        let ic = of[2] + c;
                        let base = a * st[0] + b * st[1] + c * st[2];
                        for d in 0..dm[3] {
                            let w = buf[base + d * st[3]];
                            if w == 0.0 {
                                continue;
                            }
                            let id = of[3] + d;
                            jl[ia * n + ib] += w * dt[ic * n + id];
                            kal[ia * n + ic] += w * da[ib * n + id];
                            kbl[ia * n + ic] += w * db[ib * n + id];
                        }
                    }
                }
            }
        }
    }
}

/// (ij|kl) 的 8 个置换（相对规范顺序的壳层槽位），已在壳层层面去重
///
/// 规范块 (i ≥ j, k ≥ l, ij ≥ kl) 的 8 个置换里，有 4 个是「交换 bra/ket」，
/// 另 4 个是壳层内交换。简并时（i=j 等）重复的置换被去重，这样每个不同的
/// 有序四元组恰好累加一次，与朴素的四重循环求和完全等价。
fn perms() -> [[usize; 4]; 8] {
    [[0, 1, 2, 3], [1, 0, 2, 3], [0, 1, 3, 2], [1, 0, 3, 2], [2, 3, 0, 1], [3, 2, 0, 1], [2, 3, 1, 0], [3, 2, 1, 0]]
}

/// 壳层内 libcint 顺序到高斯顺序的置换（与 `rsint.rs` 中的实现一致）
fn shell_perm(len: usize) -> Vec<usize> {
    match len {
        1 => vec![0],
        3 => vec![0, 1, 2],
        5 => vec![2, 3, 1, 4, 0],
        7 => vec![3, 4, 2, 5, 1, 6, 0],
        9 => vec![4, 5, 3, 6, 2, 7, 1, 8, 0],
        11 => vec![5, 6, 4, 7, 3, 8, 2, 9, 1, 10, 0],
        6 => vec![0, 3, 5, 1, 2, 4],
        10 => vec![0, 6, 9, 3, 1, 2, 5, 8, 7, 4],
        15 => vec![14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
        21 => vec![20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
        _ => panic!("未知壳层长度: {len}"),
    }
}
