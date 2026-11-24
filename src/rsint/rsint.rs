// 基于 rsint crate 的分子积分计算器
//
// `rsint` 是 libcint 的机械式 C→Rust 直译版本, 导出与 libcint 相同的 C ABI。
// 本模块直接以 Rust 路径调用这些函数, 因此不需要链接外部 `libcint.a`。
// 0.1.1 起内置纯 Rust 的 C 运行时 shim, build.rs 也补上了 Windows 分支,
// 可以直接在 Windows 上构建。
//
// 注意: `rsint` 目前只能在 `nightly-2025-02-01` 工具链上编译, 因此该后端由
// `use-rsint` feature 控制且默认关闭。
#![allow(non_snake_case)]

use super::Mat2e;
use crate::base::Mole;
use ndarray::{Array, Array2, Array4, ShapeBuilder, s};
use std::f64::consts::PI;
use std::os::raw::{c_double, c_int};

use ::rsint::src::autocode::intor1::{cint1e_kin_cart, cint1e_kin_sph};
use ::rsint::src::cint_bas::{CINTcgto_cart, CINTcgto_spheric};
use ::rsint::src::cint1e::{cint1e_nuc_cart, cint1e_nuc_sph, cint1e_ovlp_cart, cint1e_ovlp_sph};
use ::rsint::src::cint2e::{CINTOpt, cint2e_cart, cint2e_cart_optimizer, cint2e_sph, cint2e_sph_optimizer};
use ::rsint::src::optimizer::CINTdel_optimizer;

/// libcint `_gaussian_int(n, a) = ∫₀^∞ x^n e^{-a x²} dx` 的等价实现
fn gaussian_int(n: c_int, alpha: c_double) -> f64 {
    match n {
        0 => (PI / (4.0 * alpha)).sqrt(),
        1 => 0.5 / alpha,
        _ => (n as f64 - 1.0) / (2.0 * alpha) * gaussian_int(n - 2, alpha),
    }
}

/// libcint `CINTgto_norm` 的等价实现 (rsint 未转译 misc.c)
fn gto_norm(ang: c_int, alp: c_double) -> f64 {
    1.0 / gaussian_int(ang * 2 + 2, 2.0 * alp).sqrt()
}

/// 归一化基函数
///
/// 对每个基函数应用 `CINTgto_norm` 归一化因子。
fn norm_basis(bas: &[i32], env: &[f64]) -> Vec<f64> {
    let mut env_new = env.to_vec();
    let nbas = bas.len() / 8;
    for i in 0..nbas {
        let idx = i * 8;
        let ang = bas[idx + 1];
        let ialp = bas[idx + 5] as usize;
        let icoe = bas[idx + 6] as usize;
        let tot = icoe - ialp;
        for j in 0..tot {
            let alp = env[ialp + j];
            let coe = env[icoe + j];
            env_new[icoe + j] = coe * gto_norm(ang, alp);
        }
    }
    env_new
}

/// rsint 后端的分子积分计算器
///
/// 字段含义与 libcint 的 `atm` / `bas` / `env` 数组一致。
pub struct Rsint<'m> {
    /// 分子结构引用
    pub mole: &'m Mole,
    /// rsint 的 atm 数组
    pub atm: Vec<i32>,
    /// rsint 的 bas 数组
    pub bas: Vec<i32>,
    /// rsint 的 env 数组
    pub env: Vec<f64>,
    /// 原子轨道数量 [球谐, 笛卡尔]
    pub nato: [usize; 2],
    /// 校正系数 [球谐, 笛卡尔]
    pub fixs: [Vec<f64>; 2],
}

impl<'m> Rsint<'m> {
    /// 创建 Rsint 实例
    pub fn new(mole: &'m Mole) -> Self {
        let natm = mole.atoms.len();
        let basis = mole.basis();
        let shls = &basis.shls.shls;
        let nato = [basis.shls.nato("sph"), basis.shls.nato("car")];
        let nshl = basis.shls.len();
        let mut atm = vec![0; natm * 6];
        let mut bas = vec![0; nshl * 8];
        let mut env = vec![0.0; 20];
        let mut idx = 20;
        for i in 0..natm {
            let atom = mole.atoms.get(i);
            let nuc = atom.nuc as c_int;
            let [x, y, z] = atom.xyz;
            atm[i * 6 + 0] = nuc;
            atm[i * 6 + 1] = idx;
            env.push(x);
            env.push(y);
            env.push(z);
            idx += 3;
        }

        for i in 0..nshl {
            bas[i * 8 + 0] = shls[i].atm as c_int;
            bas[i * 8 + 1] = shls[i].ang as c_int;
            let alps = &shls[i].alps;
            let coes = &shls[i].coes;
            let num = alps.len() as c_int;
            bas[i * 8 + 2] = num; // 原始gto的数量
            bas[i * 8 + 3] = 1; // 收缩gto的数量
            bas[i * 8 + 4] = 0;
            env.extend(alps.clone());
            bas[i * 8 + 5] = idx; // 该层基函数在env中的起始位置
            idx += num;
            env.extend(coes.clone());
            bas[i * 8 + 6] = idx;
            idx += num;
        }
        let env = norm_basis(&bas, &env); // 归一化

        // 计算校正系数
        let fixs = {
            let smat_sph = &calc_int1e_sph(nato[0], &atm, &bas, &env, [1, 0, 0])[0];
            let smat_car = &calc_int1e_car(nato[1], &atm, &bas, &env, [1, 0, 0])[0];
            let vals_sph = (0..smat_sph.nrows())
                .map(|i| {
                    let val = smat_sph[(i, i)];
                    if val == 0.0 { 0.0 } else { 1.0 / val.sqrt() }
                })
                .collect();
            let vals_car = (0..smat_car.nrows())
                .map(|i| {
                    let val = smat_car[(i, i)];
                    if val == 0.0 { 0.0 } else { 1.0 / val.sqrt() }
                })
                .collect();
            [vals_sph, vals_car]
        };

        Self {
            mole,
            atm,
            bas,
            env,
            nato,
            fixs,
        }
    }

    /// 获取校正系数
    fn get_fixs(&self, form: &str) -> &Vec<f64> {
        match form {
            "sph" => &self.fixs[0],
            "car" => &self.fixs[1],
            _ => panic!("Invalid form"),
        }
    }

    /// 对单电子积分矩阵应用校正系数
    fn fix_int1e(&self, int1e: &mut Array2<f64>, form: &str) {
        let vals = self.get_fixs(form);
        let nato = vals.len();
        for i in 0..nato {
            for j in 0..nato {
                int1e[(i, j)] *= vals[i] * vals[j];
            }
        }
    }

    /// 对双电子积分张量应用校正系数
    fn fix_int2e(&self, int2e: &mut Array4<f64>, form: &str) {
        let vals = self.get_fixs(form);
        let nato = vals.len();
        for i in 0..nato {
            for j in 0..nato {
                for k in 0..nato {
                    for l in 0..nato {
                        int2e[(i, j, k, l)] *= vals[i] * vals[j] * vals[k] * vals[l];
                    }
                }
            }
        }
    }

    /// 计算单电子积分
    pub fn calc_int1e(&self, form: &str, cals: [usize; 3]) -> [Array2<f64>; 3] {
        let [mut smat, mut tmat, mut vmat] = match form {
            "sph" => calc_int1e_sph(self.nato[0], &self.atm, &self.bas, &self.env, cals),
            "car" => calc_int1e_car(self.nato[1], &self.atm, &self.bas, &self.env, cals),
            _ => panic!("Invalid form"),
        };
        // 系数矫正&顺序矫正
        if cals[0] == 1 {
            self.fix_int1e(&mut smat, form);
            smat = self.ord_int1e(&smat, form);
        }
        if cals[1] == 1 {
            self.fix_int1e(&mut tmat, form);
            tmat = self.ord_int1e(&tmat, form);
        }
        if cals[2] == 1 {
            self.fix_int1e(&mut vmat, form);
            vmat = self.ord_int1e(&vmat, form);
        }
        [smat, tmat, vmat]
    }

    /// 计算双电子积分
    pub fn calc_int2e(&self, form: &str) -> Array4<f64> {
        let mut mat2e = match form {
            "sph" => calc_int2e_sph(self.nato[0], &self.atm, &self.bas, &self.env),
            "car" => calc_int2e_car(self.nato[1], &self.atm, &self.bas, &self.env),
            _ => panic!("Invalid form"),
        };
        self.fix_int2e(&mut mat2e, form);
        let mat2e = self.ord_int2e(mat2e, form);
        mat2e
    }

    /// 计算双电子积分并直接写入压缩存储（利用 (ij|kl) 的 8 重置换对称）
    ///
    /// 与 [`Rsint::calc_int2e`] 相比：
    /// - 只遍历规范化代表元（`i >= j`、`k >= l`、`(ij) >= (kl)`），
    ///   积分计算量降为 1/8；
    /// - 不再分配 `nato⁴` 的稠密中间数组（例如 102 个基函数时就是 864 MB），
    ///   内存占用和访存量都大幅下降；
    /// - 归一化因子按 libcint 的原始顺序施加，壳层内分量直接写到高斯顺序。
    pub fn calc_int2e_packed(&self, form: &str) -> Mat2e {
        let (nmat, fix, nums) = match form {
            "sph" => (self.nato[0], &self.fixs[0], [1usize, 3, 5, 7, 9, 11, 13]),
            "car" => (self.nato[1], &self.fixs[1], [1usize, 3, 6, 10, 15, 21, 28]),
            _ => panic!("Invalid form"),
        };
        let natm = self.atm.len() / 6;
        let nbas = self.bas.len() / 8;

        // 每个壳层的起始轨道位置、分量数，以及 libcint 顺序 -> 高斯顺序的置换
        let mut locs = vec![0usize; nbas];
        let mut dims = vec![0usize; nbas];
        let mut perms: Vec<Vec<usize>> = Vec::with_capacity(nbas);
        let mut loc = 0usize;
        for i in 0..nbas {
            let ang = self.bas[i * 8 + 1] as usize;
            locs[i] = loc;
            dims[i] = nums[ang];
            loc += nums[ang];
            // len2idxs 给出 高斯位置 -> libcint 位置，这里取逆
            let idx = len2idxs(dims[i]);
            let mut perm = vec![0usize; dims[i]];
            for (gauss, libcint) in idx.iter().enumerate() {
                perm[*libcint] = gauss;
            }
            perms.push(perm);
        }
        assert_eq!(loc, nmat, "基函数数量与 nato 不一致");

        let mut emat = Mat2e::new(nmat);
        let mut opt: *mut CINTOpt = std::ptr::null_mut();
        unsafe {
            match form {
                "sph" => cint2e_sph_optimizer(
                    &mut opt,
                    self.atm.as_ptr() as *mut c_int,
                    natm as c_int,
                    self.bas.as_ptr() as *mut c_int,
                    nbas as c_int,
                    self.env.as_ptr() as *mut c_double,
                ),
                _ => cint2e_cart_optimizer(
                    &mut opt,
                    self.atm.as_ptr() as *mut c_int,
                    natm as c_int,
                    self.bas.as_ptr() as *mut c_int,
                    nbas as c_int,
                    self.env.as_ptr() as *mut c_double,
                ),
            }

            for i in 0..nbas {
                for j in 0..=i {
                    let pij = i * (i + 1) / 2 + j;
                    for k in 0..nbas {
                        for l in 0..=k {
                            let pkl = k * (k + 1) / 2 + l;
                            // 8 重对称：(ij|kl) = (ji|kl) = (ij|lk) = (kl|ij) = ...
                            // 只算满足 i>=j、k>=l、(ij)>=(kl) 的代表元
                            if pij < pkl {
                                continue;
                            }
                            let (ui, uj, uk, ul) = (locs[i], locs[j], locs[k], locs[l]);
                            let (di, dj, dk, dl) = (dims[i], dims[j], dims[k], dims[l]);
                            let mut buf = vec![0.0f64; di * dj * dk * dl];
                            let mut shls = [i as c_int, j as c_int, k as c_int, l as c_int];
                            match form {
                                "sph" => cint2e_sph(
                                    buf.as_mut_ptr(),
                                    shls.as_mut_ptr(),
                                    self.atm.as_ptr() as *mut c_int,
                                    natm as c_int,
                                    self.bas.as_ptr() as *mut c_int,
                                    nbas as c_int,
                                    self.env.as_ptr() as *mut c_double,
                                    opt,
                                ),
                                _ => cint2e_cart(
                                    buf.as_mut_ptr(),
                                    shls.as_mut_ptr(),
                                    self.atm.as_ptr() as *mut c_int,
                                    natm as c_int,
                                    self.bas.as_ptr() as *mut c_int,
                                    nbas as c_int,
                                    self.env.as_ptr() as *mut c_double,
                                    opt,
                                ),
                            };
                            let (pi, pj, pk, pl) = (&perms[i], &perms[j], &perms[k], &perms[l]);
                            for a in 0..di {
                                let fa = fix[ui + a];
                                let ga = ui + pi[a];
                                for b in 0..dj {
                                    let fab = fa * fix[uj + b];
                                    let gb = uj + pj[b];
                                    for c in 0..dk {
                                        let fabc = fab * fix[uk + c];
                                        let gc = uk + pk[c];
                                        let base = a + di * (b + dj * (c + dk * 0));
                                        for d in 0..dl {
                                            let val = buf[base + di * dj * dk * d] * fabc * fix[ul + d];
                                            emat.set(ga, gb, gc, ul + pl[d], val);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            // cint2e 与 optimizer 各自定义了布局相同的 CINTOpt, 这里做一次裸指针转换
            CINTdel_optimizer(&mut opt as *mut *mut CINTOpt as *mut *mut ::rsint::src::optimizer::CINTOpt);
        }
        emat
    }

    /// 将单电子积分的原子轨道顺序转为和高斯一致
    fn ord_int1e(&self, mat1e_raw: &Array2<f64>, form: &str) -> Array2<f64> {
        let mut mat1e_fix = mat1e_raw.clone();
        let uplo = self.mole.basis().shls.uplo(form);
        let nshl = uplo.len();
        for si in 0..nshl {
            let [ui, li] = uplo[si];
            let leni = li - ui;
            let idxi = len2idxs(leni);
            for sj in 0..nshl {
                let [uj, lj] = uplo[sj];
                let lenj = lj - uj;
                let idxj = len2idxs(lenj);
                if leni <= 3 && lenj <= 3 {
                    continue;
                }
                for i in 0..leni {
                    for j in 0..lenj {
                        mat1e_fix[(ui + i, uj + j)] = mat1e_raw[(ui + idxi[i], uj + idxj[j])]
                    }
                }
            }
        }
        mat1e_fix
    }

    /// 将双电子积分的原子轨道顺序转为和高斯一致
    fn ord_int2e(&self, mat2e_raw: Array4<f64>, form: &str) -> Array4<f64> {
        let mut mat2e_fix = mat2e_raw.clone();
        let uplo = self.mole.basis().shls.uplo(form);
        let nshl = uplo.len();
        for si in 0..nshl {
            let [ui, li] = uplo[si];
            let leni = li - ui;
            let idxi = len2idxs(leni);
            for sj in 0..nshl {
                let [uj, lj] = uplo[sj];
                let lenj = lj - uj;
                let idxj = len2idxs(lenj);
                for sk in 0..nshl {
                    let [uk, lk] = uplo[sk];
                    let lenk = lk - uk;
                    let idxk = len2idxs(lenk);
                    for sl in 0..nshl {
                        let [ul, ll] = uplo[sl];
                        let lenl = ll - ul;
                        let idxl = len2idxs(lenl);
                        for i in 0..leni {
                            for j in 0..lenj {
                                for k in 0..lenk {
                                    for l in 0..lenl {
                                        mat2e_fix[(ui + i, uj + j, uk + k, ul + l)] =
                                            mat2e_raw[(ui + idxi[i], uj + idxj[j], uk + idxk[k], ul + idxl[l])];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        mat2e_fix
    }
}

/// 计算球谐单电子积分
fn calc_int1e_sph(nato: usize, atm: &[i32], bas: &[i32], env: &[f64], cals: [usize; 3]) -> [Array2<f64>; 3] {
    let nbas = bas.len() / 8;
    let natm = atm.len() / 6;
    let mut smat = Array2::<f64>::zeros((nato, nato).f()); // 重叠矩阵
    let mut tmat = Array2::<f64>::zeros((nato, nato).f()); // 动能矩阵
    let mut vmat = Array2::<f64>::zeros((nato, nato).f()); // 势能矩阵

    let mut locs = vec![0 as usize; nbas]; // 每个基函数在矩阵中的起始位置
    let nums: [usize; 7] = [1, 3, 5, 7, 9, 11, 13]; // 每个角动量对应的数量
    let mut loc = 0usize;
    for i in 0..nbas {
        locs[i] = loc;
        let ang = bas[i * 8 + 1] as usize;
        loc += nums[ang];
    }
    unsafe {
        let mut dims = vec![0 as usize; nbas]; // 每个基函数的原子轨道数量
        for i in 0..nbas {
            dims[i] = CINTcgto_spheric(i as c_int, bas.as_ptr()) as usize;
        }
        for i in 0..nbas {
            let loci = locs[i];
            let dimi = dims[i];
            assert!(loci + dimi <= nato, "数组越界,{i},{loci},{dimi},{nato}");
            for j in 0..nbas {
                let locj = locs[j];
                let dimj = dims[j];
                assert!(locj + dimj <= nato, "数组越界,{j},{locj},{dimj},{nato}");
                let mut buf = vec![0.0f64; dimi * dimj];
                let mut shls = [i as c_int, j as c_int];
                let shape = (dimi, dimj);
                if cals[0] == 1 {
                    // 计算重叠积分
                    cint1e_ovlp_sph(
                        buf.as_mut_ptr(),
                        shls.as_mut_ptr(),
                        atm.as_ptr() as *mut c_int,
                        natm as c_int,
                        bas.as_ptr() as *mut c_int,
                        nbas as c_int,
                        env.as_ptr() as *mut c_double,
                        std::ptr::null_mut(),
                    );
                    let buf_array = Array::from_shape_vec(shape.f(), buf.clone()).unwrap();
                    smat.slice_mut(s![loci..loci + dimi, locj..locj + dimj])
                        .assign(&buf_array);
                }
                if cals[1] == 1 {
                    // 计算动能积分
                    cint1e_kin_sph(
                        buf.as_mut_ptr(),
                        shls.as_mut_ptr(),
                        atm.as_ptr() as *mut c_int,
                        natm as c_int,
                        bas.as_ptr() as *mut c_int,
                        nbas as c_int,
                        env.as_ptr() as *mut c_double,
                    );
                    let buf_array = Array::from_shape_vec(shape.f(), buf.clone()).unwrap();
                    tmat.slice_mut(s![loci..loci + dimi, locj..locj + dimj])
                        .assign(&buf_array);
                }
                if cals[2] == 1 {
                    // 计算核引积分
                    cint1e_nuc_sph(
                        buf.as_mut_ptr(),
                        shls.as_mut_ptr(),
                        atm.as_ptr() as *mut c_int,
                        natm as c_int,
                        bas.as_ptr() as *mut c_int,
                        nbas as c_int,
                        env.as_ptr() as *mut c_double,
                        std::ptr::null_mut(),
                    );
                    let buf_array = Array::from_shape_vec(shape.f(), buf.clone()).unwrap();
                    vmat.slice_mut(s![loci..loci + dimi, locj..locj + dimj])
                        .assign(&buf_array);
                }
            }
        }
    }
    [smat, tmat, vmat]
}

/// 计算笛卡尔单电子积分
fn calc_int1e_car(nato: usize, atm: &[i32], bas: &[i32], env: &[f64], cals: [usize; 3]) -> [Array2<f64>; 3] {
    let nbas = bas.len() / 8;
    let natm = atm.len() / 6;
    let mut smat = Array2::<f64>::zeros((nato, nato).f()); // 重叠矩阵
    let mut tmat = Array2::<f64>::zeros((nato, nato).f()); // 动能矩阵
    let mut vmat = Array2::<f64>::zeros((nato, nato).f()); // 势能矩阵

    let mut locs = vec![0 as usize; nbas]; // 每个基函数在矩阵中的起始位置
    let nums: [usize; 7] = [1, 3, 6, 10, 15, 21, 28]; // 每个角动量对应的数量
    let mut loc = 0usize;
    for i in 0..nbas {
        locs[i] = loc;
        let ang = bas[i * 8 + 1] as usize;
        loc += nums[ang];
    }
    unsafe {
        let mut dims = vec![0 as usize; nbas]; // 每个基函数的原子轨道数量
        for i in 0..nbas {
            dims[i] = CINTcgto_cart(i as c_int, bas.as_ptr()) as usize;
        }
        for i in 0..nbas {
            let loci = locs[i];
            let dimi = dims[i];
            assert!(loci + dimi <= nato, "数组越界,{i},{loci},{dimi},{nato}");
            for j in 0..nbas {
                let locj = locs[j];
                let dimj = dims[j];
                assert!(locj + dimj <= nato, "数组越界,{j},{locj},{dimj},{nato}");
                let mut buf = vec![0.0f64; dimi * dimj];
                let mut shls = [i as c_int, j as c_int];
                let shape = (dimi, dimj);
                if cals[0] == 1 {
                    // 计算重叠积分
                    cint1e_ovlp_cart(
                        buf.as_mut_ptr(),
                        shls.as_mut_ptr(),
                        atm.as_ptr() as *mut c_int,
                        natm as c_int,
                        bas.as_ptr() as *mut c_int,
                        nbas as c_int,
                        env.as_ptr() as *mut c_double,
                        std::ptr::null_mut(),
                    );
                    let buf_array = Array::from_shape_vec(shape.f(), buf.clone()).unwrap();
                    smat.slice_mut(s![loci..loci + dimi, locj..locj + dimj])
                        .assign(&buf_array);
                }
                if cals[1] == 1 {
                    // 计算动能积分
                    cint1e_kin_cart(
                        buf.as_mut_ptr(),
                        shls.as_mut_ptr(),
                        atm.as_ptr() as *mut c_int,
                        natm as c_int,
                        bas.as_ptr() as *mut c_int,
                        nbas as c_int,
                        env.as_ptr() as *mut c_double,
                    );
                    let buf_array = Array::from_shape_vec(shape.f(), buf.clone()).unwrap();
                    tmat.slice_mut(s![loci..loci + dimi, locj..locj + dimj])
                        .assign(&buf_array);
                }
                if cals[2] == 1 {
                    // 计算核引积分
                    cint1e_nuc_cart(
                        buf.as_mut_ptr(),
                        shls.as_mut_ptr(),
                        atm.as_ptr() as *mut c_int,
                        natm as c_int,
                        bas.as_ptr() as *mut c_int,
                        nbas as c_int,
                        env.as_ptr() as *mut c_double,
                        std::ptr::null_mut(),
                    );
                    let buf_array = Array::from_shape_vec(shape.f(), buf.clone()).unwrap();
                    vmat.slice_mut(s![loci..loci + dimi, locj..locj + dimj])
                        .assign(&buf_array);
                }
            }
        }
    }
    [smat, tmat, vmat]
}

/// 计算球谐双电子积分
fn calc_int2e_sph(nmat: usize, atm: &[i32], bas: &[i32], env: &[f64]) -> Array4<f64> {
    let natm = atm.len() / 6; // 原子数量
    let nbas = bas.len() / 8; // 基函数数量

    let mut locs = vec![0 as usize; nbas]; // 每个基函数在矩阵中的起始位置
    let nums: [usize; 7] = [1, 3, 5, 7, 9, 11, 13]; // 每个角动量对应的数量
    let mut loc = 0usize;
    for i in 0..nbas {
        locs[i] = loc;
        let ang = bas[i * 8 + 1] as usize;
        loc += nums[ang];
    }

    let mut emat = Array4::<f64>::zeros((nmat, nmat, nmat, nmat).f()); // 双电子积分矩阵
    let mut opt: *mut CINTOpt = std::ptr::null_mut();
    unsafe {
        cint2e_sph_optimizer(
            &mut opt,
            atm.as_ptr() as *mut c_int,
            natm as c_int,
            bas.as_ptr() as *mut c_int,
            nbas as c_int,
            env.as_ptr() as *mut c_double,
        );

        let mut dims = vec![0 as usize; nbas]; // 每个基函数的原子轨道数量
        for i in 0..nbas {
            dims[i] = CINTcgto_spheric(i as c_int, bas.as_ptr()) as usize;
        }

        for i in 0..nbas {
            let loci = locs[i];
            let dimi = dims[i];
            assert!(loci + dimi <= nmat, "数组越界");
            for j in 0..nbas {
                let locj = locs[j];
                let dimj = dims[j];
                assert!(locj + dimj <= nmat, "数组越界");
                for k in 0..nbas {
                    let lock = locs[k];
                    let dimk = dims[k];
                    assert!(lock + dimk <= nmat, "数组越界");
                    for l in 0..nbas {
                        let locl = locs[l];
                        let diml = dims[l];
                        assert!(locl + diml <= nmat, "数组越界");
                        let mut buf = vec![0.0f64; dimi * dimj * dimk * diml];
                        let mut shls = [i as c_int, j as c_int, k as c_int, l as c_int];
                        cint2e_sph(
                            buf.as_mut_ptr(),
                            shls.as_mut_ptr(),
                            atm.as_ptr() as *mut c_int,
                            natm as c_int,
                            bas.as_ptr() as *mut c_int,
                            nbas as c_int,
                            env.as_ptr() as *mut c_double,
                            opt,
                        );
                        let shape = (dimi, dimj, dimk, diml);
                        let buf_array = Array::from_shape_vec(shape.f(), buf).unwrap();
                        emat.slice_mut(s![
                            loci..loci + dimi,
                            locj..locj + dimj,
                            lock..lock + dimk,
                            locl..locl + diml
                        ])
                        .assign(&buf_array);
                    }
                }
            }
        }
        // cint2e 与 optimizer 各自定义了布局相同的 CINTOpt, 这里做一次裸指针转换
        CINTdel_optimizer(&mut opt as *mut *mut CINTOpt as *mut *mut ::rsint::src::optimizer::CINTOpt);
    }
    emat
}

/// 计算笛卡尔双电子积分
fn calc_int2e_car(nmat: usize, atm: &[i32], bas: &[i32], env: &[f64]) -> Array4<f64> {
    let natm = atm.len() / 6; // 原子数量
    let nbas = bas.len() / 8; // 基函数数量

    let mut locs = vec![0 as usize; nbas]; // 每个基函数在矩阵中的起始位置
    let nums: [usize; 7] = [1, 3, 6, 10, 15, 21, 28]; // 每个角动量对应的数量
    let mut loc = 0usize;
    for i in 0..nbas {
        locs[i] = loc;
        let ang = bas[i * 8 + 1] as usize;
        loc += nums[ang];
    }

    let mut emat = Array4::<f64>::zeros((nmat, nmat, nmat, nmat).f()); // 双电子积分矩阵
    let mut opt: *mut CINTOpt = std::ptr::null_mut();
    unsafe {
        cint2e_cart_optimizer(
            &mut opt,
            atm.as_ptr() as *mut c_int,
            natm as c_int,
            bas.as_ptr() as *mut c_int,
            nbas as c_int,
            env.as_ptr() as *mut c_double,
        );

        let mut dims = vec![0 as usize; nbas]; // 每个基函数的原子轨道数量
        for i in 0..nbas {
            dims[i] = CINTcgto_cart(i as c_int, bas.as_ptr()) as usize;
        }

        for i in 0..nbas {
            let loci = locs[i];
            let dimi = dims[i];
            assert!(loci + dimi <= nmat, "数组越界");
            for j in 0..nbas {
                let locj = locs[j];
                let dimj = dims[j];
                assert!(locj + dimj <= nmat, "数组越界");
                for k in 0..nbas {
                    let lock = locs[k];
                    let dimk = dims[k];
                    assert!(lock + dimk <= nmat, "数组越界");
                    for l in 0..nbas {
                        let locl = locs[l];
                        let diml = dims[l];
                        assert!(locl + diml <= nmat, "数组越界");
                        let mut buf = vec![0.0f64; dimi * dimj * dimk * diml];
                        let mut shls = [i as c_int, j as c_int, k as c_int, l as c_int];
                        cint2e_cart(
                            buf.as_mut_ptr(),
                            shls.as_mut_ptr(),
                            atm.as_ptr() as *mut c_int,
                            natm as c_int,
                            bas.as_ptr() as *mut c_int,
                            nbas as c_int,
                            env.as_ptr() as *mut c_double,
                            opt,
                        );
                        let shape = (dimi, dimj, dimk, diml);
                        let buf_array = Array::from_shape_vec(shape.f(), buf).unwrap();
                        emat.slice_mut(s![
                            loci..loci + dimi,
                            locj..locj + dimj,
                            lock..lock + dimk,
                            locl..locl + diml
                        ])
                        .assign(&buf_array);
                    }
                }
            }
        }
        CINTdel_optimizer(&mut opt as *mut *mut CINTOpt as *mut *mut ::rsint::src::optimizer::CINTOpt);
    }
    emat
}

fn len2idxs(len: usize) -> Vec<usize> {
    // 壳层内分量从 libcint 顺序到高斯顺序的置换表。
    //
    // 笛卡尔: libcint 按 `CINTcart_comp`(lx 递减、每个 lx 下 ly 递减)排列，
    // 目标顺序是本项目 `datas::btrans` 里 `CAR_*_LMNS` 的 Gaussian .fch 约定；
    // 球谐: libcint 按 m = -l..l 排列，高斯按 m = 0,+1,-1,+2,-2... 排列。
    // 长度本身就能区分 form(球谐 5/7/9/11，笛卡尔 6/10/15/21)。
    match len {
        1 => vec![0],
        3 => vec![0, 1, 2],
        // 球谐 d/f/g/h
        5 => vec![2, 3, 1, 4, 0],
        7 => vec![3, 4, 2, 5, 1, 6, 0],
        9 => vec![4, 5, 3, 6, 2, 7, 1, 8, 0],
        11 => vec![5, 6, 4, 7, 3, 8, 2, 9, 1, 10, 0],
        // 笛卡尔 d: [xx,xy,xz,yy,yz,zz] → [xx,yy,zz,xy,xz,yz]
        6 => vec![0, 3, 5, 1, 2, 4],
        // 笛卡尔 f: [xxy..] → [xxx,yyy,zzz,xyy,xxy,xxz,xzz,yzz,yyz,xyz]
        10 => vec![0, 6, 9, 3, 1, 2, 5, 8, 7, 4],
        // 笛卡尔 g/h: 高斯顺序恰好是 libcint 顺序的逆序
        15 => vec![14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
        21 => vec![20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0],
        _ => panic!("未知长度:{}", len),
    }
}
