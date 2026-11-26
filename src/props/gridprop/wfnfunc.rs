use cached::proc_macro::cached;
use ndarray::Array2;

use crate::base::Mole;
use crate::datas::btrans::sym2lmn;
use crate::datas::consts::BOHR;
use crate::datas::{FACS, elems};
use crate::props::gridprop::GridVals;
use crate::maths;
use std::f64::consts::PI;

#[cached(size = 1024)]
fn get_cut(nuc: usize, nlm: (usize, usize, isize)) -> f64 {
    let (n, l, m) = nlm;
    match (nuc, n, l, m) {
        (2, 2, 1, 0) => 11.583131071979427,
        (3, 2, 1, 0) => 7.946768369583903,
        (4, 2, 1, 0) => 6.07950166730683,
        (5, 2, 1, 0) => 4.937649206138324,
        (6, 2, 1, 0) => 4.165094792995542,
        (7, 2, 1, 0) => 3.6065789556558063,
        (8, 2, 1, 0) => 3.1834089304669204,
        (9, 2, 1, 0) => 2.8513701142465595,
        (10, 2, 1, 0) => 2.58367637157943,
        _ => 15.0,
    }
}

const COEF1: f64 = 0.5641895835477563;
const COEF2: f64 = 0.09973557010035818;
const COEF3: f64 = 0.0040214198510528385;
/// slater type function 结构化学基础 第五版 周公度 P30
pub fn stf(xyz: &[f64; 3], nuc: usize, Z: f64, nlm: &(usize, usize, isize)) -> f64 {
    let r = maths::vec_len(&xyz);
    // let cut = get_cut(nuc, *nlm);
    // if r > 5.0 {
    //     return 0.0;
    // }
    let [x, y, z] = xyz;
    let a0 = BOHR;
    let sigma = Z / a0 * r;
    match nlm {
        (1, 0, 0) => COEF1 * (Z / a0).powf(1.5) * f64::exp(-sigma),
        (2, 0, 0) => COEF1 * (Z / a0).powf(1.5) * (2.0 - sigma) * f64::exp(-sigma / 2.0),
        (2, 1, 0) => COEF2 * (Z / a0).powf(1.5) * sigma * f64::exp(-sigma / 2.0) * (z / r),
        (3, 1, 0) => COEF3 * (Z / a0).powf(1.5) * (6.0 - sigma) * sigma * f64::exp(-sigma / 3.0) * (z / r),
        _ => panic!("暂不支持的量子数:{:?}", nlm),
    }
}

/// 波函数计算器，计算原子轨道和分子轨道在空间中的波函数值及其导数
pub struct Calculator<'m> {
    pub mole: &'m Mole,
    pub atms: &'m Vec<usize>,
}

impl<'m> Calculator<'m> {
    /// 创建一个新的 Calculator 实例
    pub fn new(mole: &'m Mole) -> Self {
        let atms = mole.atoms().idxs();
        Calculator { mole: mole, atms: atms }
    }

    // 计算一个点处的所有原子轨道的波函数
    // xyzs: 原子轨道坐标，coes和alps：每个基函数的系数和指数，每个基函数的数量是不确定的
    /// 计算指定点处所有原子轨道的波函数及其一阶、二阶导数
    pub fn ato_wfn(
        &self,
        grid: &[f64; 3], // 空间中任意一点的坐标
        level: usize,
    ) -> Vec<GridVals> {
        let atos = &self.mole.basis().atos().cars;
        let nato = atos.len();
        let mut vals = vec![GridVals::non(); nato];
        let atoms = &self.mole.atoms;
        for i in 0..nato {
            // 对每个原子进行循环
            let atm = atos[i].atm;
            if self.atms.contains(&atm) {
                let [ax, ay, az] = atoms.get(atm).xyz;
                let px = grid[0] - ax;
                let py = grid[1] - ay;
                let pz = grid[2] - az;
                let coes = atos[i].coes_norm();
                let alps = &atos[i].alps;
                let xyz = [px, py, pz];
                let lmn = atos[i].lmn();
                let (val0, val1, val2) = cgf(&xyz, lmn, &alps, &coes, level);
                vals[i].set(val0, val1, val2);
            }
        }
        vals
    }

    // 计算一个点处的某个分子轨道的波函数
    /// 计算指定点处某个分子轨道的波函数及其一阶、二阶导数
    pub fn obt_wfn(
        &self,
        obt: usize, // 分子轨道索引
        grid: &[f64; 3],
        level: usize,
        // wfns: &(Vec<f64>, Vec<[f64; 3]>, Vec<[[f64; 3]; 3]>),
    ) -> GridVals {
        let cmat_car = self.mole.coefs().get_cmat("car", self.mole);
        let nato = cmat_car.nrows(); // 原子轨道数量
        let wfns = self.ato_wfn(grid, level);
        let mut obt_wfn0 = 0.0;
        let mut obt_wfn1 = [0.0; 3];
        let mut obt_wfn2 = [[0.0; 3]; 3];
        for ato in 0..nato {
            // 循环所有的原子轨道
            let coef = cmat_car[(ato, obt)];
            // 如果系数为0，则跳过
            if coef == 0.0 {
                continue;
            };
            obt_wfn0 += coef * wfns[ato].val0;

            if level >= 1 {
                for i in 0..3 {
                    obt_wfn1[i] += coef * wfns[ato].val1[i];
                }
            };

            if level >= 2 {
                for i in 0..3 {
                    for j in 0..3 {
                        obt_wfn2[i][j] += coef * wfns[ato].val2[i][j];
                    }
                }
            };
        }
        GridVals::new(obt_wfn0, obt_wfn1, obt_wfn2)
    }
}

/// 计算指定点处所有原子轨道的波函数（自由函数版本）
pub fn ato_wfn(
    mole: &Mole,     // 分子信息，提供基组与原子结构
    grid: &[f64; 3], // 空间中任意一点的坐标
    level: usize,
) -> (Vec<f64>, Vec<[f64; 3]>, Vec<[[f64; 3]; 3]>) {
    let atos = &mole.basis().atos().cars;
    let atms = mole.atoms().idxs();
    let atoms = mole.atoms();
    let nato = atos.len();
    let mut wfn0 = vec![0.0; nato];
    let mut wfn1 = vec![[0.0; 3]; nato];
    let mut wfn2 = vec![[[0.0; 3]; 3]; nato];
    for i in 0..nato {
        // 对每个原子进行循环
        let atm = &atos[i].atm;
        if !atms.contains(atm) {
            continue;
        }
        let [ax, ay, az] = atoms.get(*atm).xyz;
        let px = grid[0] - ax;
        let py = grid[1] - ay;
        let pz = grid[2] - az;
        let coes = atos[i].coes_norm(); // 使用提前归一化的系数
        let alps = &atos[i].alps;
        let xyz = [px, py, pz];
        let lmn = atos[i].lmn();
        let (val0, val1, val2) = cgf(&xyz, lmn, &alps, &coes, level);
        wfn0[i] = val0;
        wfn1[i] = val1;
        wfn2[i] = val2;
    }
    (wfn0, wfn1, wfn2)
}

#[rustfmt::skip]
/// 计算高斯基函数（Gaussian Type Function）及其一阶、二阶导数
#[inline(always)]
pub fn gtf(
    xyz: &[f64; 3], // 以原子为中心的格点坐标
    lmn: &[usize; 3],
    alp: f64,
    level: usize,
) -> (f64, [f64; 3], [[f64; 3]; 3]) {
    let [x, y, z] = *xyz;
    let l = lmn[0] as i32;
    let m = lmn[1] as i32;
    let n = lmn[2] as i32;
    let lf = l as f64;
    let mf = m as f64;
    let nf = n as f64;
    let r2 = x*x + y*y + z*z;
    let exp_x = -alp * r2; // α.r^2
    let mut wfn0 = 0.0;
    let mut wfn1 = [0.0; 3];
    let mut wfn2 = [[0.0; 3]; 3];
    // 这里是减少计算量
    if exp_x < -20.0{
        return (wfn0, wfn1, wfn2);
    }
    // let nm=1.0;
    let exp_v = f64::exp(exp_x);
    let xyz_v=x.powi(l)*y.powi(m)*z.powi(n);

    wfn0 = xyz_v * exp_v;

    if level >= 1 {
        let tx = if l >= 1 { lf * x.powi(l - 1) } else { 0.0 };
        let ty = if m >= 1 { mf * y.powi(m - 1) } else { 0.0 };
        let tz = if n >= 1 { nf * z.powi(n - 1) } else { 0.0 };

        let dx = (tx - 2.0 * alp * x.powi(l + 1)) * y.powi(m) * z.powi(n) * exp_v;
        let dy = (ty - 2.0 * alp * y.powi(m + 1)) * x.powi(l) * z.powi(n) * exp_v;
        let dz = (tz - 2.0 * alp * z.powi(n + 1)) * x.powi(l) * y.powi(m) * exp_v;

        wfn1 = [dx, dy, dz];

        if level >= 2 {
            let txx = if l >= 2 {lf * (lf - 1.0) * x.powi(l - 2)} else {0.0};
            let tyy = if m >= 2 {mf * (mf - 1.0) * y.powi(m - 2)} else {0.0};
            let tzz = if n >= 2 {nf * (nf - 1.0) * z.powi(n - 2)} else {0.0};

            let dxx = y.powi(m) * z.powi(n) * exp_v * (txx + 2.0 * alp * x.powi(l) * (-2.0 * lf + 2.0 * alp * x.powi(2) - 1.0));
            let dyy = x.powi(l) * z.powi(n) * exp_v * (tyy + 2.0 * alp * y.powi(m) * (-2.0 * mf + 2.0 * alp * y.powi(2) - 1.0));
            let dzz = x.powi(l) * y.powi(m) * exp_v * (tzz + 2.0 * alp * z.powi(n) * (-2.0 * nf + 2.0 * alp * z.powi(2) - 1.0));

            let ttx = tx - 2.0 * alp * x.powi(l+1);
            let tty = ty - 2.0 * alp * y.powi(m+1);
            let ttz = tz - 2.0 * alp * z.powi(n+1);

            let dxy = z.powi(n) * exp_v * ttx * tty;
            let dyz = x.powi(l) * exp_v * tty * ttz;
            let dxz = y.powi(m) * exp_v * ttx * ttz;

            wfn2 = [
                [dxx, dxy, dxz],
                [dxy, dyy, dyz],
                [dxz, dyz, dzz]
            ];
        };

    };
    (wfn0, wfn1, wfn2)
}

/// 计算收缩基函数（原子轨道），即基函数的线性组合
#[inline(always)]
pub fn cgf(
    xyz: &[f64; 3],
    lmn: &[usize; 3],
    alps: &Vec<f64>,
    coes: &Vec<f64>,
    level: usize,
) -> (f64, [f64; 3], [[f64; 3]; 3]) {
    let mut wfn0 = 0.0;
    let mut wfn1 = [0.0; 3];
    let mut wfn2 = [[0.0; 3]; 3];
    for i in 0..coes.len() {
        let alp = alps[i];
        let coe = coes[i];
        let (val0, val1, val2) = gtf(&xyz, &lmn, alp, level);
        wfn0 += coe * val0;
        if level >= 1 {
            for j in 0..3 {
                wfn1[j] += coe * val1[j];
            }
        };

        if level >= 2 {
            for j in 0..3 {
                for k in 0..3 {
                    wfn2[j][k] += coe * val2[j][k];
                }
            }
        };
    }
    (wfn0, wfn1, wfn2)
}
