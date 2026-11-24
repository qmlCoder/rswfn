#![allow(non_snake_case)]

use crate::{datas::FACS, maths};
use ndarray::Array2;
use special::Error;
use std::f64::{self, consts::PI};

pub struct Gto {
    lmn: [usize; 3],
    alp: f64,
    xyz: [f64; 3],
}

impl Gto {
    pub fn new(lmn: [usize; 3], alp: f64, xyz: [f64; 3]) -> Self {
        Self { lmn, alp, xyz }
    }
}

struct DM_val {
    Ga: Gto,
    Gb: Gto,
    p: f64,
    u: f64,
}

impl DM_val {
    pub fn new(Ga: Gto, Gb: Gto) -> Self {
        let p = Ga.alp + Gb.alp;
        let u = Ga.alp * Gb.alp / p;
        Self { Ga, Gb, p, u }
    }

    // McMurchie-Davidson中的E_t^{ij}系数
    pub fn eval_d(&self, i: isize, j: isize, t: isize, d: usize) -> f64 {
        let Xab = self.Ga.xyz[d] - self.Gb.xyz[d]; // 9.2.14
        let Xpa = -self.Gb.alp / self.p * Xab;
        let Xpb = self.Ga.alp / self.p * Xab;
        let Kab = f64::exp(-self.u * Xab.powi(2)); // 9.2.15
                                                   // 边界
        if i < 0 || j < 0 || t < 0 {
            return 0.0;
        }
        // 9.5.5
        if t > i + j || t < 0 {
            return 0.0;
        }
        // 9.5.8
        if i == 0 && j == 0 {
            return if t == 0 { Kab } else { 0.0 };
        }
        // McMurchie-Davidson递归关系
        if t == 0 {
            // s轨道
            if i > j {
                // 9.5.15
                return Xpa * self.eval_d(i - 1, j, 0, d) + self.eval_d(i - 1, j, 1, d);
            } else {
                // 9.5.16
                return Xpb * self.eval_d(i, j - 1, 0, d) + self.eval_d(i, j - 1, 1, d);
            }
        }
        if t > 0 {
            // 9.5.17
            return (i as f64 * self.eval_d(i - 1, j, t - 1, d) + j as f64 * self.eval_d(i, j - 1, t - 1, d)) / (2.0 * self.p * t as f64);
        }
        panic!("不符条件");
    }

    // McMurchie-Davidson中的M_t^{e}积分
    pub fn mval(&self, e: isize, t: isize) -> f64 {
        let Xpc = 0.0;
        if t < 0 || t > e {
            return 0.0;
        }
        if e == 0 {
            if t == 0 {
                return f64::sqrt(PI / self.p);
            } else {
                return 0.0;
            }
        } else {
            return t as f64 * self.mval(e - 1, t - 1) + Xpc * self.mval(e - 1, t) + self.mval(e - 1, t + 1) / (2.0 * self.p);
        }
    }
    pub fn sval(&self) -> f64 {
        let mut sval = 1.0;
        for d in 0..3 {
            let i = self.Ga.lmn[d] as isize;
            let j = self.Gb.lmn[d] as isize;
            let mut ival = 0.0; // 某个方向的重叠积分
            for t in 0..=(i + j) {
                let eval = self.eval_d(i, j, t, d);
                let mval = self.mval(0, t); // !重叠积分e为0
                ival += eval * mval;
            }
            sval *= ival;
        }
        sval
    }

    // 9.9.18-9.9.20
    pub fn rval_d(&self, t: usize, u: usize, v: usize, n: usize, Rpc: &[f64; 3]) {
        if t == 0 && u == 0 && v == 0 && n == 0 {}
    }

    // 9.9.32
    pub fn jval(&self, rc: &[f64; 3]) -> f64 {
        let a = self.Ga.alp;
        let b = self.Gb.alp;
        let p = a + b;
        let [xa, ya, za] = self.Ga.xyz;
        let [xb, yb, zb] = self.Gb.xyz;
        let xp = (a * xa + b * xb) / p;
        let yp = (a * ya + b * yb) / p;
        let zp = (a * za + b * zb) / p;
        let rp = [xp, yp, zp];
        let [xc, yc, zc] = *rc;
        let Rpc = [xp - xc, yp - yc, zp - zc];

        // let Rpc = maths::calc_dist(&rp, rc);
        0.0
    }
}

pub struct DM_mat {
    alps: Vec<Vec<f64>>,   // 每个原子轨道的 指数列表
    coes: Vec<Vec<f64>>,   // 每个原子轨道的 系数列表
    lmns: Vec<[usize; 3]>, // 每个原子轨道 的角动量分量
    xyzs: Vec<[f64; 3]>,   // 每个原子轨道的坐标
}

impl DM_mat {
    pub fn new(alps: Vec<Vec<f64>>, coes: Vec<Vec<f64>>, lmns: Vec<[usize; 3]>, xyzs: Vec<[f64; 3]>) -> Self {
        let mut coes = coes.clone();
        for i in 0..coes.len() {
            let [l, m, n] = lmns[i];
            for j in 0..coes[i].len() {
                let fac = FACS[l as usize] * FACS[m as usize] * FACS[n as usize];
                let ang = l + m + n;
                let alp = alps[i][j];
                let nm = (2.0 * alp / PI).powf(0.75) * f64::sqrt((4.0 * alp).powi(ang as i32) / fac);
                coes[i][j] = nm * coes[i][j];
            }
        }
        Self { alps, coes, lmns, xyzs }
    }

    pub fn smat(&self) -> Array2<f64> {
        let nato = self.xyzs.len(); // 原子轨道数量
        let mut smat = Array2::zeros((nato, nato)); // 重叠矩阵
        for i in 0..nato {
            let ni = self.alps[i].len();
            for j in 0..nato {
                let nj = self.alps[j].len();

                for k in 0..ni {
                    for l in 0..nj {
                        let Ga = Gto::new(self.lmns[i], self.alps[i][k], self.xyzs[i]);
                        let Gb = Gto::new(self.lmns[j], self.alps[j][l], self.xyzs[j]);
                        let sval = DM_val::new(Ga, Gb).sval();
                        smat[(i, j)] += self.coes[i][k] * self.coes[j][l] * sval;
                    }
                }
            }
        }
        smat
    }
}

// 误差函数
pub fn erf(x: f64) -> f64 {
    if x.abs() < 0.5 {
        let mut val = 0.0;
        for k in 0..6_isize {
            // 泰勒展开
            val += (-1.0f64).powi(k as i32) * x.powi(2 * k as i32 + 1) / (maths::fac(k) * (2 * k + 1)) as f64;
        }
        return 2.0 * val / PI.sqrt();
    }
    if x.abs() > 4.0 {
        return 1.0 - erfc(x);
    }
    panic!("还不能计算的区间");
}

pub fn erfc(x: f64) -> f64 {
    if x.abs() < 0.5 {
        return 1.0 - erf(x);
    }
    if x.abs() > 4.0 {
        let mut val = 1.0;
        for k in 1..=2_isize {
            val += (-1.0_f64).powi(k as i32) * maths::fac2(2 * k - 1) as f64 / (2.0 * x.powi(2)).powi(k as i32) as f64;
        }
        val *= (-x.powi(2)).exp() / (x * x.sqrt());
        return val;
    }
    panic!("还不能计算的区间");
}
