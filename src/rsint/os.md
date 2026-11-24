// The Obara-Saika scheme
#![allow(non_snake_case)]

use crate::{base::Mole, datas::FACS};
use ndarray::Array2;
use std::f64::consts::PI;

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
struct OS_val {
    Ga: Gto,
    Gb: Gto,
    p: f64,
    u: f64,
}

impl OS_val {
    pub fn new(Ga: Gto, Gb: Gto) -> Self {
        let p = Ga.alp + Gb.alp;
        let u = Ga.alp * Gb.alp / p;
        Self { Ga, Gb, p, u }
    }

    pub fn sval_d(&self, i: isize, j: isize, d: usize) -> f64 {
        if i < 0 || j < 0 {
            return 0.0;
        }
        let Xab = self.Ga.xyz[d] - self.Gb.xyz[d];
        let Xpa = -self.Gb.alp / self.p * Xab;
        let Xpb = self.Ga.alp / self.p * Xab;
        let fi = i as f64;
        let fj = j as f64;
        if i == 0 && j == 0 {
            return f64::sqrt(PI / self.p) * f64::exp(-self.u * Xab.powi(2)); // 9.3.10
        }
        if i > 0 {
            // 9.3.8
            return Xpa * self.sval_d(i - 1, j, d) + ((fi - 1.0) * self.sval_d(i - 2, j, d) + fj * self.sval_d(i - 1, j - 1, d)) / (2.0 * self.p);
        } else {
            // 9.3.9
            return Xpb * self.sval_d(i, j - 1, d) + (fi * self.sval_d(i - 1, j - 1, d) + (fj - 1.0) * self.sval_d(i, j - 2, d)) / (2.0 * self.p);
        }
    }

    pub fn tval_d(&self, i: isize, j: isize, d: usize) -> f64 {
        if i < 0 || j < 0 {
            return 0.0;
        }
        let Xab = self.Ga.xyz[d] - self.Gb.xyz[d];
        let Xpa = -self.Gb.alp / self.p * Xab;
        let Xpb = self.Ga.alp / self.p * Xab;
        let a = self.Ga.alp;
        let b = self.Gb.alp;
        if i == 0 && j == 0 {
            // 9.3.42
            let S00 = self.sval_d(0, 0, d);
            // let term = a - 2.0 * a.powi(2) * (Xab.powi(2) + 2.0 / self.p);
            let term = a - 2.0 * a.powi(2) * (Xpa.powi(2) + 1.0 / (2.0 * self.p));
            let val = term * S00;
            return val;
        }
        if i > 0 {
            // 9.3.40
            let i = i - 1;
            let tzz = self.tval_d(i, j, d);
            let tnz = self.tval_d(i - 1, j, d);
            let tzn = self.tval_d(i, j - 1, d);
            let spz = self.sval_d(i + 1, j, d);
            let snz = self.sval_d(i - 1, j, d);
            let val = Xpa * tzz + (i as f64 * tnz + j as f64 * tzn) / (2.0 * self.p) + (b / self.p) * (2.0 * a * spz - i as f64 * snz);
            return val;
        } else {
            // 9.3.41
            let j = j - 1;
            let tzz = self.tval_d(i, j, d);
            let tnz = self.tval_d(i - 1, j, d);
            let tzn = self.tval_d(i, j - 1, d);
            let szp = self.sval_d(i, j + 1, d);
            let szn = self.sval_d(i, j - 1, d);
            let val = Xpb * tzz + (i as f64 * tnz + j as f64 * tzn) / (2.0 * self.p) + (a / self.p) * (2.0 * b * szp - j as f64 * szn);
            return val;
        }
    }

    pub fn sval(&self) -> f64 {
        let mut sval = 1.0;
        for d in 0..3 {
            let i = self.Ga.lmn[d] as isize;
            let j = self.Gb.lmn[d] as isize;
            sval *= self.sval_d(i, j, d);
        }
        sval
    }

    pub fn tval(&self) -> f64 {
        let i = self.Ga.lmn[0] as isize;
        let k = self.Ga.lmn[1] as isize;
        let m = self.Ga.lmn[2] as isize;
        let j = self.Gb.lmn[0] as isize;
        let l = self.Gb.lmn[1] as isize;
        let n = self.Gb.lmn[2] as isize;
        let V0 = self.tval_d(i, j, 0) * self.sval_d(k, l, 1) * self.sval_d(m, n, 2);
        let V1 = self.sval_d(i, j, 0) * self.tval_d(k, l, 1) * self.sval_d(m, n, 2);
        let V2 = self.sval_d(i, j, 0) * self.sval_d(k, l, 1) * self.tval_d(m, n, 2);

        let val = V0 + V1 + V2;
        val
    }
}

pub struct OS_mat {
    ato_alps: Vec<Vec<f64>>,   // 每个原子轨道的 指数列表
    ato_coes: Vec<Vec<f64>>,   // 每个原子轨道的 系数列表
    ato_lmns: Vec<[usize; 3]>, // 每个原子轨道 的角动量分量
    ato_atms: Vec<usize>,      // 每个原子轨道对应的原子索引
    atm_xyzs: Vec<[f64; 3]>,   // 每个原子的坐标
    atm_crgs: Vec<f64>,
}

impl OS_mat {
    pub fn new(mole: &Mole) -> Self {
        let map_ato = mole.basis().map_ato(mole);
        Self {
            ato_alps: map_ato.alps.clone(),
            ato_coes: map_ato.coes_norm(),
            ato_lmns: map_ato.lmns.clone(),
            ato_atms: map_ato.atms.clone(),
            atm_xyzs: mole.atoms().xyzs().clone(),
            atm_crgs: mole.atoms().crgs().clone(),
        }
    }

    pub fn smat(&self) -> Array2<f64> {
        let nato = self.ato_atms.len(); // 原子轨道数量
        let mut smat = Array2::zeros((nato, nato)); // 重叠矩阵
        for i in 0..nato {
            let ni = self.ato_alps[i].len();
            for j in 0..nato {
                let nj = self.ato_alps[j].len();

                for k in 0..ni {
                    for l in 0..nj {
                        let Ga = Gto::new(self.ato_lmns[i], self.ato_alps[i][k], self.atm_xyzs[self.ato_atms[i] as usize]);
                        let Gb = Gto::new(self.ato_lmns[j], self.ato_alps[j][l], self.atm_xyzs[self.ato_atms[j] as usize]);
                        let sval = OS_val::new(Ga, Gb).sval();
                        smat[(i, j)] += self.ato_coes[i][k] * self.ato_coes[j][l] * sval;
                    }
                }
            }
        }
        smat
    }

    pub fn tmat(&self) -> Array2<f64> {
        let nato = self.ato_atms.len(); // 原子轨道数量
        let mut tmat = Array2::zeros((nato, nato)); // 重叠矩阵
        for i in 0..nato {
            let ni = self.ato_alps[i].len();
            for j in 0..nato {
                let nj = self.ato_alps[j].len();

                for k in 0..ni {
                    for l in 0..nj {
                        let Ga = Gto::new(self.ato_lmns[i], self.ato_alps[i][k], self.atm_xyzs[self.ato_atms[i] as usize]);
                        let Gb = Gto::new(self.ato_lmns[j], self.ato_alps[j][l], self.atm_xyzs[self.ato_atms[j] as usize]);
                        let tval = OS_val::new(Ga, Gb).tval();
                        tmat[(i, j)] += self.ato_coes[i][k] * self.ato_coes[j][l] * tval;
                    }
                }
            }
        }
        tmat
    }
}
