#![allow(non_snake_case)]
use crate::base::Mole;
use crate::datas::rys_root;
use crate::datas::FACS;
use crate::maths;
use ndarray::Array2;
use ndarray::Array4;
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

    // 高斯乘积定理
    pub fn prod(&self, gto_j: &Gto) -> (f64, Gto) {
        let gto_i = self;
        let [xi, yi, zi] = gto_i.xyz;
        let [xj, yj, zj] = gto_j.xyz;
        let ai = gto_i.alp;
        let aj = gto_j.alp;
        let am = ai + aj;
        let xm = (ai * xi + aj * xj) / am;
        let ym = (ai * yi + aj * yj) / am;
        let zm = (ai * zi + aj * zj) / am;
        let xyz = [xm, ym, zm];
        let [li, mi, ni] = gto_i.lmn;
        let [lj, mj, nj] = gto_j.lmn;
        let lmn = [li + lj, mi + mj, ni + nj];
        let gto_m = Gto::new(lmn, am, xyz);
        let r2 = (xi - xj).powi(2) + (yi - yj).powi(2) + (zi - zj).powi(2);
        let K = f64::exp(-ai * aj / am * r2);
        (K, gto_m)
    }
}

struct Ry_val_1e {
    Ga: Gto,
    Gb: Gto,
    p: f64,
    u: f64,
}

impl Ry_val_1e {
    pub fn new(Ga: Gto, Gb: Gto) -> Self {
        let p = Ga.alp + Gb.alp; // 9.2.11
        let u = Ga.alp * Gb.alp / p; // 9.2.12
        Self { Ga, Gb, p, u }
    }

    pub fn eval_d(&self, i: isize, j: isize, t: isize, d: usize) -> f64 {
        let Xab = self.Ga.xyz[d] - self.Gb.xyz[d]; // 9.2.14
        let Xpa = -self.Gb.alp / self.p * Xab;
        let Xpb = self.Ga.alp / self.p * Xab;
        let Kab = f64::exp(-self.u * Xab.powi(2)); // 9.2.15
                                                   // 边界
        if i < 0 || j < 0 {
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

    // 公式 9.11.36
    pub fn hval_d(&self, t: isize, d: usize, xc: f64, s2: f64) -> f64 {
        let xa = self.Ga.xyz[d];
        let xb = self.Gb.xyz[d];
        let a = self.Ga.alp;
        let b = self.Gb.alp;
        let xp = (a * xa + b * xb) / self.p; // 9.2.13 电荷中心坐标
        let xpc = xp - xc;
        if t < 0 {
            return 0.0;
        }
        if t == 0 {
            return 1.0;
        }
        if t == 1 {
            return -2.0 * self.p * s2 * xpc;
        }
        let term = xpc * self.hval_d(t - 1, d, xc, s2) + (t - 1) as f64 * self.hval_d(t - 2, d, xc, s2);
        -2.0 * self.p * s2 * term
    }

    // 公式 9.11.43
    pub fn ival_d(&self, i: isize, j: isize, d: usize, xc: f64, s2: f64) -> f64 {
        let mut sum = 0.0;
        for t in 0..=(i + j) {
            let val = self.eval_d(i, j, t, d) * self.hval_d(t, d, xc, s2);
            sum += val
        }
        sum
    }

    // 核引积分，单个原子核 公式9.11.42
    pub fn vval_c(&self, rc: &[f64; 3]) -> f64 {
        let [xa, ya, za] = self.Ga.xyz;
        let [xb, yb, zb] = self.Gb.xyz;
        let a = self.Ga.alp;
        let b = self.Gb.alp;
        let xp = (a * xa + b * xb) / (a + b);
        let yp = (a * ya + b * yb) / (a + b);
        let zp = (a * za + b * zb) / (a + b);
        let rp = [xp, yp, zp]; // 重叠中心坐标
        let Rpc = maths::calc_dist(&rp, rc); // 重叠中心到原子核距离
        let alpha = self.p * Rpc.powi(2); // 9.11.34
        let mut gamma = 0;
        for d in 0..3 {
            let i = self.Ga.lmn[d] as isize;
            let j = self.Gb.lmn[d] as isize;
            gamma += i + j;
        }
        let gamma = (gamma as f64 / 2.0).floor() as usize + 1;
        let mut sum = 0.0;
        let roots = rys_root::rys_roots(gamma, alpha);
        for kappa in 0..gamma {
            let [s, w] = roots[kappa];
            let mut val = w;
            assert!(w > 0.0, "权重需为正值");
            for d in 0..3 {
                let i = self.Ga.lmn[d] as isize;
                let j = self.Gb.lmn[d] as isize;
                let xc = rc[d];
                val *= self.ival_d(i, j, d, xc, s.powi(2));
            }
            sum += val;
        }
        let vval = 2.0 * PI / self.p * sum; // 9.11.42
        vval
    }

    pub fn vval(&self, xyzs: &Vec<[f64; 3]>, crgs: &Vec<f64>) -> f64 {
        let mut vval = 0.0;
        let natm = crgs.len();
        for i in 0..natm {
            let xyz = xyzs[i];
            let crg = crgs[i];
            let val = -crg * self.vval_c(&xyz);
            vval += val;
        }
        vval
    }
}

struct Ry_val_2e {
    Ga: Gto,
    Gb: Gto,
    Gc: Gto,
    Gd: Gto,
    Gp: Gto,
    Gq: Gto,
    Kp: f64,
    Kq: f64,
}

impl Ry_val_2e {
    pub fn new(Ga: Gto, Gb: Gto, Gc: Gto, Gd: Gto) -> Self {
        let (Kp, Gp) = Ga.prod(&Gb);
        let (Kq, Gq) = Gc.prod(&Gd);
        Self { Ga, Gb, Gc, Gd, Gp, Gq, Kp, Kq }
    }

    // 公式 9.11.39
    pub fn hval_d(&self, t: isize, d: usize, s2: f64) -> f64 {
        let p = self.Gp.alp; // 9.9.28
        let q = self.Gq.alp; // 9.9.30
        let alp = p * q / (p + q);
        // let alp = 1.0 / (p + q);

        let Xpq = self.Gp.xyz[d] - self.Gq.xyz[d]; // 9.2.13 两个新的基函数的中心距离
        if t < 0 {
            return 0.0;
        }
        if t == 0 {
            return 1.0;
        }
        if t == 1 {
            return -2.0 * alp * s2 * Xpq;
        }
        // 9.11.39
        let term = Xpq * self.hval_d(t - 1, d, s2) + (t - 1) as f64 * self.hval_d(t - 2, d, s2);
        return -2.0 * alp * s2 * term;
    }

    pub fn eval_d(&self, i: isize, j: isize, t: isize, d: usize, c: &str) -> f64 {
        let (Ga, Gb, Gp) = match c {
            "ab" => (&self.Ga, &self.Gb, &self.Gp),
            "cd" => (&self.Gc, &self.Gd, &self.Gq),
            _ => panic!("c参数错误"),
        };
        let Xab = Ga.xyz[d] - Gb.xyz[d]; // 9.2.14

        let p = Gp.alp;

        let Xpa = Gp.xyz[d] - Ga.xyz[d];
        let Xpb = Gp.xyz[d] - Gb.xyz[d];

        let u = Ga.alp * Gb.alp / p;
        let Kab = f64::exp(-u * Xab.powi(2)); // 9.2.15
                                              // 边界
        if i < 0 || j < 0 {
            return 0.0;
        }
        // 9.5.5
        if t > i + j || t < 0 {
            return 0.0;
        }
        // 9.5.8
        if i == 0 && j == 0 && t == 0 {
            return if t == 0 { Kab } else { 0.0 };
        }
        // McMurchie-Davidson递归关系
        if t == 0 {
            // s轨道
            if i > j {
                // 9.5.15
                return Xpa * self.eval_d(i - 1, j, 0, d, c) + self.eval_d(i - 1, j, 1, d, c);
            } else {
                // 9.5.16
                return Xpb * self.eval_d(i, j - 1, 0, d, c) + self.eval_d(i, j - 1, 1, d, c);
            }
        }
        if t > 0 {
            // 9.5.17
            return (i as f64 * self.eval_d(i - 1, j, t - 1, d, c) + j as f64 * self.eval_d(i, j - 1, t - 1, d, c)) / (2.0 * p * t as f64);
        }
        panic!("不符条件");
    }

    // 公式 9.11.49 - 9.11.51
    pub fn ival_d(&self, i: isize, j: isize, k: isize, l: isize, d: usize, s2: f64) -> f64 {
        let mut sum = 0.0;
        for t in 0..=(i + j) {
            for u in 0..=(k + l) {
                let sign = if u % 2 == 0 { 1.0 } else { -1.0 };
                // 这里u是公式9.11.49中的tau
                sum += sign * self.eval_d(i, j, t, d, "ab") * self.eval_d(k, l, u, d, "cd") * self.hval_d(t + u, d, s2);
            }
        }
        sum
    }

    // 双电子积分 公式 9.11.48
    pub fn gval(&self) -> f64 {
        let p = self.Ga.alp + self.Gb.alp; // 9.9.28
        let q = self.Gc.alp + self.Gd.alp; // 9.9.30
        let mu = p * q / (p + q);
        let Rpq = maths::calc_dist(&self.Gp.xyz, &self.Gq.xyz);
        let alpha = mu * Rpq.powi(2);

        let mut L = 0;
        for d in 0..3 {
            let la = self.Ga.lmn[d] as isize;
            let lb = self.Gb.lmn[d] as isize;
            let lc = self.Gc.lmn[d] as isize;
            let ld = self.Gd.lmn[d] as isize;
            L += la + lb + lc + ld;
        }
        let gamma = (L as f64 / 2.0).floor() as usize + 1;
        let roots = rys_root::rys_roots(gamma, alpha);
        assert_eq!(roots.len(), gamma);
        let mut sum = 0.0;
        for kappa in 0..gamma {
            let [s, w] = roots[kappa];
            let mut mul = w;
            for d in 0..3 {
                let i = self.Ga.lmn[d] as isize;
                let j = self.Gb.lmn[d] as isize;
                let k = self.Gc.lmn[d] as isize;
                let l = self.Gd.lmn[d] as isize;
                mul *= self.ival_d(i, j, k, l, d, s.powi(2));
            }
            sum += mul
        }
        let k = 2.0 * PI.powf(2.5) / (p * q * f64::sqrt(p + q));
        k * sum
    }
}

pub struct Ry_mat {
    ato_alps: Vec<Vec<f64>>,   // 每个原子轨道的 指数列表
    ato_coes: Vec<Vec<f64>>,   // 每个原子轨道的 系数列表
    ato_lmns: Vec<[usize; 3]>, // 每个原子轨道 的角动量分量
    ato_atms: Vec<usize>,      // 每个原子轨道对应的原子索引
    atm_xyzs: Vec<[f64; 3]>,   // 每个原子的坐标
    atm_crgs: Vec<f64>,
}

impl Ry_mat {
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

    pub fn vmat(&self) -> Array2<f64> {
        let nato = self.ato_atms.len(); // 原子轨道数量
        let mut vmat = Array2::zeros((nato, nato)); // 核引矩阵
        for i in 0..nato {
            for j in 0..nato {
                let ni = self.ato_alps[i].len();
                let nj = self.ato_alps[j].len();

                for a in 0..ni {
                    for b in 0..nj {
                        let Ga = Gto::new(self.ato_lmns[i], self.ato_alps[i][a], self.atm_xyzs[self.ato_atms[i] as usize]);
                        let Gb = Gto::new(self.ato_lmns[j], self.ato_alps[j][b], self.atm_xyzs[self.ato_atms[j] as usize]);
                        let vval = Ry_val_1e::new(Ga, Gb).vval(&self.atm_xyzs, &self.atm_crgs);
                        vmat[(i, j)] += self.ato_coes[i][a] * self.ato_coes[j][b] * vval;
                    }
                }
            }
        }
        vmat
    }

    // 双电子库伦积分
    pub fn emat(&self) -> Array4<f64> {
        let nato = self.ato_atms.len(); // 原子轨道数量
        let mut emat = Array4::zeros((nato, nato, nato, nato));
        for i in 0..nato {
            for j in 0..nato {
                for k in 0..nato {
                    for l in 0..nato {
                        let ni = self.ato_alps[i].len();
                        let nj = self.ato_alps[j].len();
                        let nk = self.ato_alps[k].len();
                        let nl = self.ato_alps[l].len();
                        let mut sum = 0.0;
                        for a in 0..ni {
                            for b in 0..nj {
                                for c in 0..nk {
                                    for d in 0..nl {
                                        let Ga = Gto::new(self.ato_lmns[i], self.ato_alps[i][a], self.atm_xyzs[self.ato_atms[i] as usize]);
                                        let Gb = Gto::new(self.ato_lmns[j], self.ato_alps[j][b], self.atm_xyzs[self.ato_atms[j] as usize]);
                                        let Gc = Gto::new(self.ato_lmns[k], self.ato_alps[k][c], self.atm_xyzs[self.ato_atms[k] as usize]);
                                        let Gd = Gto::new(self.ato_lmns[l], self.ato_alps[l][d], self.atm_xyzs[self.ato_atms[l] as usize]);
                                        let gval = Ry_val_2e::new(Ga, Gb, Gc, Gd).gval();
                                        sum += gval * self.ato_coes[i][a] * self.ato_coes[j][b] * self.ato_coes[k][c] * self.ato_coes[l][d];
                                    }
                                }
                            }
                        }
                        emat[(i, j, k, l)] = sum;
                    }
                }
            }
        }
        emat
    }
}
