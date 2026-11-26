#![allow(non_snake_case)]

use crate::base::Mole;

use crate::props::gridprop::{density, dftgrid};
use crate::maths;

/// 碎片间静电相互作用能计算器
pub struct Calculator<'m> {
    pub mole: &'m Mole,
}

impl<'m> Calculator<'m> {
    /// 创建一个新的 Calculator 实例
    pub fn new(mole: &'m Mole) -> Self {
        Calculator { mole }
    }

    // 计算原子(片段)间静电相互作用能 electrostatic interaction energy between atoms
    /// 计算两个原子（片段）之间的静电相互作用能
    ///
    /// 返回 (电子-电子排斥, 核-电子吸引, 核-核排斥)
    pub fn EIEBA(&self, fragA: &Vec<usize>, fragB: &Vec<usize>) -> (f64, f64, f64) {
        let gridCaler = dftgrid::Calculator::new(&self.mole);
        let densCaler = density::Calculator::new(&self.mole);
        let (gridsA, weitsA) = gridCaler.frag_grids(&fragA);
        let (gridsB, weitsB) = gridCaler.frag_grids(&fragA);
        let atms = self.mole.atoms().idxs();
        let rhosA: Vec<f64> = gridsA.iter().map(|grid| densCaler.mol_rho_dm(grid, 0).val0).collect();
        let rhosB: Vec<f64> = gridsB.iter().map(|grid| densCaler.mol_rho_dm(grid, 0).val0).collect();
        let mut eng_ee = 0.0; // 电子-电子排斥
        let mut eng_en = 0.0; // 核-电子吸引
        let mut eng_nn = 0.0;

        for i in 0..rhosA.len() {
            // 电子-电子排斥
            for j in 0..rhosB.len() {
                let dist_ee = maths::calc_dist(&gridsA[i], &gridsB[j]); // 两个格点之间的距离
                if dist_ee < 1e-6 {
                    continue;
                }
                eng_ee += rhosA[i] * weitsA[i] * rhosB[j] * weitsB[j] / dist_ee
            }
        }

        let atoms = self.mole.atoms();
        // A核-B电子 吸引
        for g in 0..rhosB.len() {
            for a in fragA {
                let [x, y, z] = atoms.get(*a).xyz;
                let c = atoms.get(*a).nuc as f64;
                let dist_ne = maths::calc_dist(&[x, y, z], &gridsB[g]);
                if dist_ne < 1e-6 {
                    continue;
                }
                eng_en -= rhosB[g] * weitsB[g] * c / dist_ne;
            }
        }

        // B核-A电子 吸引
        for g in 0..rhosA.len() {
            for a in fragB {
                let [x, y, z] = atoms.get(*a).xyz;
                let c = atoms.get(*a).nuc as f64;
                let dist_ne = maths::calc_dist(&[x, y, z], &gridsA[g]);
                if dist_ne < 1e-6 {
                    continue;
                }
                eng_en -= rhosA[g] * weitsA[g] * c / dist_ne;
            }
        }

        // A核-B核 排斥
        for i in fragA {
            let [xi, yi, zi] = atoms.get(*i).xyz;
            let ci = atoms.get(*i).nuc as f64;
            for j in fragB {
                let [xj, yj, zj] = atoms.get(*j).xyz;
                let cj = atoms.get(*j).nuc as f64;
                let dist_ij = maths::calc_dist(&[xi, yi, zi], &[xj, yj, zj]);
                eng_nn += ci * cj / dist_ij;
            }
        }

        (eng_ee, eng_en, eng_nn)
    }
}
