#![allow(non_snake_case)]

/**
 * 拓扑分析，我理解就是寻找分子中的一些点
*/
use crate::base::Mole;
use crate::maths;
use crate::props::gridprop::{self, GridVals, density};

use ndarray::{Array2, NewAxis};

use rayon::prelude::*;

/// 拓扑分析计算器，用于寻找分子中的临界点
pub struct Calculator<'m> {
    pub mole: &'m Mole,
}

impl<'m> Calculator<'m> {
    /// 创建一个新的 Calculator 实例
    pub fn new(mole: &'m Mole) -> Self {
        Calculator { mole }
    }

    // 提供一个初猜点，从该初猜点寻找临界点
    /// 从一个初猜点出发，使用 Newton-Raphson 方法寻找临界点
    pub fn get_point(
        &self,
        grid: &[f64; 3],                                              // 位置初猜
        func: &Box<dyn Fn(&[f64; 3]) -> GridVals + Send + Sync + 'm>, // 传入一个函数，这个函数要能计算二阶导
    ) -> ([f64; 3], f64, i32) {
        let mut grid = *grid; // 初猜位置
        let mut Hessian = Array2::zeros((3, 3));
        // let DHessian: Array2<f64> = Array2::eye(3) * 10.0;
        let mut rho = 0.0;
        for _i in 0..100 {
            // 最多迭代100次
            let rhos = func(&grid);
            rho = rhos.val0;
            Hessian = Array2::from_shape_vec((3, 3), rhos.val2s().to_vec()).expect("矩阵生成失败");
            // Hessian += &DHessian;
            let det = maths::det(&Hessian);

            if det.abs() < 1e-10 {
                break;
            }
            let [dx, dy, dz] = rhos.val1;
            // 梯度太小就结束
            if (dx.powi(2) + dy.powi(2) + dz.powi(2)).sqrt() < 1e-8 {
                break;
            }
            let Hessinv = maths::inv(&Hessian);

            let gradent = Array2::from_shape_vec((3, 1), rhos.val1.to_vec()).expect("矩阵生成失败");
            let delta = Hessinv.dot(&gradent).flatten().to_vec();
            let mut mx = delta[0];
            let mut my = delta[1];
            let mut mz = delta[2];
            grid[0] -= mx;
            grid[1] -= my;
            grid[2] -= mz;
        }
        let eig_vals = maths::eigh(&Hessian).0;
        let mut n_pos = 0;
        let mut n_neg = 0;
        for val in eig_vals.iter() {
            if *val > 1e-8 {
                n_pos += 1;
            }
            if *val < -1e-8 {
                n_neg += 1;
            }
        }
        let cpt = n_pos - n_neg; // 临界点类型
        (grid, rho, cpt)
    }

    /// 获取分子中的临界点初猜
    pub fn get_grids(&self) -> Vec<[f64; 3]> {
        let atoms = self.mole.atoms();
        let xyzs = atoms.xyzs();
        let mut grids = Vec::new();
        // 每个原子中心格点
        for xyz in xyzs.iter() {
            grids.push(xyz.clone());
        }
        let natm = xyzs.len();
        // 键中心格点
        for i in 0..natm - 1 {
            let [xi, yi, zi] = xyzs[i];
            for j in i + 1..natm {
                let [xj, yj, zj] = xyzs[j];
                if !atoms.is_bond(i, j) {
                    continue;
                }
                let cent = [(xi + xj) / 2.0, (yi + yj) / 2.0, (zi + zj) / 2.0];
                grids.push(cent);
            }
        }
        // 网格格点
        // let mut grid_caler = gridprop::CubeGrid::new();
        // let [p0, p1] = self.mole.border(2.0);
        // grid_caler.set_v1(&p0, &p1, 0.1, 0.0);
        // grids.extend(grid_caler.get().1);
        grids
    }

    pub fn get_func(&self, func: &str) -> Box<dyn Fn(&[f64; 3]) -> GridVals + Send + Sync + 'm> {
        match func {
            "mol_rho_dm" => {
                let rho_caler = gridprop::density::Calculator::new(self.mole);
                Box::new(move |grid| rho_caler.mol_rho_dm(grid, 2))
            }
            "pi_rho_pocv" => {
                let pi_caler = gridprop::density::Calculator::new(self.mole);
                Box::new(move |grid| pi_caler.pi_rho_pocv(grid, 2))
            }
            "pi_rho_mocv" => {
                let pi_caler = gridprop::density::Calculator::new(self.mole);
                Box::new(move |grid| pi_caler.pi_rho_mocv(grid, 2))
            }
            _ => panic!("未知函数:{}", func),
        }
    }

    pub fn analysis(&self, func: &str) -> Vec<([f64; 3], f64, i32)> {
        let func = self.get_func(func);
        let grids = self.get_grids();
        log::debug!("grids: {}", grids.len());

        // 并行计算所有候选临界点（密度过滤 + 临界点搜索）
        let mut candidates: Vec<([f64; 3], f64, i32)> = grids
            .par_iter()
            .filter_map(|grid| {
                let val = func(grid);
                if val.val0 < 1e-1 {
                    return None;
                }
                let res = self.get_point(grid, &func);
                if res.1 < 1e-1 {
                    return None;
                }
                Some(res)
            })
            .collect();

        // 串行去重（维持原有的邻居判重逻辑）
        let mut cpts = vec![];
        for res in candidates.into_iter() {
            let mut right = true;
            for (pos, val, _) in cpts.iter() {
                let dist = maths::vec_dist(&res.0, pos);
                if dist < 1e-1 {
                    right = false;
                    break;
                }
            }
            if right {
                cpts.push(res);
            }
        }
        cpts
    }
}
