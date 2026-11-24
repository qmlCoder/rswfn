use ndarray::Array2;
use rayon::prelude::*;
use std::f64::{EPSILON, consts::PI};

use crate::base::Vec3;
use crate::march;

pub mod density;
pub mod dftgrid;
pub mod topolog;
pub mod wfnfunc;

/// 一维线性格点
pub struct LineGrid {
    pub shape: [usize; 1],
    pub grids: Vec<[f64; 3]>,
}

impl LineGrid {
    /// 创建一个空的 LineGrid 实例
    pub fn new() -> LineGrid {
        LineGrid {
            shape: [0],
            grids: vec![],
        }
    }

    /// 根据起点、终点和步长生成一维线性格点
    pub fn set(&mut self, p0: &Vec3, p1: &Vec3, step: f64) {
        let vector = p1 - p0;
        let length = vector.length();
        let count = (length / step) as usize;
        let mut grids: Vec<[f64; 3]> = Vec::with_capacity(count as usize * 3);
        let vector = vector.normalize();
        for i in 0..count {
            let point = p0 + (&vector * i as f64) * step;
            let [px, py, pz] = point.into_arr();
            grids.push([px, py, pz]);
        }
        self.shape = [count];
        self.grids = grids;
    }

    /// 获取格点形状和数据
    pub fn get(&self) -> ([usize; 1], &Vec<[f64; 3]>) {
        (self.shape, &self.grids)
    }
}

/// 二维矩形格点
pub struct RectGrid {
    pub shape: [usize; 2],
    pub grids: Vec<[f64; 3]>,
}

impl RectGrid {
    /// 创建一个空的 RectGrid 实例
    pub fn new() -> Self {
        RectGrid {
            shape: [0, 0],
            grids: vec![],
        }
    }

    /// 根据中心点、法向量、方向向量、尺寸和步长生成平面矩形格点
    pub fn set_v1(&mut self, center: &Vec3, normal: &Vec3, vector: &Vec3, size: [f64; 2], step: [f64; 2]) {
        let vz = normal.normalize();
        let vx = vector.normalize();
        let shape = [
            (size[0] / step[0] + EPSILON).floor() as usize + 1,
            (size[1] / step[1] + EPSILON).floor() as usize + 1,
        ];

        let vy = vx.cross(&vz).normalize();
        let vx = vz.cross(&vy).normalize();

        let p0 = center - (&vx * shape[0] as f64 * step[0] / 2.0 + &vy * shape[1] as f64 * step[1] / 2.0);
        let count = shape[0] * shape[1];
        let mut grids = Vec::with_capacity(count as usize * 3);

        for i in 0..shape[0] {
            for j in 0..shape[1] {
                let point = &p0 + &vx * i as f64 * step[0] + &vy * j as f64 * step[1];
                let [px, py, pz] = point.into_arr();
                grids.push([px, py, pz]);
            }
        }
        self.shape = shape;
        self.grids = grids;
    }

    // 根据三个点获确定一个平面的格点
    /// 根据三个点确定一个平面并生成矩形格点
    pub fn set_v2(&mut self, p0: &Vec3, p1: &Vec3, p2: &Vec3, step: [f64; 2]) {
        let v01 = p1 - p0;
        let v02 = p2 - p0;
        let vx = v01.normalize();
        let vn = v02.normalize();
        let angle = v01.angle(&v02);
        let xp = v02.project_onto(&vx); // v02在vx上的投影
        let p00: Vec3;
        let p01: Vec3;
        if angle < PI / 2.0 {
            if xp.length() < v01.length() {
                p00 = p0.clone();
                p01 = p1.clone();
            } else {
                p00 = p0.clone();
                p01 = p0 + xp;
            }
        } else {
            p00 = p0 + xp;
            p01 = p1.clone()
        }
        let vz = vx.cross(&vn).normalize();
        let vy = vz.cross(&vx).normalize();
        let yp = v02.project_onto(&vy);
        let p11 = p0 + yp;
        let lx = (&p01 - &p00).length();
        let ly = (&p11 - &p00).length();
        let nx = (lx / step[0]).ceil() as usize;
        let ny = (ly / step[1]).ceil() as usize;
        self.set_v3(&p00, &p01, &p11, [nx, ny]);
    }

    /// 根据三个角点和尺寸生成矩形格点
    pub fn set_v3(&mut self, p00: &Vec3, p01: &Vec3, p11: &Vec3, size: [usize; 2]) {
        let vx = p01 - p00;
        let vy = p11 - p00;
        let [nx, ny] = size;
        let dx = vx / nx as f64;
        let dy = vy / ny as f64;
        let mut grids = vec![];
        for i in 0..nx {
            for j in 0..ny {
                let pos = p00 + &dx * i as f64 + &dy * j as f64;
                let [px, py, pz] = pos.into_arr();
                grids.push([px, py, pz]);
            }
        }
        let shape = [nx, ny];
        self.shape = shape;
        self.grids = grids;
    }

    /// 获取格点形状和数据
    pub fn get(&self) -> ([usize; 2], &Vec<[f64; 3]>) {
        (self.shape, &self.grids)
    }

    /// 根据三个角点计算平面方程参数 [a, b, c, d]，满足 ax + by + cz + d = 0
    pub fn paras(&self, p00: &Vec3, p01: &Vec3, p11: &Vec3) -> [f64; 4] {
        let vx = (p01 - p00).normalize();
        let vy = (p11 - p00).normalize();
        let vz = vx.cross(&vy);
        let [a, b, c] = vz.into_arr();
        let d = -vz.dot(&p00);
        [a, b, c, d]
    }
}

/// 三维立方体格点
pub struct CubeGrid {
    pub shape: [usize; 3],
    pub start: [f64; 3],
    pub steps: [f64; 3],
}

impl CubeGrid {
    /// 创建一个空的 CubeGrid 实例
    pub fn new() -> Self {
        let shape = [0, 0, 0];
        let start = [0.0, 0.0, 0.0];
        let steps = [0.0, 0.0, 0.0];
        CubeGrid { shape, start, steps }
    }

    /// 根据两个对角点和步长设置立方体格点参数
    pub fn set_v1(&mut self, p0: &[f64; 3], p1: &[f64; 3], step: f64, bord: f64) -> &Self {
        let p0 = p0.clone().map(|v| v - bord);
        let p1 = p1.clone().map(|v| v + bord);
        let [x0, y0, z0] = p0.into();
        let [x1, y1, z1] = p1.into();

        let nx = ((x1 - x0) / step).ceil() as usize;
        let ny = ((y1 - y0) / step).ceil() as usize;
        let nz = ((z1 - z0) / step).ceil() as usize;
        let shape = [nx as usize, ny as usize, nz as usize];
        self.shape = shape;
        self.start = p0.clone();
        self.steps = [step, step, step];
        self
    }

    /// 获取格点形状和所有格点坐标
    pub fn get(&self) -> ([usize; 3], Vec<[f64; 3]>) {
        let [nx, ny, nz] = self.shape;
        let [x0, y0, z0] = self.start;
        let [sx, sy, sz] = self.steps;
        let mut grids = Vec::with_capacity(nx * ny * nz);
        for i in 0..nx {
            for j in 0..ny {
                for k in 0..nz {
                    grids.push([x0 + sx * i as f64, y0 + sy * j as f64, z0 + sz * k as f64]);
                }
            }
        }
        (self.shape, grids)
    }
}

#[derive(Clone, Copy)]
/// 格点上的标量、梯度和 Hessian 矩阵值
pub struct GridVals {
    pub val0: f64,
    pub val1: [f64; 3],
    pub val2: [[f64; 3]; 3],
}

impl GridVals {
    /// 创建一个新的 GridVals 实例
    pub fn new(val0: f64, val1: [f64; 3], val2: [[f64; 3]; 3]) -> Self {
        GridVals { val0, val1, val2 }
    }

    pub fn set(&mut self, val0: f64, val1: [f64; 3], val2: [[f64; 3]; 3]) {
        self.val0 = val0;
        self.val1 = val1;
        self.val2 = val2;
    }

    /// 创建一个全零的 GridVals 实例
    pub fn non() -> Self {
        GridVals::new(0.0, [0.0; 3], [[0.0; 3]; 3])
    }

    /// 计算梯度向量的模长
    pub fn grad_norm(&self) -> f64 {
        let [gx, gy, gz] = self.val2;
        let norm = (gx[0].powi(2) + gy[0].powi(2) + gz[0].powi(2)).sqrt();
        norm
    }

    /// 计算梯度向量模长的平方
    pub fn grad_norm2(&self) -> f64 {
        let [gx, gy, gz] = self.val2;
        let norm2 = gx[0].powi(2) + gy[0].powi(2) + gz[0].powi(2);
        norm2
    }

    pub fn val2s(&self) -> [f64; 9] {
        let [gx, gy, gz] = self.val2;
        [gx[0], gx[1], gx[2], gy[0], gy[1], gy[2], gz[0], gz[1], gz[2]]
    }
}
