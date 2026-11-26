use crate::base::{Mole, Vec3};
use crate::maths;
use std::f64::consts::PI;

/// 键方向计算器，用于生成键轴周围的方向向量
pub struct Calculator<'m> {
    pub mole: &'m Mole,
}

impl<'m> Calculator<'m> {
    /// 创建一个新的 Calculator 实例
    pub fn new(mole: &'m Mole) -> Self {
        Calculator { mole }
    }

    /// 生成从原子 ai 到 aj 键轴周围的 36 个归一化方向向量
    pub fn verts(&self, ai: usize, aj: usize) -> Vec<Vec3> {
        let mut vectors = vec![];
        let vr = Vec3::rand();
        let bond_vec = self.mole.atoms().vec(ai, aj).normalize(); // 键轴单位向量
        let norm_vec = bond_vec.cross(&vr);
        for i in 0..36 {
            let vector = norm_vec.rotate_arrow(&bond_vec, i as f64 * PI / 18.0);
            vectors.push(vector.normalize());
        }
        vectors
    }
}
