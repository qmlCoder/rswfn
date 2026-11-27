#![allow(non_snake_case)]
//! 原子局部方向与法向量计算。

use core::f64;
use std::sync::Mutex;
use std::{collections::HashMap, f64::consts::PI, vec};

use ndarray::Array2;

use crate::base::{Mole, Stm, Vec3};
use crate::datas;
use crate::maths;
use crate::props::gridprop;

/// 原子局部方向与法向量计算器。
///
/// 提供原子法向量、局部坐标系（LCS）以及反应方向生成等功能，
/// 用于后续的 π 轨道分析和方向自由价计算。
pub struct Calculator<'m> {
    /// 分子引用
    pub mole: &'m Mole,
    /// 构造局部坐标系时使用的随机辅助向量
    pub rand_vz: Vec3,
    normals: Mutex<HashMap<usize, Option<Vec3>>>, // 记录每个原子的法向量
}

impl<'m> Calculator<'m> {
    /// 创建一个新的方向计算器。
    ///
    /// 初始化随机辅助向量，若分子为线性则预计算所有原子的法向量。
    ///
    /// # Parameters
    /// - `mole` — 分子引用
    pub fn new(mole: &'m Mole) -> Self {
        // let rand_vz = maths::rand_vect();
        let rand_vz = Vec3::rand();
        let mut normals = HashMap::new();
        if mole.is_line() {
            let bond_vec = mole.atoms.vec(0, 1);
            let rand_vec = Vec3::rand();
            let norm_vec = bond_vec.cross(&rand_vec).normalize();
            for i in 0..mole.atoms.len() {
                normals.insert(i, Some(norm_vec.clone()));
            }
        }
        Calculator {
            mole: mole,
            rand_vz: rand_vz,
            normals: Mutex::new(normals),
        }
    }

    /// 获取原子的法向量。
    ///
    /// 优先从缓存获取，若无法直接计算则通过 BFS 搜索相邻原子的法向量。
    /// 对于只有一个相邻原子的情况，保证法向量垂直于键轴。
    ///
    /// # Parameters
    /// - `atm` — 原子索引
    ///
    /// # Returns
    /// 法向量，若无法确定则返回 `None`
    pub fn normal_vector(&self, atm: usize) -> Option<Vec3> {
        let atoms = self.mole.atoms();
        let nebs = atoms.nebs(atm);
        let mut vz = None;

        match self.get_normal(atm) {
            Some(normal) => vz = Some(normal),
            // 如果没法直接计算的话则迭代搜索
            None => {
                if nebs.len() == 4 {
                    vz = None;
                } else if nebs.len() == 5 {
                    vz = None;
                } else {
                    let mut searchd = vec![atm]; // 已经搜索过的
                    let mut searchs = atoms.nebs(atm).clone();
                    while searchs.len() > 0 {
                        let idx = searchs.pop().expect("弹出列表项失败");
                        searchd.push(idx);
                        let normal = self.get_normal(idx);
                        match normal {
                            Some(normal) => vz = Some(normal),
                            None => {
                                for neb in atoms.nebs(idx) {
                                    if searchd.contains(&neb) {
                                        continue;
                                    };
                                    if searchs.contains(&neb) {
                                        continue;
                                    };
                                    searchs.push(neb);
                                }
                            }
                        }
                    }
                }
            }
        };
        // 保证只有一个原子的时候法向量垂直于键轴
        vz = match vz {
            Some(vz) => {
                if nebs.len() == 1 {
                    let vx = atoms.vec(nebs[0], atm).normalize();
                    let vz = if vx.angle(&vz).abs() < 0.1 { Vec3::rand() } else { vz };
                    let vy = vz.cross(&vx).normalize();
                    Some(vx.cross(&vy).normalize())
                } else {
                    Some(vz)
                }
            }
            None => None,
        };
        // assert!()
        vz
    }

    /// 计算原子的法向量（含缓存）。
    ///
    /// 根据相邻原子数采用不同策略:
    /// - 1 个相邻原子: 无法确定，返回 `None`
    /// - 2 个相邻原子: 通过两键向量的叉积得到法向量
    /// - 3 个相邻原子: 通过三点确定的平面计算法向量
    /// - 6 个相邻原子: 碳硼烷特例，取与其他键最小排斥方向
    ///
    /// # Parameters
    /// - `atm` — 原子索引
    ///
    /// # Returns
    /// 法向量，若无法确定则返回 `None`
    pub fn get_normal(&self, atm: usize) -> Option<Vec3> {
        let mut normals = self.normals.lock().expect("锁定失败");
        normals
            .entry(atm)
            .or_insert_with(|| {
                let atoms = self.mole.atoms();
                let nebs = atoms.nebs(atm);
                if nebs.len() == 1 {
                    return None;
                }
                if nebs.len() == 2 {
                    let ai = nebs[0];
                    let aj = nebs[1];
                    let vi = atoms.vec(atm, ai).normalize();
                    let vj = atoms.vec(atm, aj).normalize();
                    let angle = vi.angle(&vj);
                    if (angle - PI).abs() < 0.01 {
                        return None;
                    } else {
                        let mut normal = vi.cross(&vj).normalize();
                        // let [dx, dy, dz] = self.rand_vz.clone().into();
                        // let angle = maths::vec_angle(&normal.into(), &[dx, dy, dz]);
                        let angle = normal.angle(&self.rand_vz);
                        if angle.abs() > 0.5 * PI {
                            normal = normal * -1.0;
                        }
                        return Some(normal);
                    }
                }
                if nebs.len() == 3 {
                    let pi: Vec3 = atoms.get(nebs[0]).vec();
                    let pj: Vec3 = atoms.get(nebs[1]).vec();
                    let pk: Vec3 = atoms.get(nebs[2]).vec();
                    let vij: Vec3 = (&pj - &pi).normalize();
                    let vik: Vec3 = (&pk - &pi).normalize();
                    let mut normal = vij.cross(&vik).normalize();
                    let [A, B, C] = normal.into_arr();
                    let [cx, cy, cz] = ((pi + pj + pk) / 3.0).into_arr();
                    let D = -A * cx - B * cy - C * cz;
                    let [ax, ay, az] = atoms.get(atm).xyz;
                    let dist = f64::abs(A * ax + B * ay + C * az + D) / f64::sqrt(A.powi(2) + B.powi(2) + C.powi(2)); // 点到平面的距离
                    if dist < 0.1 {
                        // 是一个平面
                        let rand_vect = self.rand_vz.clone();
                        // let angle = maths::vec_angle(&normal.into(), &rand_vect.into());
                        let angle = normal.angle(&rand_vect);
                        if angle > 0.5 * PI {
                            normal *= -1.0;
                        }
                        return Some(normal);
                    } else {
                        let vect = atoms.get(atm).vec() - Vec3::new(cx, cy, cz);
                        let angle = normal.angle(&vect);
                        if angle > 0.5 * PI {
                            normal *= -1.0;
                        }
                        return Some(normal);
                    }
                }
                if nebs.len() == 6 {
                    // 碳硼烷的特例，计算每一个键和其它原子的最小排斥方向
                    let mut min_val = f64::MIN;
                    let mut min_idx = nebs[0];
                    for i in 0..nebs.len() {
                        let mut val = f64::MAX;
                        for j in 0..nebs.len() {
                            if i == j {
                                continue;
                            }
                            let angle = atoms.ang(nebs[i], atm, nebs[j]);
                            if angle < val {
                                val = angle;
                            }
                        }
                        if val > min_val {
                            min_val = val;
                            min_idx = nebs[i];
                        }
                    }
                    let vect = atoms.vec(atm, min_idx);
                    return Some(vect.normalize());
                } else {
                    return Some(self.min_exclusion(atm));
                }
                return None;
            })
            .clone()
    }

    /// 计算法向量修正矩阵。
    ///
    /// 对于高配位数（≥5）的相邻原子对，根据其法向量夹角的余弦值
    /// 对键级矩阵进行修正。
    ///
    /// # Returns
    /// `natm × natm` 的修正矩阵
    pub fn dmat_fix(&self) -> Array2<f64> {
        let mut dirs = vec![];
        let atoms = self.mole.atoms();
        let natm = atoms.len();
        for atm in 0..natm {
            let normal = self.get_normal(atm as usize);
            dirs.push(normal);
        }
        let mut dmat_fix = Array2::<f64>::ones((natm, natm));
        for ai in 0..natm {
            for aj in 0..natm {
                if !atoms.is_bond(ai as usize, aj as usize) {
                    continue;
                }
                if atoms.nebs(ai as usize).len() < 5 && atoms.nebs(aj as usize).len() < 5 {
                    continue;
                }
                let c = match (dirs[ai].clone(), dirs[aj].clone()) {
                    (Some(di), Some(dj)) => {
                        let ang = di.angle(&dj);
                        ang.cos()
                    }
                    _ => 1.0,
                };
                dmat_fix[(ai, aj)] = c;
            }
        }
        dmat_fix
    }

    /// 计算原子所在的局部坐标系 local coord system
    pub fn LCS(&mut self, atm: usize, neb: Option<usize>) -> Option<Stm> {
        let atoms = self.mole.atoms();
        let atom = atoms.get(atm);
        if !atom.stm.is_none() {
            return atom.stm.clone();
        }
        let nebs = atoms.nebs(atm);

        let normal: Option<Vec3> = match self.normal_vector(atm) {
            Some(normal) => Some(normal),
            None => {
                let idx1;
                let idx2;
                if nebs.len() == 1 {
                    idx1 = nebs[0];
                    idx2 = atm;
                } else if nebs.len() == 2 {
                    idx1 = nebs[0];
                    idx2 = nebs[1];
                } else {
                    return None;
                }

                let v_bond = (atoms.get(idx1).vec() - atoms.get(idx2).vec()).normalize();
                let vr = Vec3::rand();
                Some(v_bond.cross(&vr).normalize())
            }
        };
        let stm = match normal {
            None => None,
            Some(normal) => {
                let v_norm = normal.normalize(); // 法向量
                let [v_norm, v_join, v_bond] = match neb {
                    Some(neb) => {
                        // 如果有相邻原子
                        let v_bond = (atoms.get(neb).vec() - atom.vec()).normalize();
                        let v_norm = if v_bond.angle(&v_norm).abs() < 0.1 { Vec3::rand() } else { v_norm };
                        let v_join = v_norm.cross(&v_bond).normalize();
                        let v_norm = v_bond.cross(&v_join).normalize(); // 重新正交化
                        [v_norm, v_join, v_bond]
                    }
                    None => {
                        if nebs.len() == 1 {
                            let v_bond: Vec3 = atoms.vec(atm, nebs[0]).normalize(); // 这个是确定的，键轴的方向
                            // 如果vx和vz在同一个方向，则重新整个vz
                            let v_norm = if v_bond.angle(&v_norm).abs() < 0.01 { Vec3::rand() } else { v_norm };
                            let v_join = v_norm.cross(&v_bond).normalize();
                            let v_norm = v_bond.cross(&v_join).normalize(); // 重新正交化
                            // log::info!("原子:{},v_bond:{},v_join:{},v_norm:{}", atm + 1, v_bond, v_join, v_norm);
                            [v_norm, v_join, v_bond]
                        } else if nebs.len() == 2 {
                            let idx1 = nebs[0];
                            let idx2 = nebs[1];
                            let v_bond = atoms.vec(idx1, idx2).normalize();
                            let v_join = v_norm.cross(&v_bond).normalize();
                            let v_bond = v_join.cross(&v_norm).normalize(); // 重新正交化
                            [v_norm, v_join, v_bond]
                        } else {
                            let cent = atoms.cent(&nebs); // 相邻原子的中心
                            let vect = cent - atom.vec(); // 中心到原子的向量
                            let dist = vect.length();
                            let vr = if dist > 0.01 {
                                vect.normalize() // 还不一定垂直于z轴
                            } else {
                                Vec3::rand()
                            };
                            let vr = if vr.angle(&v_norm) < 0.01 { Vec3::rand() } else { vr };
                            let v_join = v_norm.cross(&vr).normalize();
                            let v_bond = v_join.cross(&v_norm).normalize();
                            [v_norm, v_join, v_bond]
                        }
                    }
                };
                let ex: [f64; 3] = v_norm.into_arr();
                let ey: [f64; 3] = v_join.into_arr();
                let ez: [f64; 3] = v_bond.into_arr();
                let stm = Stm::new(ex, ey, ez);
                // log::info!("atm:{}", atm);
                if !stm.check() {
                    log::error!("局部坐标系有误");
                }
                Some(stm)
            }
        };
        // let mut stm = self.RLCS(atm);
        // stm.ey = [0.0, 1.0, 0.0];
        // stm.ez = [-1.0, 0.0, 0.0];
        // Some(stm)
        stm
    }

    /// 随机局部坐标系
    pub fn RLCS(&self, atm: usize) -> Stm {
        let vx = match self.normal_vector(atm) {
            Some(normal) => normal,
            None => {
                log::warn!("原子{}没有法向量", atm);
                Vec3::rand()
            }
        };
        let vy = Vec3::rand();
        let vz = vx.cross(&vy).normalize();
        let vy = vz.cross(&vx).normalize();
        let ex = vx.into_arr();
        let ey = vy.into_arr();
        let ez = vz.into_arr();
        Stm::new(ex, ey, ez)
    }

    /// 多重局部坐标系
    pub fn MLCS(&mut self, atm: usize) -> Vec<Stm> {
        let nebs = self.mole.atoms().nebs(atm);
        let mut mlcs = vec![];
        for idx in nebs.iter() {
            if let Some(stm) = self.LCS(*idx, Some(atm)) {
                mlcs.push(stm);
            }
        }
        mlcs
    }

    /// 生成原子周围的反应方向采样点。
    ///
    /// 根据相邻原子数量采用不同策略:
    /// - 1 个相邻原子: 在键反向半球均匀采样
    /// - 2 个相邻原子: 在键角平分面周围均匀采样
    /// - 3 个相邻原子: 基于平面法向量和键方向生成
    ///
    /// # Parameters
    /// - `atm` — 原子索引
    ///
    /// # Returns
    /// 单位方向向量列表
    pub fn reactions(&self, atm: usize) -> Vec<Vec3> {
        let atoms = self.mole.atoms();
        let nebs = atoms.nebs(atm);
        let mut vectors = Vec::new();
        if nebs.len() == 1 {
            let vr: Vec3 = Vec3::rand(); // 随机向量
            let bond_vect: Vec3 = atoms.vec(atm, nebs[0]).normalize(); // 键向量
            let vert_vect: Vec3 = bond_vect.cross(&vr); // 垂直于键的向量
            let atom_cent = atoms.get(nebs[0]).vec(); // 原子的中心

            let point = &atom_cent + &vert_vect;

            let axis_i = bond_vect.cross(&vert_vect);
            let axis_j = bond_vect.clone();

            for i in 0..4 {
                let angle = i as f64 * PI / 8.0;
                let diri = vert_vect.rotate_arrow(&axis_i, angle);
                for j in 0..36 {
                    let angle = j as f64 * PI / 18.0;
                    let dirj = diri.rotate_arrow(&axis_j, angle);
                    // let vec = pj - &atom_cent;
                    vectors.push(dirj.normalize());
                }
            }
            vectors.push(bond_vect.clone() * -1.0);
            return vectors;
            // return self.gt_bond_angle(atm);
        }
        if nebs.len() == 2 {
            let ai = nebs[0];
            let aj = nebs[1];
            let vi = atoms.vec(atm, ai).normalize();
            let vj = atoms.vec(atm, aj).normalize();
            let angle = vi.angle(&vj); // 两向量之间的夹角
            let cent = Vec3::zero();
            if atoms.is_line(atm) {
                // 线性分子
                let vr = Vec3::rand();
                let bond_vec = atoms.vec(ai, aj).normalize();
                let vert_vec = bond_vec.cross(&vr);

                for i in 0..36 {
                    // let vector = maths::rotate_point(&vert_vec, &cent, &bond_vec, i as f64 * PI / 18.0);
                    let vector = vert_vec.rotate_arrow(&bond_vec, i as f64 * PI / 18.0);
                    vectors.push(vector.normalize());
                }
            } else {
                let angle = PI - angle; // 有效角度
                let vert_vec = vi.cross(&vj).normalize();
                let mut count = (angle / (PI / 18.0)) as usize;
                if count < 3 {
                    count = 3;
                }
                let dang = angle / (count as f64);
                for i in 0..1 + count {
                    // 水平旋转
                    let rang = dang * i as f64; // 旋转角度，跳过90度
                    let axis = vi.rotate_arrow(&vert_vec, -rang).normalize();
                    for j in 0..18 {
                        // 垂直旋转
                        // 垂直方向都是18个
                        // let vector = maths::rotate_point(&vert_vec, &cent, &axis, i as f64 * PI / 18.0);
                        if i != 0 && (j == 0 || j == 17) {
                            continue;
                        }
                        let vector = vert_vec.rotate_arrow(&axis, j as f64 * PI / 18.0);
                        vectors.push(vector.normalize());
                    }
                }
            }
            return vectors;
        }
        if nebs.len() == 3 {
            let atom = atoms.get(atm);
            let ai = nebs[0];
            let aj = nebs[1];
            let ak = nebs[2];
            let pi = atoms.get(ai).vec();
            let pj = atoms.get(aj).vec();
            let pk = atoms.get(ak).vec();

            let vi = atoms.vec(atm, ai);
            let vj = atoms.vec(atm, aj);
            let vk = atoms.vec(atm, ak);

            let plan_params = maths::get_plan_by_3points(&pi, &pj, &pk);
            let distance = atom.vec().dist_to_plan(&plan_params);
            let is_plane = distance.abs() < 1e-2;
            let [a, b, c, _d] = plan_params;
            let mut normal = Vec3::new(a, b, c).normalize();
            let angle = normal.angle(&self.rand_vz);
            if angle.abs() > 0.5 * PI {
                normal = normal * -1.0;
            }
            if is_plane {
                // 如果是平面
                vectors.push(normal.clone());
                vectors.push(normal * -1.0);
            } else {
                let mean_vec = (vi + vj + vk) / 3.0;
                let angle = mean_vec.angle(&normal);
                if angle < PI / 2.0 {
                    normal = normal * -1.0;
                }
                vectors.push(normal.clone());
                for neb in nebs.iter() {
                    let neb_vec = atoms.vec(atm, *neb).normalize(); // vx
                    let vy = normal.cross(&neb_vec).normalize(); // vy=vz@vx
                    let vz = neb_vec.cross(&vy).normalize();
                    vectors.push(vz);
                }
            }
            return vectors;
        }
        return vectors;
    }

    /// 通过格点排斥法获取原子周围的可达方向。
    ///
    /// 在格点中筛选出与所有相邻原子键向量夹角均大于 90° 的方向。
    ///
    /// # Parameters
    /// - `atm` — 原子索引
    ///
    /// # Returns
    /// 满足排斥条件的方向向量列表
    pub fn gt_bond_angle(&self, atm: usize) -> Vec<Vec3> {
        // let dir_caler=atomprop::direction::Calculator::new(self.mole);
        let gridCaler = gridprop::dftgrid::Calculator::new(self.mole);
        let (grids, _weits) = gridCaler.atm_grids(atm, 1, 434);
        let atom = self.mole.atoms().get(atm);
        let nebs = self.mole.atoms().nebs(atm);
        let atm_xyz = atom.vec();
        let mut dirs = vec![];
        for g in 0..grids.len() {
            let grid = Vec3::from_arr(grids[g]);
            let mut meet = true; // 该点是否满足条件
            let vec_ag = (&grid - &atm_xyz).normalize();
            for neb in nebs.iter() {
                let neb_xyz = self.mole.atoms().get(*neb).vec();
                let vec_ab = (&neb_xyz - &atm_xyz).normalize();
                // let angle = maths::vec_angle(&vec_ab.into(), &vec_ag.into());
                let angle = vec_ab.angle(&vec_ag);
                if angle < 0.5 * PI {
                    meet = false;
                }
            }
            if meet {
                dirs.push(vec_ag);
            }
        }
        dirs
    }

    /// 为线性分子（键角约 180°）生成反应方向。
    ///
    /// 若原子处于线性环境则委托 `reactions` 生成方向；否则返回空列表。
    ///
    /// # Parameters
    /// - `atm` — 原子索引
    ///
    /// # Returns
    /// 单位方向向量列表
    pub fn linear_circle(&self, atm: usize) -> Vec<Vec3> {
        let dirs: Vec<Vec3> = Vec::new();
        let atoms = self.mole.atoms();
        let nebs = atoms.nebs(atm);
        if nebs.len() == 2 {
            let ai = nebs[0];
            let aj = nebs[1];
            let vi = atoms.vec(atm, ai);
            let vj = atoms.vec(atm, aj);
            let angle = vi.angle(&vj);
            if (angle - PI).abs() < 1e-3 {
                return self.reactions(atm);
            } else {
                return dirs;
            }
        } else {
            return dirs;
        }
    }

    /// 获取与所有键轴最小排斥力的方向
    pub fn min_exclusion(&self, atm: usize) -> Vec3 {
        fn exclusion(dir: &Vec3, bond_dirs: &[Vec3]) -> f64 {
            let angles: Vec<f64> = bond_dirs.iter().map(|bond_dir| dir.angle(bond_dir)).collect();
            // 计算角度的RMSE
            let rmse = angles.iter().fold(0.0, |acc, angle| acc + angle.powi(2)) / angles.len() as f64;
            rmse
        }
        // 方向初猜
        let init_dirs: Vec<Vec3> = datas::lebedev::ld0014().0.iter().map(|dir| Vec3::from_arr(*dir)).collect();
        let bond_dirs: Vec<Vec3> = {
            let atoms = self.mole.atoms();
            let nebs = atoms.nebs(atm);
            nebs.iter().map(|neb| atoms.vec(atm, *neb)).collect()
        };
        let mut idx = 0;
        let mut score = 0.0;
        for (i, dir) in init_dirs.iter().enumerate() {
            let val = exclusion(dir, &bond_dirs);
            if val > score {
                score = val;
                idx = i;
            }
        }
        // 获取初猜方向
        let mut init_dir = init_dirs[idx].clone();
        let mut step = 0.01;
        let z_axis = Vec3::new(0.0, 0.0, 1.0);
        let x_axis = Vec3::new(1.0, 0.0, 0.0);

        while step > 0.0001 {
            let next_dirs = vec![
                init_dir.clone().rotate_arrow(&z_axis, step),
                init_dir.clone().rotate_arrow(&z_axis, -step),
                init_dir.clone().rotate_arrow(&x_axis, step),
                init_dir.clone().rotate_arrow(&x_axis, -step),
            ];
            let vals: Vec<f64> = next_dirs.iter().map(|dir| exclusion(dir, &bond_dirs)).collect();
            let (max_idx, max_val) = maths::arg_max(&vals);
            if max_val > score {
                init_dir = next_dirs[max_idx].clone();
            } else {
                if step > 0.0001 {
                    step /= 10.0;
                } else {
                    break;
                }
            }
            score = max_val;
        }
        init_dir
    }
}
