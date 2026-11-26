use core::fmt;
use ndarray::Array2;
use std::collections::HashMap;
use std::f64::consts::PI;
use std::sync::{Mutex, OnceLock};

use crate::base::{Mole, Vec3};
use crate::datas::consts::BOHR;
use crate::datas::elems;
use crate::{datas, maths, utils};

// 原子局部坐标系，局部坐标系的每一列是一个基向量，从全局坐标系到局部坐标系的变换方式为 R_l = LT @ R_g
#[derive(Debug, Clone)]
pub struct Stm {
    pub ex: [f64; 3],
    pub ey: [f64; 3],
    pub ez: [f64; 3],
    pub mat: Array2<f64>, // 局部坐标系的矩阵，每一列是一个单位基向量
    pub non: bool,
}

impl Stm {
    /// 检查局部坐标系是否为正交矩阵（即旋转矩阵）
    ///
    /// # Returns
    /// - `true` 如果该矩阵是正交矩阵（逆矩阵等于转置矩阵）
    pub fn check(&self) -> bool {
        let mat_t = self.mat.t().to_owned();
        let mat_i = maths::inv(&self.mat);
        let same = utils::same_matrix(&mat_i, &mat_t);
        same
    }

    /// 从三个基向量创建局部坐标系
    ///
    /// # Parameters
    /// - `ex`: x方向的单位基向量，长度为3的数组
    /// - `ey`: y方向的单位基向量，长度为3的数组
    /// - `ez`: z方向的单位基向量，长度为3的数组
    pub fn new(ex: [f64; 3], ey: [f64; 3], ez: [f64; 3]) -> Self {
        let [xx, xy, xz] = ex;
        let [yx, yy, yz] = ey;
        let [zx, zy, zz] = ez;
        let mat = Array2::from_shape_vec((3, 3), vec![xx, yx, zx, xy, yy, zy, xz, yz, zz]).unwrap();
        let non = false;
        Self { ex, ey, ez, mat, non }
    }

    /// 创建一个空的局部坐标系（单位矩阵），标记为未设置
    pub fn blank() -> Self {
        let mut stm = Stm::new([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
        stm.non = true;
        stm
    }

    /// 返回局部坐标系的所有9个分量值
    ///
    /// 顺序为：ex_x, ex_y, ex_z, ey_x, ey_y, ey_z, ez_x, ez_y, ez_z
    pub fn vals(&self) -> Vec<f64> {
        let [xx, xy, xz] = self.ex;
        let [yx, yy, yz] = self.ey;
        let [zx, zy, zz] = self.ez;
        vec![xx, xy, xz, yx, yy, yz, zx, zy, zz]
    }

    /// 从3x3矩阵创建局部坐标系
    ///
    /// 矩阵的每一列对应一个基向量
    ///
    /// # Parameters
    /// - `mat`: 3x3矩阵，每列为单位基向量
    pub fn from_mat(mat: Array2<f64>) -> Self {
        let mut es = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                es[i][j] = mat[(j, i)];
            }
        }
        let ex = es[0];
        let ey = es[1];
        let ez = es[2];
        Self {
            ex,
            ey,
            ez,
            mat,
            non: false,
        }
    }

    /// 局部坐标系B在局部坐标系A中的朝向
    pub fn relative(&self, stmb: &Stm) -> Stm {
        let mata = self.mat.clone();
        let matb = stmb.mat.clone();
        let dmat = mata.t().dot(&matb);
        Stm::from_mat(dmat)
    }

    /// 返回局部坐标系的三个基向量
    pub fn dirs(&self) -> [[f64; 3]; 3] {
        [self.ex, self.ey, self.ez]
    }

    /// 从另一个局部坐标系复制所有数据到当前坐标系
    ///
    /// # Parameters
    /// - `stm`: 源局部坐标系
    pub fn copy(&mut self, stm: &Stm) {
        self.ex = stm.ex;
        self.ey = stm.ey;
        self.ez = stm.ex;
        self.mat = stm.mat.clone();
        self.non = stm.non;
    }

    /// 生成一个与单位向量 `axis` 正交的随机单位向量
    ///
    /// 若 10 次随机尝试均与 `axis` 几乎平行（概率极小），则退化为一个确定的垂直方向。
    fn rand_perp(axis: &[f64; 3]) -> [f64; 3] {
        // 随机方向去掉 axis 方向的分量，归一化后即为与 axis 正交的单位向量
        for _ in 0..10 {
            let rand = Vec3::rand().into_arr();
            let dot = rand[0] * axis[0] + rand[1] * axis[1] + rand[2] * axis[2];
            let orth = [
                rand[0] - dot * axis[0],
                rand[1] - dot * axis[1],
                rand[2] - dot * axis[2],
            ];
            if maths::vec_len(&orth) > 1e-6 {
                return maths::vec_norm(&orth);
            }
        }
        // 随机方向与 axis 几乎平行时，退化为一个确定的垂直方向
        let other = if axis[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        maths::vec_norm(&maths::vec_cross(axis, &other))
    }

    /// ex方向是指定的，另外两个方向是任意的
    ///
    /// `ey` 由随机方向对 `ex` 正交化得到，`ez = ex × ey`，因此返回的局部坐标系
    /// 始终是一组右手正交基。
    pub fn rand_ex(ex: [f64; 3]) -> Self {
        assert!(maths::vec_len(&ex) > 1e-12, "ex 不能是零向量");
        let ex = maths::vec_norm(&ex);
        let ey = Self::rand_perp(&ex);
        let ez = maths::vec_cross(&ex, &ey);
        Stm::new(ex, ey, ez)
    }

    /// ez方向是指定的，另外两个方向是任意的
    ///
    /// `ex` 由随机方向对 `ez` 正交化得到，`ey = ez × ex`，因此返回的局部坐标系
    /// 始终是一组右手正交基。
    pub fn rand_ez(ez: [f64; 3]) -> Self {
        assert!(maths::vec_len(&ez) > 1e-12, "ez 不能是零向量");
        let ez = maths::vec_norm(&ez);
        let ex = Self::rand_perp(&ez);
        let ey = maths::vec_cross(&ez, &ex);
        Stm::new(ex, ey, ez)
    }
}

impl fmt::Display for Stm {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let [xx, xy, xz] = self.ex;
        let [yx, yy, yz] = self.ey;
        let [zx, zy, zz] = self.ez;
        write!(f, "{}  ", if self.non { "F" } else { "T" })?;
        write!(f, "ex: ({:>8.4}{:>8.4}{:>8.4})  |  ", xx, xy, xz)?;
        write!(f, "ey: ({:>8.4}{:>8.4}{:>8.4})  |  ", yx, yy, yz)?;
        write!(f, "ez: ({:>8.4}{:>8.4}{:>8.4})", zx, zy, zz)?;
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub struct Atom {
    pub idx: usize,    // 原子在分子中的索引
    pub nuc: usize,    // 原子核电荷数
    pub crg: f64,      // 有效核电荷数，赝势基组需要
    pub xyz: [f64; 3], // 原子坐标
    pub rad: f64,      // 原子半径
    pub sym: String,   // 原子符号
    pub stm: Option<Stm>,
}

impl Atom {
    /// 创建一个新的原子
    ///
    /// # Parameters
    /// - `nuc`: 原子核电荷数
    /// - `xyz`: 原子坐标 [x, y, z]
    /// - `idx`: 原子在分子中的索引
    /// - `crg`: 有效核电荷数（赝势基组需要）
    pub fn new(nuc: usize, xyz: [f64; 3], idx: usize, crg: f64) -> Self {
        let mut rad = 0.0;
        let mut sym = None;
        for elem in datas::elems::ELEMS.iter() {
            if nuc as isize == elem.idx {
                rad = elems::atm2rad(nuc);
                sym = Some(elem.sym.clone());
                break;
            }
        }
        Self {
            nuc: nuc,
            crg: crg,
            xyz: xyz,
            idx: idx,
            rad: rad,
            sym: sym.unwrap(),
            stm: None,
        }
    }

    /// 将原子坐标转换为 Vec3 向量
    pub fn vec(&self) -> Vec3 {
        let [x, y, z] = self.xyz;
        Vec3::new(x, y, z)
    }

    /// 返回该原子在元素周期表中的行号
    pub fn irow(&self) -> usize {
        elems::atm2row(self.nuc)
    }
}

#[derive(Debug)]
pub struct Atoms {
    pub atoms: Vec<Atom>,
    pub dists: Array2<f64>,                  // 距离矩阵
    pub neles: [usize; 2],                   // α电子数和β电子数
    pub bonds: Mutex<Vec<[usize; 2]>>,       // 键列表
    nebs: Mutex<HashMap<usize, Vec<usize>>>, // 每个原子的相邻原子
    angs: Mutex<HashMap<[usize; 3], f64>>,   // 键角
    vecs: Mutex<HashMap<[usize; 2], Vec3>>,
    is_lines: Mutex<HashMap<usize, bool>>,
    is_bonds: Mutex<HashMap<[usize; 2], bool>>,
    nucs: OnceLock<Vec<usize>>,
    crgs: OnceLock<Vec<f64>>,
    idxs: OnceLock<Vec<usize>>,
    xyzs: OnceLock<Vec<[f64; 3]>>,
    stms: OnceLock<Vec<Option<Stm>>>,
}

impl Atoms {
    /// 从原子序数列表、坐标列表和有效核电荷数列表创建 Atoms 结构
    ///
    /// 自动计算距离矩阵、键矩阵以及α/β电子数
    ///
    /// # Parameters
    /// - `atms`: 原子序数列表
    /// - `xyzs`: 原子坐标列表
    /// - `crgs`: 有效核电荷数列表
    pub fn new(atms: &Vec<usize>, xyzs: &Vec<[f64; 3]>, crgs: &Vec<f64>) -> Self {
        let natm = atms.len();
        let mut atoms = Vec::with_capacity(natm);
        let mut nele = 0;
        for i in 0..natm {
            let nuc = atms[i];
            let xyz = xyzs[i];
            let idx = i;
            let crg = crgs[i];
            let atom = Atom::new(nuc, xyz, idx, crg);
            atoms.push(atom);
            nele += nuc;
        }
        let nela = nele / 2 + nele % 2;
        let nelb = nele / 2;
        let neles = [nela, nelb];
        let mut dists = Array2::zeros((natm, natm));
        let mut bonds = vec![];
        for i in 0..natm.saturating_sub(1) {
            for j in i + 1..natm {
                let dist = maths::calc_dist(&xyzs[i], &xyzs[j]);
                dists[(i, j)] = dist;
                dists[(j, i)] = dist;
                if dist < (atoms[i].rad + atoms[j].rad) * 1.1 {
                    bonds.push([i, j]);
                }
            }
        }
        Atoms {
            atoms: atoms,
            dists: dists,
            neles: neles,
            bonds: Mutex::new(bonds),
            nebs: Mutex::new(HashMap::new()),
            angs: Mutex::new(HashMap::new()),
            vecs: Mutex::new(HashMap::new()),
            is_lines: Mutex::new(HashMap::new()),
            is_bonds: Mutex::new(HashMap::new()),
            nucs: OnceLock::new(),
            crgs: OnceLock::new(),
            idxs: OnceLock::new(),
            xyzs: OnceLock::new(),
            stms: OnceLock::new(),
        }
    }

    /// 添加一个原子到分子中
    ///
    /// # Parameters
    /// - `atm`: 原子序数
    /// - `xyz`: 原子坐标 [x, y, z]
    /// - `crg`: 有效核电荷数
    pub fn add(&mut self, atm: usize, xyz: [f64; 3], crg: f64) {
        let idx = self.atoms.len() as usize;
        self.atoms.push(Atom::new(atm, xyz, idx, crg));
    }

    /// 获取指定索引的原子，如果索引越界则 panic
    ///
    /// # Parameters
    /// - `idx`: 原子索引
    pub fn get(&self, idx: usize) -> &Atom {
        &self.atoms.get(idx as usize).unwrap_or_else(|| {
            log::error!("原子索引超出范围:{}/{}", idx, self.atoms.len());
            panic!("原子索引超出范围")
        })
    }

    /// 返回原子总数
    pub fn len(&self) -> usize {
        return self.atoms.len();
    }

    /// 计算两个原子之间的距离
    pub fn dist(&self, i: usize, j: usize) -> f64 {
        let [xi, yi, zi] = self.atoms[i as usize].xyz;
        let [xj, yj, zj] = self.atoms[j as usize].xyz;
        ((xi - xj).powi(2) + (yi - yj).powi(2) + (zi - zj).powi(2)).sqrt()
    }

    /// 返回距离矩阵
    ///
    /// 矩阵元素 (i, j) 表示原子 i 和原子 j 之间的距离
    pub fn dmat(&self) -> Array2<f64> {
        let natm = self.len();
        let mut dmat = Array2::zeros((natm, natm));
        for i in 0..natm.saturating_sub(1) {
            for j in i + 1..natm {
                dmat[(i, j)] = self.dist(i, j);
                dmat[(j, i)] = dmat[(i, j)];
            }
        }
        dmat
    }

    /// 返回所有原子的坐标列表（惰性初始化）
    pub fn xyzs(&self) -> &Vec<[f64; 3]> {
        self.xyzs
            .get_or_init(|| self.atoms.iter().map(|atom| atom.xyz).collect())
    }

    /// 返回所有原子的坐标列表，单位转换为埃（Å）
    pub fn xyzs_ang(&self) -> Vec<[f64; 3]> {
        self.atoms
            .iter()
            .map(|atom| {
                let [x, y, z] = atom.xyz;
                [x / BOHR, y / BOHR, z / BOHR]
            })
            .collect()
    }

    /// 返回所有原子的半径列表
    pub fn rads(&self) -> Vec<f64> {
        let mut rads = Vec::new();
        for atom in self.atoms.iter() {
            rads.push(atom.rad);
        }
        rads
    }

    /// 返回所有原子的元素符号列表
    pub fn syms(&self) -> Vec<String> {
        let mut syms = Vec::new();
        for atom in self.atoms.iter() {
            syms.push(atom.sym.clone());
        }
        syms
    }

    /// 每个原子的核电荷数
    pub fn nucs(&self) -> &Vec<usize> {
        let nucs = self
            .nucs
            .get_or_init(|| self.atoms.iter().map(|atom| atom.nuc).collect::<Vec<usize>>());
        nucs
    }

    /// 每个原子的有效核电荷数
    pub fn crgs(&self) -> &Vec<f64> {
        let crgs = self
            .crgs
            .get_or_init(|| self.atoms.iter().map(|atom| atom.crg).collect::<Vec<f64>>());
        crgs
    }

    /// 返回所有原子的局部坐标系列表（惰性初始化）
    pub fn stms(&self) -> &Vec<Option<Stm>> {
        self.stms
            .get_or_init(|| self.atoms.iter().map(|atom| atom.stm.clone()).collect())
    }

    /// 原子索引范围 (0..natm-1)
    pub fn idxs(&self) -> &Vec<usize> {
        let idxs = self
            .idxs
            .get_or_init(|| self.atoms.iter().map(|atom| atom.idx).collect::<Vec<usize>>());
        idxs
    }

    /// 获取指定原子的邻接原子列表
    ///
    /// # Parameters
    /// - `atm`: 原子索引
    pub fn nebs(&self, atm: usize) -> Vec<usize> {
        let mut nebs = self.nebs.lock().expect("锁定失败");
        let nebs = nebs
            .entry(atm)
            .or_insert_with(|| {
                let mut nebs = Vec::new();
                for i in 0..self.len() {
                    if i == atm as usize {
                        continue;
                    }
                    if self.is_bond(atm, i) {
                        nebs.push(i);
                    }
                }
                nebs
            })
            .clone();
        // log::info!("获取邻接原子{}:{:?}", atm, nebs);
        nebs
    }

    /// ai->aj 的向量
    pub fn vec(&self, ai: usize, aj: usize) -> Vec3 {
        let mut vects = self.vecs.lock().expect("锁定失败");
        vects
            .entry([ai, aj])
            .or_insert_with(|| {
                let [xi, yi, zi] = self.get(ai).xyz;
                let [xj, yj, zj] = self.get(aj).xyz;
                Vec3::new(xj - xi, yj - yi, zj - zi)
            })
            .clone()
    }

    /// 计算三个原子之间的键角
    ///
    /// 计算 ai-aj-ak 三个原子之间的键角，以 aj 为顶点
    ///
    /// # Parameters
    /// - `ai`: 第一个原子索引
    /// - `aj`: 第二个原子索引（顶点）
    /// - `ak`: 第三个原子索引
    pub fn ang(&self, ai: usize, aj: usize, ak: usize) -> f64 {
        let mut angles = self.angs.lock().expect("锁定失败");
        *angles.entry([ai, aj, ak]).or_insert_with(|| {
            let veci = self.vec(aj, ai);
            let vecj = self.vec(aj, ak);
            // maths::vec_angle(&veci, &vecj)
            veci.angle(&vecj)
        })
    }

    /// 判断指定原子是否为线性结构
    ///
    /// # Parameters
    /// - `atm`: 原子索引
    pub fn is_line(&self, atm: usize) -> bool {
        let mut is_lines = self.is_lines.lock().expect("锁定失败");
        *is_lines.entry(atm).or_insert_with(|| {
            let nebs = self.nebs(atm);
            if nebs.len() == 1 {
                let neb = nebs[0];
                if self.nebs(neb).len() == 1 {
                    return true;
                } else if self.nebs(neb).len() == 2 {
                    let nebs = self.nebs(neb);
                    let nebs = nebs.iter().filter(|idx| *idx != &atm).copied().collect::<Vec<usize>>();
                    let idx = nebs[0];
                    let angle = self.ang(atm, neb, idx);
                    return (angle - PI).abs() < PI / 18.0;
                } else {
                    return false;
                }
            }
            if nebs.len() == 2 {
                let ai = nebs[0];
                let aj = nebs[1];
                let vi = self.vec(atm, ai);
                let vj = self.vec(atm, aj);
                let angle = vi.angle(&vj);
                return (angle - PI).abs() < PI / 18.0;
            }
            return false;
        })
    }

    /// 判断两个原子之间是否成键
    pub fn is_bond(&self, ai: usize, aj: usize) -> bool {
        let [ai, aj] = if ai < aj { [ai, aj] } else { [aj, ai] };
        let bonds = self.bonds.lock().expect("锁定失败");
        bonds.iter().any(|bond| bond[0] == ai && bond[1] == aj)
    }

    /// 返回成键列表（0 基原子序号对）
    ///
    /// 读取文件时按距离判据自动生成，`set_bonds` 可以整体覆盖。
    pub fn bonds(&self) -> Vec<[usize; 2]> {
        self.bonds.lock().expect("锁定失败").clone()
    }

    /// 手动设置原子间的成键关系
    ///
    /// 设置后会清除所有缓存的邻接原子和线性判断结果
    ///
    /// # Parameters
    /// - `bonds`: 成键关系映射，键为 [原子i, 原子j]，值为是否成键
    pub fn set_bonds(&self, bonds: Vec<[usize; 2]>) {
        log::debug!("设置成键: {bonds:?}");
        let mut bonds_store = self.bonds.lock().expect("锁定失败");
        *bonds_store = bonds;
        let mut is_bonds = self.is_bonds.lock().expect("锁定失败");
        is_bonds.clear();
        let mut is_lines = self.is_lines.lock().expect("锁定失败");
        is_lines.clear();
        let mut nebs = self.nebs.lock().expect("锁定失败");
        nebs.clear();
    }

    /// 获取指定原子列表的中心坐标
    ///
    /// # Parameters
    /// - `atms`: 原子索引列表
    pub fn cent(&self, atms: &Vec<usize>) -> Vec3 {
        let count = atms.len();
        let mut cent = [0.0, 0.0, 0.0];
        for atm in atms.iter() {
            let [px, py, pz] = self.get(*atm).xyz;
            cent[0] += px;
            cent[1] += py;
            cent[2] += pz;
        }
        cent[0] /= count as f64;
        cent[1] /= count as f64;
        cent[2] /= count as f64;
        Vec3::from_arr(cent)
    }

    /// 键长
    pub fn bond_len(&self, ai: usize, aj: usize) -> f64 {
        return self.dist(ai, aj);
    }

    /// 键角
    pub fn bond_ang(&self, ai: usize, aj: usize, ak: usize) -> f64 {
        let vecji = self.vec(aj, ai);
        let vecjk = self.vec(aj, ak);
        vecji.angle(&vecjk)
    }

    /// 二面角
    ///
    /// 计算 ai-aj-ak-al 四个连续原子的二面角（绕 aj-ak 键），返回弧度值
    /// 范围 (-π, π]，正值表示顺时针（右手定则沿 aj→ak 方向看）
    ///
    /// # Parameters
    /// - `ai`: 第一个原子索引
    /// - `aj`: 第二个原子索引
    /// - `ak`: 第三个原子索引
    /// - `al`: 第四个原子索引
    pub fn face_ang(&self, ai: usize, aj: usize, ak: usize, al: usize) -> f64 {
        // v1 = ri - rj, v2 = rk - rj（键轴）, v3 = rl - rk
        let v1 = self.vec(ai, aj).normalize();
        let v2 = self.vec(aj, ak).normalize();
        let v3 = self.vec(ak, al).normalize();

        // 平面 ai-aj-ak 的法向量 n1 = (ri-rj) × (rk-rj)
        let n1 = v1.cross(&v2);
        // 平面 aj-ak-al 的法向量 n2 = (rk-rj) × (rl-rk)
        let n2 = v2.cross(&v3);

        let n1_len = n1.length();
        let n2_len = n2.length();

        // 三点共线时无法定义二面角
        if n1_len < 1e-10 || n2_len < 1e-10 {
            return 0.0;
        }

        // cos(φ) = (n1·n2) / (|n1|·|n2|)
        let cos_phi = n1.dot(&n2) / (n1_len * n2_len);
        // 用 v2·(n1×n2) 确定符号（右手定则沿 v2 方向）
        let n1_cross_n2 = n1.cross(&n2);
        let sin_phi = v2.dot(&n1_cross_n2) / (n1_len * n2_len * v2.length());

        f64::atan2(sin_phi, cos_phi)
    }

    /// 判断两个原子结构是否相等
    ///
    /// 比较原子符号和坐标，坐标误差小于 1e-2 视为相等
    ///
    /// # Parameters
    /// - `other`: 另一个 Atoms 结构
    pub fn equal(&self, other: &Atoms) -> bool {
        let atoms_0 = &self.atoms;
        let atoms_1 = &other.atoms;
        if atoms_0.len() != atoms_1.len() {
            log::info!("原子数量不相等");
            return false;
        }
        for i in 0..atoms_0.len() {
            let atom_0 = atoms_0.get(i as usize).unwrap();
            let atom_1 = atoms_1.get(i as usize).unwrap();
            if atom_0.sym != atom_1.sym {
                log::info!("元素符号不相等");
                return false;
            }
            let [x0, y0, z0] = atom_0.xyz;
            let [x1, y1, z1] = atom_1.xyz;
            let dist = ((x1 - x0).powi(2) + (y1 - y0).powi(2) + (z1 - z0).powi(2)).sqrt();
            if dist > 1e-2 {
                log::info!("原子坐标不相等");
                return false;
            }
        }
        return true;
    }
}

impl Clone for Atoms {
    fn clone(&self) -> Self {
        let bonds = self.bonds.lock().unwrap().clone();
        Atoms {
            atoms: self.atoms.clone(),
            dists: self.dists.clone(),
            neles: self.neles.clone(),
            bonds: Mutex::new(bonds),
            angs: Mutex::new(HashMap::new()),
            nebs: Mutex::new(HashMap::new()),
            vecs: Mutex::new(HashMap::new()),
            is_lines: Mutex::new(HashMap::new()),
            is_bonds: Mutex::new(HashMap::new()),
            idxs: OnceLock::new(),
            crgs: OnceLock::new(),
            nucs: OnceLock::new(),
            xyzs: OnceLock::new(),
            stms: OnceLock::new(),
        }
    }
}

impl fmt::Display for Atoms {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let natm = self.atoms.len();
        writeln!(f, "{:>3}{:>14}{:>14}{:>14}{:>10}", "atm", "x", "y", "z", "crg")?;
        for i in 0..natm {
            let atom = self.get(i);
            let [x, y, z] = atom.xyz;
            writeln!(f, "{:>3}{:>14.8}{:>14.8}{:>14.8}{:>10.4}", atom.sym, x, y, z, atom.crg)?;
        }
        for i in 0..natm {
            let atom = self.get(i);
            if !atom.stm.is_none() {
                writeln!(f, "{}", atom.stm.as_ref().unwrap())?;
            }
        }
        Ok(())
    }
}

unsafe impl Send for Atoms {}
unsafe impl Sync for Atoms {}
