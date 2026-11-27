use ndarray::Array2;

use crate::base::atoms::Atoms;
use crate::base::basis::Basis;
use crate::base::coefs::Coefs;
use crate::datas::elems;
use crate::reader::{FileType, Reader, get_fileType, read_lines};
use crate::{maths, reader};

use core::fmt;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::str::FromStr;

// use mem_dbg::MemDbg;

/// 分子结构，包含原子结构、基组和分子轨道系数
#[derive(Debug, Clone)]
pub struct Mole {
    /// 分子来源文件路径（若由文件读入则记录）
    pub path: Option<String>,
    /// 原子结构
    pub atoms: Atoms,
    /// 基组（可选）
    pub basis: Option<Basis>,
    /// 分子轨道系数（可选）
    pub coefs: Option<Coefs>,
}

impl Mole {
    /// 创建新的分子结构
    ///
    /// # Parameters
    /// - `atoms`: 原子结构
    /// - `basis`: 基组（可选）
    /// - `coefs`: 分子轨道系数（可选）
    pub fn new(atoms: Atoms, basis: Option<Basis>, coefs: Option<Coefs>) -> Mole {
        let mut mole = Self {
            path: None,
            atoms,
            basis,
            coefs,
        };
        mole
    }

    /// 从读取器读取分子信息
    ///
    /// # Parameters
    /// - `reader`: 实现了 Reader trait 的文件读取器
    pub fn from_read(reader: &dyn Reader) -> Mole {
        let atoms = reader.get_atoms();
        let basis = reader.get_basis();
        let coefs = reader.get_coefs();
        let mole = Mole::new(atoms, basis, coefs);
        mole
    }

    /// 从文件路径读取分子信息
    ///
    /// 自动根据文件扩展名选择合适的读取器
    /// 支持格式：gjf, log, fch, mol, sdf, xyz, mde, mwfn, dmo, inp, cub
    ///
    /// # Parameters
    /// - `path`: 文件路径
    pub fn from_file(path: &str) -> Mole {
        let file_type = get_fileType(path);
        let mut mole = match file_type {
            FileType::Gjf => Mole::from_read(&reader::GjfReader::new(path)),
            FileType::Log => Mole::from_read(&reader::LogReader::new(path)),
            FileType::Fch => Mole::from_read(&reader::FchReader::new(path)),
            FileType::Mol => Mole::from_read(&reader::MolReader::new(path)),
            FileType::Sdf => Mole::from_read(&reader::SdfReader::new(path)),
            FileType::Xyz => Mole::from_read(&reader::XyzReader::new(path)),
            FileType::Mde => Mole::from_read(&reader::MdeReader::new(path)),
            FileType::Mfn => Mole::from_read(&reader::MfnReader::new(path)),
            FileType::Dmo => Mole::from_read(&reader::DmoReader::new(path)),
            FileType::Inp => Mole::from_read(&reader::InpReader::new(path)),
            FileType::Cub => Mole::from_read(&reader::CubReader::new(path)),
        };
        mole.path = Some(path.to_string());
        mole
    }

    /// 返回分子中的原子结构
    pub fn atoms(&self) -> &Atoms {
        &self.atoms
    }

    /// 返回分子中的基组（如果存在），否则 panic
    pub fn basis(&self) -> &Basis {
        match &self.basis {
            Some(basis) => basis,
            None => panic!("该分子没有基组"),
        }
    }

    /// 返回分子中的轨道系数（如果存在），否则 panic
    pub fn coefs(&self) -> &Coefs {
        match &self.coefs {
            Some(coefs) => coefs,
            None => panic!("该分子没有系数"),
        }
    }

    /// 判断分子是否为开壳层（多重度不为 1）
    pub fn is_open(&self) -> bool {
        self.multip() != 1
    }

    /// 返回分子的化学式（如 "C2H6O"）
    pub fn formula(&self) -> String {
        let mut formula = String::new();
        let syms = self.atoms().syms();
        for elem in elems::ELEMS.iter() {
            let count = syms.iter().filter(|sym| **sym == elem.sym).count();
            if count > 0 {
                formula += &format!("{}{}", elem.sym, count);
            }
        }
        formula
    }

    /// 返回分子总核电荷数
    pub fn ncrg(&self) -> usize {
        let mut ncrg = 0.0;
        let crgs = self.atoms().crgs();
        for crg in crgs.iter() {
            ncrg += *crg;
        }
        ncrg as usize
    }

    /// 返回分子总电荷
    ///
    /// 电荷 = 总核电荷数 - 总电子数
    pub fn charge(&self) -> isize {
        let [nela, nelb] = self.atoms.neles;
        let ncrg = self.ncrg();
        ncrg as isize - nela as isize - nelb as isize
    }

    /// 返回分子多重度（自旋多重度）
    ///
    /// 多重度 = α电子数 - β电子数 + 1
    pub fn multip(&self) -> usize {
        let [nela, nelb] = self.atoms.neles;
        nela - nelb + 1
    }

    /// 计算格点到分子的最短距离
    ///
    /// # Parameters
    /// - `grid`: 三维格点坐标 [x, y, z]
    pub fn dist(&self, grid: &[f64; 3]) -> f64 {
        // 格点到分子的距离
        let [gx, gy, gz] = *grid;
        let atoms = self.atoms();
        let mut min_dist = f64::MAX;
        for i in 0..atoms.len() {
            let [ax, ay, az] = atoms.get(i as usize).xyz;
            let dist = ((gx - ax).powi(2) + (gy - ay).powi(2) + (gz - az).powi(2)).sqrt();
            if dist < min_dist {
                min_dist = dist;
            }
        }
        min_dist
    }

    /// 查找分子中的所有环（最小环集合，SSSR）
    ///
    /// 过滤掉环内存在多余键的非最小环，仅保留原子数等于键数的环
    pub fn rings(&self) -> Vec<Vec<usize>> {
        let atoms = self.atoms();
        let natm = atoms.len();
        let mut adjacency_list = HashMap::new();
        for i in 0..natm {
            for j in 0..natm {
                if atoms.is_bond(i as usize, j as usize) {
                    adjacency_list.entry(i as usize).or_insert(Vec::new()).push(j as usize);
                }
            }
        }
        let rings = find_sssr(&adjacency_list);
        // 如果环内有相邻的原子，则排除
        let rings = rings
            .iter()
            .filter(|ring| {
                let natm = ring.len();
                let mut nbond = 0;
                for i in 0..natm - 1 {
                    for j in i + 1..natm {
                        if atoms.is_bond(ring[i], ring[j]) {
                            nbond += 1;
                        }
                    }
                }
                nbond == natm // 原子的数量等于键的数量，说明是小环
            })
            .cloned()
            .collect();
        rings
    }

    /// 判断两个分子是否完全相等
    ///
    /// 比较原子结构、基组信息和轨道系数
    ///
    /// # Parameters
    /// - `other`: 另一个 Mole 结构
    pub fn equal(&self, other: &Mole) -> bool {
        if self.atoms().equal(&other.atoms()) == false {
            log::info!("分子结构不相等");
            return false;
        }
        if self.basis().equal(&other.basis()) == false {
            log::info!("分子基组不相等");
            return false;
        }
        if self.coefs().equal(&other.coefs()) == false {
            log::info!("分子系数不相等");
            return false;
        }
        return true;
    }

    /// 返回分子的空间边界
    ///
    /// # Returns
    /// 包含两个点的数组，[最小坐标, 最大坐标]，每个点为 [x, y, z]
    pub fn border(&self, bord: f64) -> [[f64; 3]; 2] {
        let xyzs = self.atoms().xyzs();
        let xs: Vec<f64> = xyzs.iter().map(|xyz| xyz[0]).collect();
        let ys: Vec<f64> = xyzs.iter().map(|xyz| xyz[1]).collect();
        let zs: Vec<f64> = xyzs.iter().map(|xyz| xyz[2]).collect();
        let minx = maths::arg_min(&xs).1;
        let miny = maths::arg_min(&ys).1;
        let minz = maths::arg_min(&zs).1;
        let maxx = maths::arg_max(&xs).1;
        let maxy = maths::arg_max(&ys).1;
        let maxz = maths::arg_max(&zs).1;
        let p0 = [minx - bord, miny - bord, minz - bord];
        let p1 = [maxx + bord, maxy + bord, maxz + bord];
        [p0, p1]
    }

    /// 判断分子中所有原子是否为线性结构
    pub fn is_line(&self) -> bool {
        let mut is_line = true;
        let natm = self.atoms.len();
        for i in 0..natm {
            if !self.atoms.is_line(i) {
                is_line = false;
            }
        }
        is_line
    }
}

impl fmt::Display for Mole {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let charge = self.charge();
        let multip = self.multip();
        writeln!(f, "{:>3}{:>3}", charge, multip)?;
        write!(f, "{}", self.atoms())?;
        if let Some(basis) = &self.basis {
            write!(f, "{}", basis)?;
        }
        if let Some(coefs) = &self.coefs {
            let atos = &self.basis().atos().raws;
            let coefs = format!("{}", self.coefs());
            let nato = atos.len();
            for (l, line) in coefs.split("\n").enumerate() {
                let idx = l % (nato + 2);
                let ato_str = if idx >= 2 {
                    let atm = atos[idx - 2].atm;
                    let elm = &self.atoms.get(atm).sym;
                    let shl = atos[idx - 2].shl;
                    let sym = &atos[idx - 2].sym;
                    let loc = atos[idx - 2].loc;
                    format!(
                        "{:>3}{:>3}{:>2}{:>3}{:>5}{:>3}",
                        idx - 1,
                        atm + 1,
                        elm,
                        shl + 1,
                        sym,
                        loc,
                    )
                } else {
                    format!("                   ")
                };
                writeln!(f, "{}{}", ato_str, line)?;
            }
        }
        Ok(())
    }
}

/// 查找图中的最小环集合（SSSR）
///
/// 使用深度优先搜索（DFS）遍历所有环，并去重
///
/// # Parameters
/// - `adjacency_list`: 邻接表表示的图
///
/// # Returns
/// 所有环的列表，每个环由原子索引列表表示
pub fn find_sssr(adjacency_list: &HashMap<usize, Vec<usize>>) -> Vec<Vec<usize>> {
    let n = adjacency_list.len();
    let mut rings: Vec<Vec<usize>> = Vec::new();

    fn dfs_find_rings(
        current: usize,
        start: usize,
        path: &mut Vec<usize>,
        visited: &mut HashSet<usize>,
        adjacency_list: &HashMap<usize, Vec<usize>>,
        rings: &mut Vec<Vec<usize>>,
    ) {
        path.push(current);
        visited.insert(current);

        if let Some(neighbors) = adjacency_list.get(&current) {
            for &neighbor in neighbors {
                if neighbor == start && path.len() >= 3 {
                    // 找到环
                    let ring = path.clone();
                    if is_new_ring(&ring, rings) {
                        rings.push(ring);
                    }
                } else if !visited.contains(&neighbor) {
                    dfs_find_rings(neighbor, start, path, visited, adjacency_list, rings);
                }
            }
        }

        path.pop();
        visited.remove(&current);
    }

    fn is_new_ring(candidate_ring: &[usize], existing_rings: &[Vec<usize>]) -> bool {
        let candidate_set: HashSet<usize> = candidate_ring.iter().cloned().collect();
        for ring in existing_rings {
            let ring_set: HashSet<usize> = ring.iter().cloned().collect();
            if candidate_set == ring_set {
                return false;
            }
        }
        true
    }

    for start_atom in 0..n as usize {
        let mut path = Vec::new();
        let mut visited = HashSet::new();
        dfs_find_rings(
            start_atom,
            start_atom,
            &mut path,
            &mut visited,
            adjacency_list,
            &mut rings,
        );
    }

    rings
}
