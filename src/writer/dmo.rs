use crate::datas::consts::BOHR;
use crate::writer::Writer;
use ndarray::Array2;
use std::fs::File;
use std::io::prelude::*;
/// Dmo/Dmi格式文件写入器
pub struct DmoWriter {
    /// 原子符号列表
    pub syms: Vec<String>,
    /// 原子坐标（bohr）
    pub xyzs: Vec<[f64; 3]>,
    /// 原子轨道列表，每个元素为 `(原子索引, 壳层索引, 轨道符号, 轨道方向向量)`
    pub atos: Vec<(usize, usize, String, [f64; 3])>,
    /// 键连矩阵（原子轨道 × 原子轨道）
    pub BM: Array2<f64>,
    /// 轨道系数矩阵（原子轨道 × 分子轨道）
    pub CM: Array2<f64>,
    /// 各分子轨道的能量
    pub es: Vec<f64>,
}

impl DmoWriter {
    /// 创建一个内容为空的写入器
    pub fn blank() -> Self {
        Self {
            syms: vec![],
            xyzs: vec![],
            atos: vec![],
            BM: Array2::zeros((0, 0)),
            CM: Array2::zeros((0, 0)),
            es: vec![],
        }
    }
    /// 将数据保存到文件
    ///
    /// # Parameters
    /// - `path`: 目标文件路径
    pub fn save(&self, path: &str) {
        let mut file = File::create(path).expect("创建文件失败");
        let natm = self.syms.len();
        // 保存几何结构
        writeln!(&mut file, "#atms").unwrap();

        for atm in 0..natm {
            let sym = &self.syms[atm];
            let [x, y, z] = self.xyzs[atm];
            write!(
                &mut file,
                "{:>2}{:>14.8}{:>14.8}{:>14.8}\n",
                sym,
                x / BOHR,
                y / BOHR,
                z / BOHR
            )
            .unwrap();
        }
        let atos = &self.atos;
        let nato = atos.len();
        // 保存轨道朝向
        if !atos.is_empty() {
            writeln!(&mut file, "#atos").unwrap();

            for ato in atos.iter() {
                let [dx, dy, dz] = ato.3;
                let atm = ato.0;
                let sym = &ato.2;
                write!(
                    &mut file,
                    "{:>5}{:>3}{:>3}{:>14.8}{:>14.8}{:>14.8}\n",
                    atm + 1,
                    2,
                    sym,
                    dx,
                    dy,
                    dz
                )
                .unwrap();
            }
        }

        // 保存键连矩阵
        let BM = &self.BM;
        if !BM.is_empty() {
            writeln!(&mut file, "#BM").unwrap();
            for i in 0..nato {
                for j in i..nato {
                    let ai = atos[i].0;
                    let aj = atos[j].0;
                    let bv = BM[(i, j)];
                    if bv != 0.0 {
                        write!(
                            &mut file,
                            "{:>5}{:>5}{:>5}{:>5}{:>14.8}\n",
                            ai + 1,
                            aj + 1,
                            i + 1,
                            j + 1,
                            bv
                        )
                        .unwrap();
                    }
                }
            }
        }

        // 保存系数矩阵
        let CM = &self.CM;
        if !CM.is_empty() {
            writeln!(&mut file, "#CM").unwrap();
            let nato = CM.nrows();
            let nobt = CM.ncols();
            let engs = &self.es;
            for b in (0..nobt).step_by(10) {
                write!(&mut file, "{:>10}", "").unwrap();
                for obt in b..(b + 10).min(nobt) {
                    write!(&mut file, "{:>10}", obt + 1).unwrap();
                }
                write!(&mut file, "\n").unwrap();
                write!(&mut file, "{:>10}", "").unwrap();
                for obt in b..(b + 10).min(nobt) {
                    write!(&mut file, "{:>10.4}", engs[obt]).unwrap();
                }
                write!(&mut file, "\n").unwrap();
                for ato in 0..nato {
                    let atm = atos[ato].0;
                    let sym = &atos[ato].1;
                    write!(&mut file, "{:>5}{:<5}", atm + 1, sym).unwrap();
                    for obt in b..(b + 10).min(nobt) {
                        write!(&mut file, "{:>10.4}", CM[(ato, obt)]).unwrap();
                    }
                    write!(&mut file, "\n").unwrap();
                }
            }
        }
    }
}

impl Writer for DmoWriter {
    fn read_mole(&mut self, mole: &crate::base::Mole) {
        self.syms = mole.atoms().syms().clone();
        self.xyzs = mole.atoms().xyzs().clone();
    }
}
