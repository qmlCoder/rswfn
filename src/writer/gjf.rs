use crate::datas::consts::BOHR;
use crate::{datas::temps, writer::Writer};

/// Gjf格式文件(.gjf)写入器
pub struct GjfWriter {
    /// 检查点文件名（chk）
    pub chk: String, // 检查点文件名
    /// 作业类型
    pub job: String, // 作业类型
    /// 电荷
    pub charge: isize, // 电荷和自旋
    /// 自旋多重度
    pub multip: usize, // 自旋多重度
    /// 原子符号列表
    pub syms: Vec<String>,
    /// 原子坐标（bohr）
    pub xyzs: Vec<[f64; 3]>,
}

impl GjfWriter {
    /// 创建GjfWriter实例
    pub fn new() -> Self {
        Self {
            chk: "".to_string(),
            job: "".to_string(),
            charge: 0,
            multip: 1,
            syms: Vec::new(),
            xyzs: Vec::new(),
        }
    }

    /// 生成gjf文件内容字符串
    pub fn build(&self) -> String {
        let mut geom = String::new();
        for i in 0..self.syms.len() {
            let sym = &self.syms[i];
            let [x, y, z] = &self.xyzs[i];
            geom += &format!("{:>3}           {:>14.8}{:>14.8}{:>14.8}\n", sym, x / BOHR, y / BOHR, z / BOHR);
        }
        temps::GJF
            .to_string()
            .replace("{chk}", &self.chk)
            .replace("{job}", &self.job)
            .replace("{ele.0}", format!("{:>2}", self.charge).as_str())
            .replace("{ele.1}", format!("{:>2}", self.multip).as_str())
            .replace("{geome}", &geom)
            .to_string()
    }

    /// 保存到文件
    pub fn save(&self, path: &str) {
        let text = self.build();
        std::fs::write(path, text).expect("Unable to write file");
    }
}

impl Writer for GjfWriter {
    fn read_mole(&mut self, mole: &crate::base::Mole) {
        let syms = mole.atoms.syms();
        let xyzs = mole.atoms.xyzs();
        self.charge = mole.charge();
        self.multip = mole.multip();
        self.syms = syms;
        self.xyzs = xyzs.clone();
    }
}
