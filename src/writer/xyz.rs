use crate::datas::consts::BOHR;
use crate::datas::temps;
use crate::writer::Writer;
use std;

/// Xyz格式文件(.xyz)写入器
pub struct XyzWriter {
    /// 原子符号列表
    pub syms: Vec<String>,
    /// 原子坐标（bohr）
    pub xyzs: Vec<[f64; 3]>,
    /// 标题（注释行）
    pub title: String,
}

impl XyzWriter {
    /// 创建XyzWriter实例
    pub fn new() -> Self {
        Self {
            syms: Vec::new(),
            xyzs: Vec::new(),
            title: "generate by rswfn".to_string(),
        }
    }

    /// 生成xyz文件内容字符串
    pub fn build(&self) -> String {
        let mut geom = String::new();
        for i in 0..self.syms.len() {
            let sym = &self.syms[i];
            let [x, y, z] = &self.xyzs[i];
            geom += &format!("{:>3}           {:>14.8}{:>14.8}{:>14.8}\n", sym, x / BOHR, y / BOHR, z / BOHR);
        }
        let natm = self.syms.len();
        temps::XYZ.replace("{natm}", &format!("{}", natm)).replace("{title}", &self.title).replace("{geome}", &geom).to_string()
    }

    /// 保存到文件
    pub fn save(&self, path: &str) {
        let text = self.build();
        std::fs::write(path, text).expect("Unable to write file");
    }
}

impl Writer for XyzWriter {
    fn read_mole(&mut self, mole: &crate::base::Mole) {
        let syms = mole.atoms.syms();
        let xyzs = mole.atoms.xyzs();
        self.syms = syms;
        self.xyzs = xyzs.clone();
    }
}
