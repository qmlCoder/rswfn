use std::fs::File;
use std::io::Write;

/// Obj格式文件(.obj)写入器,用于3D模型导出
pub struct ObjWriter {
    /// 顶点坐标列表
    pub verts: Vec<[f64; 3]>,
    /// 三角面片列表（每项为三个顶点的索引）
    pub faces: Vec<[usize; 3]>,
}

impl ObjWriter {
    /// 创建ObjWriter实例
    pub fn new(verts: Vec<[f64; 3]>, faces: Vec<[usize; 3]>) -> ObjWriter {
        ObjWriter { verts, faces }
    }

    /// 保存到文件
    pub fn save(&self, path: &str) {
        let mut file = File::create(path).expect("创建文件失败");
        for vert in self.verts.iter() {
            let [x, y, z] = *vert;
            file.write_all(format!("v {:>10.4}{:>10.4}{:>10.4}\n", x, y, z).as_bytes()).expect("写入失败");
        }
        for face in self.faces.iter() {
            let [a, b, c] = *face;
            file.write_all(format!("f {:>4} {:>4} {:>4}\n", a + 1, b + 1, c + 1).as_bytes()).expect("写入失败");
        }
        log::info!("文件写入成功");
    }
}

// impl Writer for ObjWriter {
//     fn from_mole(&mut self, mole: &crate::base::Mole) {}
// }
