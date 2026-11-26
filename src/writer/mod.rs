mod cub;
mod dmo;
mod fch;
mod gjf;
mod mfn;
mod obj;
mod xyz;

pub use cub::CubWriter;
pub use dmo::DmoWriter;
pub use fch::FchWriter;
pub use gjf::GjfWriter;
pub use mfn::MfnWriter;
pub use obj::ObjWriter;
pub use xyz::XyzWriter;

/// 文件写入器trait
///
/// 定义将分子信息写入各种文件格式的统一接口。
pub trait Writer {
    /// 从分子结构读取数据
    fn read_mole(&mut self, mole: &crate::base::Mole);
}
