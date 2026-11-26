use crate::base::Mole;
use crate::maths;
use crate::props::moleprop;
use crate::props::orbtprop;

use ndarray::Array2;

/// 布居分析计算器
pub struct Calculator<'m> {
    pub mole: &'m Mole,
}

impl<'m> Calculator<'m> {
    /// 创建一个新的 Calculator 实例
    pub fn new(mole: &'m Mole) -> Self {
        Calculator { mole }
    }
}
