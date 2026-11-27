//! 键能量的可加分解（闭壳层、杂化泛函）。

use crate::base::Mole;

use super::energy_v1;

/// 键能量的可加分解结果（闭壳层，杂化泛函）
///
/// 把完整电子能量按密度拆成"键"与"其余"两部分：只保留两个原子的系数得到
/// `D_bond`，其余部分取补密度 `D_rest = D − D_bond`，于是 `D = D_bond + D_rest`
/// 精确成立。由于 `J`、`K` 对 `D` 线性，库仑/精确交换可精确拆开；`E_xc` 非线性，
/// 用 `E_xc[ρ]` 与 `E_xc[ρ_bond]`、`E_xc[ρ_rest]` 的差值归入交叉项，从而保证
///
/// `E_bond + E_rest + E_cross = E_elec`（严格成立）。
#[derive(Debug, Clone, Copy)]
pub struct BondEnergy {
    /// 键自能量 `E_bond`：两原子（A∪B）密度块自身的动能/核吸引、库仑/交换与 XC
    pub e_bond: f64,
    /// 其余部分能量 `E_rest`
    pub e_rest: f64,
    /// 键-其余交叉项 `E_cross`（含 XC 交叉/残差）
    pub e_cross: f64,
    /// 交叉项中的库仑部分 `Tr(D_bond J_rest)`（键电子与其余电子的静电排斥，通常为正）
    pub e_cross_coul: f64,
    /// 交叉项中的精确交换部分 `−aₓ/2·Tr(D_bond K_rest)`（通常为负）
    pub e_cross_exch: f64,
    /// 交叉项中的 XC 残差 `E_xc[ρ]−E_xc[ρ_bond]−E_xc[ρ_rest]`（通常为小负）
    pub e_cross_xc: f64,
    /// 完整电子能量 `E_elec = E_bond + E_rest + E_cross`
    pub e_elec: f64,
    /// 核-核排斥能 `E_nuc`（整个分子）
    pub e_nuc: f64,
    /// 总能量 `E_total = E_elec + E_nuc`（可与 SCF Done 对照）
    pub e_total: f64,
}

impl BondEnergy {
    /// 键对电子能量的可加贡献：`E_bond + E_cross`
    ///
    /// 注意：这是"电子能量"的份额（大负数），不是键解离能；与总能量可比需再
    /// 按约定分摊核排斥能。
    pub fn bond(&self) -> f64 {
        self.e_bond + self.e_cross
    }
}

/// 键能量计算器
pub struct Calculator<'m> {
    /// 分子引用
    pub mole: &'m Mole,
}

impl<'m> Calculator<'m> {
    /// 创建一个新的 Calculator 实例
    pub fn new(mole: &'m Mole) -> Self {
        Calculator { mole }
    }

    /// 键（两原子）对电子能量的可加贡献 `E_bond + E_cross`
    ///
    /// 完整分解见 [`Calculator::bond_parts`]。
    pub fn bond_energy(&self, atm1: usize, atm2: usize) -> Result<f64, String> {
        Ok(self.bond_parts(atm1, atm2)?.bond())
    }

    /// 对指定原子对的完整电子能量做 bond / rest / cross 分解
    ///
    /// 保留两个原子的分子轨道系数（把其余原子的原子轨道系数置零），得到键密度
    /// `D_bond`，并取补密度 `D_rest = D − D_bond`，再按 KS 能量泛函做加法分解。
    pub fn bond_parts(&self, atm1: usize, atm2: usize) -> Result<BondEnergy, String> {
        energy_v1::bond_parts(self.mole, atm1, atm2)
    }
}
