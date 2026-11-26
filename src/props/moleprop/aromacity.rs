use rayon::vec;

use crate::base::Mole;
use crate::props::bondprop::order;

/// 芳香性计算器
pub struct Calculator<'m> {
    pub mole: &'m Mole,
}

impl<'m> Calculator<'m> {
    /// 创建一个新的 Calculator 实例
    pub fn new(mole: &'m Mole) -> Self {
        Calculator { mole }
    }

    /// 使用pi键级标准差计算芳香性
    pub fn PISD(&self, rings: Option<&Vec<Vec<usize>>>) -> Result<Vec<f64>, String> {
        let atoms = self.mole.atoms();
        let natm = atoms.len();
        let rings = match rings {
            Some(r) => r.clone(),
            None => vec![atoms.idxs().clone()], // 只有一个环，且为所有原子
        };
        let caler = order::Calculator::new(self.mole);
        let (_dirs, omat) = caler.pi_pocv("mayer")?;
        let mut vals = Vec::new();
        let syms = self.mole.atoms().syms();
        for ring in rings {
            let mut val = 0.0;
            for i in 0..natm - 1 {
                for j in i + 1..natm {
                    if atoms.is_bond(i as usize, j as usize)
                        && ring.contains(&(i as usize))
                        && ring.contains(&(j as usize))
                        && syms[i] != "H" // 排除C-H键
                        && syms[j] != "H"
                    {
                        let order = omat[(i, j)];
                        val += order
                    }
                }
            }
            vals.push(val);
        }
        Ok(vals)
    }
}
