#![allow(non_snake_case)]

use ndarray::Array2;
use std::collections::HashMap;

use crate::base::{Mole, Stm, Vec3};

use crate::props::atomprop;
use crate::props::bondprop;

/// 分子活性（福井函数、自由价）计算器。
///
/// 封装分子参考，提供福井函数和自由价相关的计算方法。
pub struct Calculator<'m> {
    pub mole: &'m Mole,
}

impl<'m> Calculator<'m> {
    /// 创建一个新的活性计算器。
    ///
    /// # Parameters
    /// - `mole` — 分子引用
    pub fn new(mole: &'m Mole) -> Self {
        Calculator { mole }
    }

    /**
    计算原子的福井函数 (Fukui function) 及双描述符 (dual descriptor)。
    # Parameters
    - `mole_n` — 多一个电子的分子（阴离子）
    - `mole_p` — 少一个电子的分子（阳离子）
    - `ctype` — 电荷类型: `"mulliken"` | `"lowdin"` | `"hirshfeld"`
    # Returns
    每个原子 7 个值: `[e0, en, ep, f-, f+, f0, Δf]`，分别对应:
    - `e0` — 中性分子的电子布局
    - `en` — 阴离子的电子布局
    - `ep` — 阳离子的电子布局
    - `f-` = e0 - ep（亲核攻击指数，f⁻, nucleophilic）
    - `f+` = en - e0（亲电攻击指数，f⁺, electrophilic）
    - `f0` = (f⁻ + f⁺) / 2（自由基攻击指数，f⁰, radical）
    - `Δf` = f⁺ - f⁻（双描述符, dual descriptor）
    */
    pub fn fukui(&self, mole_n: &Mole, mole_p: &Mole, ctype: &str) -> Vec<[f64; 7]> {
        let caler_0 = atomprop::charge::Calculator::new(self.mole);
        let caler_n = atomprop::charge::Calculator::new(mole_n);
        let caler_p = atomprop::charge::Calculator::new(mole_p);

        let eles_0 = caler_0.eles(ctype);
        let eles_n = caler_n.eles(ctype);
        let eles_p = caler_p.eles(ctype);

        let nele = eles_n.len();
        let mut vals = Vec::new();
        for i in 0..nele {
            let e0 = eles_0[i];
            let en = eles_n[i];
            let ep = eles_p[i];
            let _fp = en - e0;
            let _fn = e0 - ep;
            let f0 = (_fn + _fp) / 2.0;
            let df = _fp - _fn;
            vals.push([e0, en, ep, _fn, _fp, f0, df]);
        }
        vals
    }

    /// 计算 π 轨道的福井函数及双描述符。
    ///
    /// 通过 π-POÇV 方法得到 π 电子布居，再计算福井函数。
    ///
    /// # Parameters
    /// - `mole_n` — 多一个电子的分子（阴离子）
    /// - `mole_p` — 少一个电子的分子（阳离子）
    /// - `ctype` — 电荷类型: `"mulliken"` | `"lowdin"` | `"hirshfeld"`
    ///
    /// # Returns
    /// 返回 `(dirs, vals)`，其中 `dirs` 是每个原子的 π 方向映射，
    /// `vals` 是每个原子 7 个值: `[e0, en, ep, f⁻, f⁺, f⁰, Δf]`。
    pub fn fukui_pi_pocv(
        &self,
        mole_n: &Mole,
        mole_p: &Mole,
        ctype: &str,
    ) -> Result<(HashMap<usize, Vec3>, Vec<[f64; 7]>), String> {
        let caler_0 = atomprop::charge::Calculator::new(self.mole);
        let caler_n = atomprop::charge::Calculator::new(&mole_n);
        let caler_p = atomprop::charge::Calculator::new(&mole_p);

        let (dirs, eles_n) = caler_n.pi_pocv(ctype)?;
        let (dirs, eles_p) = caler_p.pi_pocv(ctype)?;
        let (dirs, eles_0) = caler_0.pi_pocv(ctype)?;

        let nele = eles_n.len();
        let mut vals = Vec::new();
        for i in 0..nele {
            let e0 = eles_0[i];
            let en = eles_n[i];
            let ep = eles_p[i];
            let _fp = en - e0;
            let _fn = e0 - ep;
            let f0 = (_fn + _fp) / 2.0;
            let df = _fp - _fn;
            vals.push([e0, en, ep, _fn, _fp, f0, df]);
        }
        Ok((dirs, vals))
    }

    pub fn fukui_pi_mocv(
        &self,
        mole_n: &Mole,
        mole_p: &Mole,
        ctype: &str,
    ) -> Result<(HashMap<usize, Stm>, Vec<usize>, Vec<[f64; 7]>), String> {
        let caler_n = atomprop::charge::Calculator::new(&mole_n);
        let caler_p = atomprop::charge::Calculator::new(&mole_p);
        let caler_0 = atomprop::charge::Calculator::new(&self.mole);
        let (stms, atos, eles_0) = caler_0.pi_mocv(ctype)?;
        let (stms, atos, eles_n) = caler_n.pi_mocv(ctype)?;
        let (stms, atos, eles_p) = caler_p.pi_mocv(ctype)?;
        let nele = eles_n.len();
        let mut vals = Vec::new();
        for i in 0..nele {
            let e0 = eles_0[i];
            let en = eles_n[i];
            let ep = eles_p[i];
            let _fp = en - e0;
            let _fn = e0 - ep;
            let f0 = (_fn + _fp) / 2.0;
            let df = _fp - _fn;
            vals.push([e0, en, ep, _fn, _fp, f0, df]);
        }
        Ok((stms, atos, vals))
    }

    /// 计算单个原子在指定方向的福井函数。
    ///
    /// # Parameters
    /// - `atm` — 原子索引
    /// - `dir` — 方向向量
    /// - `mole_n` — 阴离子
    /// - `mole_p` — 阳离子
    /// - `ctype` — 电荷类型
    ///
    /// # Returns
    /// 7 个值: `[e0, en, ep, f⁻, f⁺, f⁰, Δf]`
    pub fn fukui_dir(
        &self,
        atm: usize,
        dir: &Vec3,
        mole_n: &Mole,
        mole_p: &Mole,
        ctype: &str,
    ) -> Result<[f64; 7], String> {
        let caler_0 = atomprop::charge::Calculator::new(self.mole);
        let caler_n = atomprop::charge::Calculator::new(&mole_n);
        let caler_p = atomprop::charge::Calculator::new(&mole_p);

        let mut dirs_map = HashMap::new();
        dirs_map.insert(atm, dir.clone());

        let eles_0 = caler_0.pocv(&dirs_map, false, false, ctype)?;
        let eles_n = caler_n.pocv(&dirs_map, false, false, ctype)?;
        let eles_p = caler_p.pocv(&dirs_map, false, false, ctype)?;

        let e0 = eles_0[atm as usize];
        let en = eles_n[atm as usize];
        let ep = eles_p[atm as usize];
        let _fp = en - e0;
        let _fn = e0 - ep;
        let f0 = (_fn + _fp) / 2.0;
        let df = _fp - _fn;
        Ok([e0, en, ep, _fn, _fp, f0, df])
    }

    /// 计算所有原子的自由价（free valence）。
    ///
    /// # Parameters
    /// - `otype` — 键级类型
    ///
    /// # Returns
    /// 每个原子 3 个值: `[标准价, 实际键级和, 自由价 = 标准价 - 键级和]`
    pub fn freev(&self, otype: &str) -> Vec<[f64; 3]> {
        let caler = bondprop::order::Calculator::new(self.mole);
        let omat = caler.order(otype);
        let natm = self.mole.atoms.len();
        let mut crgs = self.mole.atoms.crgs();
        let mut vals = vec![[0.0; 3]; natm];
        for ai in 0..natm {
            vals[ai][0] = element_valence(crgs[ai] as usize)[0] as f64;
            for aj in self.mole.atoms.nebs(ai) {
                vals[ai][1] += omat[(ai, aj)];
            }
            vals[ai][2] = vals[ai][0] - vals[ai][1]
        }
        vals
    }

    /// 计算所有原子的 π 自由价。
    ///
    /// 通过 π-POÇV 方法得到 π 键级矩阵，再计算 π 自由价。
    ///
    /// # Parameters
    /// - `otype` — 键级类型
    ///
    /// # Returns
    /// 返回 `(dirs, vals)`，其中 `dirs` 是每个原子的 π 方向映射，
    /// `vals` 是每个原子 3 个值: `[标准π价, π键级和, π自由价]`
    pub fn freev_pi(&self, otype: &str) -> Result<(HashMap<usize, Vec3>, Vec<[f64; 3]>), String> {
        let caler = bondprop::order::Calculator::new(&self.mole);
        let (dirs, omat) = caler.pi_pocv(otype)?;
        let natm = self.mole.atoms.len();
        let mut crgs = self.mole.atoms.crgs();
        let mut vals = vec![[0.0; 3]; natm];
        for ai in 0..natm {
            vals[ai][0] = element_valence(crgs[ai] as usize)[1] as f64;
            for aj in self.mole.atoms.nebs(ai) {
                vals[ai][1] += omat[(ai, aj)];
            }
            vals[ai][2] = vals[ai][0] - vals[ai][1]
        }
        Ok((dirs, vals))
    }

    /// 计算单个原子在指定方向的自由价，衡量某个方向的剩余成键能力
    ///
    /// # Parameters
    /// - `atm` — 原子索引
    /// - `dir` — 方向向量
    /// - `nebs` — 相邻原子索引列表
    /// - `otype` — 键级类型
    /// - `mtype` — 方法类型 pocv|mocv
    ///
    /// # Returns
    /// `(标准价, 各相邻原子在此方向的键级, 方向自由价)`
    pub fn freev_dir(
        &self,
        atm: usize,
        dir: &Vec3,
        nebs: &Vec<usize>,
        otype: &str,
        mtype: &str,
    ) -> Result<(f64, Vec<f64>, f64), String> {
        let crg = self.mole.atoms.crgs()[atm];
        let STAND = element_valence(crg as usize)[1];
        let caler = bondprop::order::Calculator::new(&self.mole);
        let vals = caler.bound(atm, dir, nebs, otype, mtype)?;
        let value = STAND - vals.iter().sum::<f64>();
        Ok((STAND, vals, value))
    }
}

/// 返回元素的价数和 π 价数。
///
/// # Parameters
/// - `crg` — 元素核电荷数（原子序数）
///
/// # Returns
/// `[价数, π价数]`
pub fn element_valence(crg: usize) -> [f64; 2] {
    match crg {
        6 => [4.0, 3.0],
        _ => [4.0, 3.0],
    }
}
