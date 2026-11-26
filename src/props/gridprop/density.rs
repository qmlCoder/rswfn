#![allow(non_snake_case)]
use std::cell::RefCell;
use std::sync::OnceLock;

use crate::base::{Atos, Mole};
use crate::datas::rdens;
use crate::props::gridprop::GridVals;
use crate::props::gridprop::wfnfunc;
use crate::maths;
use crate::props::moleprop;
use crate::rsint::EspCalculator;
use ndarray::Array2;

/// 电子密度计算器，支持分子、原子、轨道等不同级别的电子密度计算
///
/// 同时提供静电势的计算（核项 + 电子项，单位 Hartree），见
/// [`Calculator::nuc_esp`]、[`Calculator::ele_esp`] 与 [`Calculator::tot_esp`]，
/// 原理见 `docs/静电势.md`。
pub struct Calculator<'m> {
    pub mole: &'m Mole,
    pub atms: &'m Vec<usize>,
    wfn_caler: wfnfunc::Calculator<'m>,
    mole_pi_pocv: OnceLock<Mole>,
    mole_pi_mocv: OnceLock<Mole>,
    mole_spin: OnceLock<(Mole, Mole)>,
    /// 与格点无关的预计算量（惰性初始化）
    rho_plan: OnceLock<RhoPlan<'m>>,
    /// 静电势的解析积分器（惰性初始化）
    esp_plan: OnceLock<EspCalculator>,
}

/// 与格点无关的预计算量，避免在每个格点上重复计算。
struct RhoPlan<'m> {
    /// car 形式的原子轨道
    atos: &'m Atos,
    /// car 形式的密度矩阵
    dmat: &'m Array2<f64>,
    /// 实际参与计算的基函数在 `atos.cars` 中的索引（升序）
    idxs: Vec<usize>,
    /// `idxs[k]` 所属基函数的中心坐标
    xyzs: Vec<[f64; 3]>,
    /// `idxs[k]` 的截断距离平方，超过该值该基函数恒为 0
    r2cut: Vec<f64>,
}

thread_local! {
    /// 线程内复用的临时缓冲，避免逐格点分配堆内存
    static SIG_BUF: RefCell<Vec<usize>> = RefCell::new(Vec::new());
    static V0_BUF: RefCell<Vec<f64>> = RefCell::new(Vec::new());
    static GV_BUF: RefCell<Vec<GridVals>> = RefCell::new(Vec::new());
}

impl<'m> Calculator<'m> {
    /// 创建一个新的 Calculator 实例
    pub fn new(mole: &'m Mole) -> Self {
        let atms = mole.atoms().idxs();
        let mut wfn_caler = wfnfunc::Calculator::new(mole);
        wfn_caler.atms = atms;
        Calculator {
            mole: mole,
            atms: atms,
            wfn_caler: wfn_caler,
            mole_pi_pocv: OnceLock::new(),
            mole_pi_mocv: OnceLock::new(),
            mole_spin: OnceLock::new(),
            rho_plan: OnceLock::new(),
            esp_plan: OnceLock::new(),
        }
    }

    /// 获取仅包含 π 电子成分的分子对象（懒加载）
    pub fn mole_pi_pocv(&self) -> &Mole {
        self.mole_pi_pocv.get_or_init(|| {
            moleprop::variants::Calculator::new(self.mole)
                .mole_pi_pocv()
                .expect("π轨道投影失败")
        })
    }

    pub fn mole_pi_mocv(&self) -> &Mole {
        self.mole_pi_mocv.get_or_init(|| {
            moleprop::variants::Calculator::new(self.mole)
                .mole_pi_mocv()
                .expect("π轨道投影失败")
        })
    }

    /// 获取自旋分离后的 α 和 β 分子对象（懒加载）
    pub fn mole_spin(&self) -> &(Mole, Mole) {
        self.mole_spin
            .get_or_init(|| moleprop::variants::Calculator::new(self.mole).mole_spin())
    }

    // 计算一个点处的电子密度，这些只计算一个点的函数不需要使用numpy数组传参，但是如果有格点的话就需要了
    #[rustfmt::skip]
    /// 根据分子轨道系数矩阵计算一点处的电子密度及其梯度、Hessian
    pub fn mol_rho_cm(&self,grid: &[f64; 3],level: usize) -> GridVals {
        let mut mol_rho0 = 0.0;
        let mut mol_rho1 = [0.0;3];
        let mut mol_rho2 = [[0.0;3];3];
        let occs = &self.mole.coefs().obt_occs;
        for (obt,occ) in occs.iter().enumerate() {
            let occ = *occ as f64;
            if occ == 0.0 {continue;}
            let wfns=self.wfn_caler.obt_wfn(obt as usize,grid, level);

            mol_rho0 += wfns.val0.powf(2.0)*occ;

            if level >= 1 {
                mol_rho1[0] += 2.0 * wfns.val0 * wfns.val1[0]*occ;
                mol_rho1[1] += 2.0 * wfns.val0 * wfns.val1[1]*occ;
                mol_rho1[2] += 2.0 * wfns.val0 * wfns.val1[2]*occ;
            }
            if level >= 2 {
                mol_rho2[0][0] += 2.0 * (wfns.val1[0].powi(2) + wfns.val0 * wfns.val2[0][0])*occ;
                mol_rho2[1][1] += 2.0 * (wfns.val1[1].powi(2) + wfns.val0 * wfns.val2[1][1])*occ;
                mol_rho2[2][2] += 2.0 * (wfns.val1[2].powi(2) + wfns.val0 * wfns.val2[2][2])*occ;

                mol_rho2[0][1] += 2.0 * (wfns.val1[0] * wfns.val1[1] + wfns.val0 * wfns.val2[0][1])*occ;
                mol_rho2[0][2] += 2.0 * (wfns.val1[0] * wfns.val1[2] + wfns.val0 * wfns.val2[0][2])*occ;
                mol_rho2[1][2] += 2.0 * (wfns.val1[1] * wfns.val1[2] + wfns.val0 * wfns.val2[1][2])*occ;

            }
        }
        if level>=2{
            mol_rho2[1][0] = mol_rho2[0][1];
            mol_rho2[2][0] = mol_rho2[0][2];
            mol_rho2[2][1] = mol_rho2[1][2];
        }
        GridVals::new(mol_rho0, mol_rho1, mol_rho2)
    }

    /// 构建（或复用）与格点无关的预计算量。
    fn rho_plan(&self) -> &RhoPlan<'m> {
        self.rho_plan.get_or_init(|| {
            let atos = self.mole.basis().atos();
            let dmat = self.mole.coefs().get_dmat("car", self.mole);
            let atoms = &self.mole.atoms;
            let cars = &atos.cars;
            let mut idxs = Vec::new();
            let mut xyzs = Vec::new();
            let mut r2cut = Vec::new();
            for (i, ato) in cars.iter().enumerate() {
                if !self.atms.contains(&ato.atm) {
                    continue;
                }
                // `gtf` 中当 alpha * r^2 > 20 时该原高斯基元恒为 0，
                // 因此只要最弥散（alpha 最小）的原高斯基元都为 0，
                // 整个收缩基函数在该点即为 0，可以直接跳过。
                let amin = ato.alps.iter().cloned().fold(f64::INFINITY, f64::min);
                idxs.push(i);
                xyzs.push(atoms.get(ato.atm).xyz);
                r2cut.push(20.0 / amin);
            }
            RhoPlan {
                atos,
                dmat,
                idxs,
                xyzs,
                r2cut,
            }
        })
    }

    /// 根据密度矩阵计算一点处的电子密度及其梯度、Hessian。
    ///
    /// 只在该格点邻域内取值非零的基函数上求和，截断判据与 `gtf` 内部的
    /// `exp(-alpha * r^2)` 截断一致，因此结果与原实现相同；但单点计算量
    /// 只取决于邻域内基函数的数目，而不再随分子整体增大而增长。
    pub fn mol_rho_dm(&self, grid: &[f64; 3], level: usize) -> GridVals {
        let plan = self.rho_plan();
        let [gx, gy, gz] = *grid;
        if level == 0 {
            SIG_BUF
                .with(|sig| V0_BUF.with(|v0| rho0_core(plan, gx, gy, gz, &mut sig.borrow_mut(), &mut v0.borrow_mut())))
        } else {
            SIG_BUF.with(|sig| {
                GV_BUF.with(|wfns| rho_core(plan, gx, gy, gz, level, &mut sig.borrow_mut(), &mut wfns.borrow_mut()))
            })
        }
    }

    // 计算自旋电子密度
    /// 计算自旋电子密度（α 和 β 两个分量）
    pub fn mol_rho_spin(&self, grid: &[f64; 3], level: usize) -> [GridVals; 2] {
        let (mole_a, mole_b) = self.mole_spin();
        let mut cal_a = Calculator::new(mole_a);
        let mut cal_b = Calculator::new(mole_b);
        cal_a.atms = self.atms;
        cal_b.atms = self.atms;
        let val_a = cal_a.mol_rho_dm(grid, level);
        let val_b = cal_b.mol_rho_dm(grid, level);
        [val_a, val_b]
    }

    // 计算所有原子轨道的电子密度
    /// 计算所有原子轨道的电子密度
    pub fn ato_rho(&self, grid: &[f64; 3], level: usize) -> Vec<GridVals> {
        let pmat_car = self.mole.coefs().get_dmat("car", self.mole);
        let nato = pmat_car.shape()[0];
        let ato_wfns = self.wfn_caler.ato_wfn(grid, level);
        let mut ato_rhos = vec![GridVals::non(); nato];
        for i in 0..nato {
            for j in 0..nato {
                if pmat_car[(i, j)] == 0.0 {
                    continue;
                };
                ato_rhos[i].val0 += ato_wfns[i].val0 * ato_wfns[j].val0 * pmat_car[(i, j)];
                if level >= 1 {
                    for k in 0..3 {
                        ato_rhos[i].val1[k] += (ato_wfns[i].val0 * ato_wfns[j].val1[k]
                            + ato_wfns[j].val0 * ato_wfns[i].val1[k])
                            * pmat_car[(i, j)];
                    }
                }
                if level >= 2 {
                    for k in 0..3 {
                        for l in 0..3 {
                            ato_rhos[i].val2[k][l] += (ato_wfns[i].val2[k][l] * ato_wfns[j].val0
                                + ato_wfns[i].val1[k] * ato_wfns[j].val1[l]
                                + ato_wfns[i].val1[l] * ato_wfns[j].val1[k]
                                + ato_wfns[i].val0 * ato_wfns[j].val2[k][l])
                                * pmat_car[(i, j)];
                        }
                    }
                }
            }
        }
        ato_rhos
    }

    // 计算预电子密度
    /// 计算预电子密度（孤立原子密度叠加）
    pub fn pro_rho(&self, grid: &[f64; 3]) -> f64 {
        let natm = self.mole.atoms().len();
        let mut pro_rho = 0.0;
        for atm in self.atms {
            let atom = self.mole.atoms().get(*atm);
            let xyz = atom.xyz;
            let dist = maths::calc_dist(&xyz, &grid);
            pro_rho += rdens::get_rad_dens_v1(atom.nuc, dist);
        }
        pro_rho
    }

    #[rustfmt::skip]
    /// 计算指定分子轨道的电子密度
    pub fn obt_rho(&self, grid: &[f64; 3], obt: usize, level: usize) -> GridVals {
        let mut obt_rho0 = 0.0;
        let mut obt_rho1 = [0.0; 3];
        let mut obt_rho2 = [[0.0; 3]; 3];
        let occ = self.mole.coefs().obt_occs[obt as usize];
        // 循环所有的分子轨道
        // 先构建该分子轨道的波函数
        let obt_wfns = self.wfn_caler.obt_wfn(obt, grid, level);
        obt_rho0 += obt_wfns.val0.powf(2.0) * (occ as f64);

        if level >= 1 {
            obt_rho1[0] = 2.0 * obt_wfns.val0 * obt_wfns.val1[0] * (occ as f64);
            obt_rho1[1] = 2.0 * obt_wfns.val0 * obt_wfns.val1[1] * (occ as f64);
            obt_rho1[2] = 2.0 * obt_wfns.val0 * obt_wfns.val1[2] * (occ as f64);
        }
        if level >= 2 {
            obt_rho2[0][0] = 2.0 * (obt_wfns.val1[0] * obt_wfns.val1[0] + obt_wfns.val0 * obt_wfns.val2[0][0]) * (occ as f64);
            obt_rho2[1][1] = 2.0 * (obt_wfns.val1[1] * obt_wfns.val1[1] + obt_wfns.val0 * obt_wfns.val2[1][1]) * (occ as f64);
            obt_rho2[2][2] = 2.0 * (obt_wfns.val1[2] * obt_wfns.val1[2] + obt_wfns.val0 * obt_wfns.val2[2][2]) * (occ as f64);
            obt_rho2[0][1] = 2.0 * (obt_wfns.val1[0] * obt_wfns.val1[1] + obt_wfns.val0 * obt_wfns.val2[0][1]) * (occ as f64);
            obt_rho2[0][2] = 2.0 * (obt_wfns.val1[0] * obt_wfns.val1[2] + obt_wfns.val0 * obt_wfns.val2[0][2]) * (occ as f64);
            obt_rho2[1][2] = 2.0 * (obt_wfns.val1[1] * obt_wfns.val1[2] + obt_wfns.val0 * obt_wfns.val2[1][2]) * (occ as f64);
            obt_rho2[1][0] = obt_rho2[0][1];
            obt_rho2[2][0] = obt_rho2[0][2];
            obt_rho2[2][1] = obt_rho2[1][2];
        }
        GridVals::new(obt_rho0, obt_rho1, obt_rho2)
    }

    // Reduced Density Gradient
    /// 计算约化密度梯度（Reduced Density Gradient）
    pub fn RDG(&self, grid: &[f64; 3]) -> f64 {
        let RDG_maxrho = 0.05;
        let grid_vals = self.mol_rho_dm(grid, 1);
        let rho = grid_vals.val0.max(1e-6);
        let [gx, gy, gz] = grid_vals.val1;
        let norm = (gx.powi(2) + gy.powi(2) + gz.powi(2)).sqrt();
        // let coef=1.0/(2.0*(3.0*PI.powi(2)).powf(1.0/3.0));
        let coef = 0.1616204596739955;
        if rho > RDG_maxrho && RDG_maxrho != 0.0 {
            return 100.0;
        } else {
            return coef * norm / rho.powf(4.0 / 3.0);
        }
    }

    /// 计算 IRI（Interaction Region Indicator）
    pub fn IRI(&self, grid: &[f64; 3], alp: Option<f64>) -> f64 {
        let grid_vals = self.mol_rho_dm(grid, 1);
        let alp = alp.unwrap_or(1.1);
        let rho = grid_vals.val0.max(1e-6);
        let [gx, gy, gz] = grid_vals.val1;
        let norm = (gx.powi(2) + gy.powi(2) + gz.powi(2)).sqrt();
        norm / rho.powf(alp)
    }

    /// 计算电子密度的拉普拉斯量（∇²ρ）
    pub fn laplace(&self, grid: &[f64; 3]) -> f64 {
        let grid_vals = self.mol_rho_dm(grid, 2);
        let mut val = 0.0;

        for i in 0..3 {
            val += grid_vals.val2[i][i];
        }
        val
    }

    /// 计算 sign(λ₂)ρ，用于 NCI 分析
    pub fn signl2rho(&self, grid: &[f64; 3]) -> f64 {
        let grid_vals = self.mol_rho_dm(grid, 2);
        let hessian =
            Array2::from_shape_vec((3, 3), grid_vals.val2.iter().cloned().flatten().collect()).expect("矩阵转换失败");
        let eig_vals = maths::eigh(&hessian).0;
        let lambda2 = eig_vals[1];
        let singl2 = if lambda2 < 0.0 {
            -1.0
        } else if lambda2 > 0.0 {
            1.0
        } else {
            0.0
        };
        singl2 * grid_vals.val0
    }

    /// 计算 π 电子密度 (POCV)
    pub fn pi_rho_pocv(&self, grid: &[f64; 3], level: usize) -> GridVals {
        let mole_pocv = self.mole_pi_pocv();
        let rho_caler = Calculator::new(mole_pocv);
        rho_caler.mol_rho_dm(grid, level)
    }

    /// 计算 π 电子密度（MOCV）
    pub fn pi_rho_mocv(&self, grid: &[f64; 3], level: usize) -> GridVals {
        let mole_mocv = self.mole_pi_mocv();
        let rho_caler = Calculator::new(mole_mocv);
        rho_caler.mol_rho_dm(grid, level)
    }

    // 变形电子密度
    /// 计算变形电子密度（Δρ = ρ_promolecular - ρ_molecular）
    pub fn del_rho(&self, grid: &[f64; 3]) -> f64 {
        self.pro_rho(grid) - self.mol_rho_dm(grid, 0).val0
    }

    // 福井函数
    /// 计算福井函数（Fukui function）：f⁻, f⁺, f⁰, Δf
    pub fn fukui(&self, grid: &[f64; 3], moleN: &Mole, moleP: &Mole) -> [f64; 4] {
        let calerN = Calculator::new(moleN); // N+1个电子的分子
        let calerP = Calculator::new(moleP); // N-1个电子的分子
        let rho0 = self.mol_rho_dm(grid, 0).val0;
        let rhoN = calerN.mol_rho_dm(grid, 0).val0;
        let rhoP = calerP.mol_rho_dm(grid, 0).val0;
        let fN = rhoN - rho0;
        let fP = rho0 - rhoP;
        let f0 = (rhoN - rhoP) / 2.0;
        let df = fP - fN;
        [fN, fP, f0, df]
    }

    /// 计算一点处的拉格朗日动能密度，总量以及算个方向的分量
    pub fn lagkin(&self, grid: &[f64; 3]) -> [f64; 4] {
        let occs = &self.mole.coefs().obt_occs;
        let mut vals = [0.0; 4];
        for (obt, occ) in occs.iter().enumerate() {
            if *occ == 0.0 {
                continue;
            }
            let [dx, dy, dz] = self.wfn_caler.obt_wfn(obt as usize, grid, 1).val1;
            let occ = *occ as f64;
            let gx = occ * dx.powi(2);
            let gy = occ * dy.powi(2);
            let gz = occ * dz.powi(2);
            vals[1] += gx;
            vals[2] += gy;
            vals[3] += gz;
            vals[0] += gx + gy + gz;
        }
        for val in vals.iter_mut() {
            *val /= 2.0;
        }
        vals
    }

    /// 计算一点处的哈密顿动能密度，总量以及算个方向的分量
    pub fn hamkin(&self, grid: &[f64; 3]) -> f64 {
        let occs = &self.mole.coefs().obt_occs;
        let mut val = 0.0;
        for (obt, occ) in occs.iter().enumerate() {
            if *occ == 0.0 {
                continue;
            }
            let occ = *occ as f64;
            let obt_wfns = self.wfn_caler.obt_wfn(obt as usize, grid, 2);
            let kx = occ * obt_wfns.val0 * obt_wfns.val2[0][0];
            let ky = occ * obt_wfns.val0 * obt_wfns.val2[1][1];
            let kz = occ * obt_wfns.val0 * obt_wfns.val2[2][2];
            val += kx + ky + kz;
        }
        val * -0.5
    }

    /// 动能密度（kinetic energy density）
    pub fn KED(&self, grid: &[f64; 3]) -> f64 {
        let obts = self.mole.coefs().obt_idxs(true); // 占据分子轨道索引
        let mut val = 0.0;
        let occs = &self.mole.coefs().obt_occs;
        for obt in obts.iter() {
            let obt_wfns = self.wfn_caler.obt_wfn(*obt, grid, 1);
            let [gx, gy, gz] = obt_wfns.val1;
            val += (gx.powi(2) + gy.powi(2) + gz.powi(2)) * occs[*obt];
        }
        val
    }

    /// 电子局域化函数（Electron Localization Function）
    pub fn ELF(&self, grid: &[f64; 3]) -> f64 {
        let open = self.mole.is_open();
        let fcs = [4.557799872345596, 2.871234000188191]; // Thomas-Fermi constants =  (3/10)*(6*Pi^2)**(2/3) (3/10)*(3*Pi^2)**(2/3)
        let mut ked = self.KED(grid) / 2.0;
        let alp = 1.6666666666666; // 5/3
        let (Dr, D0) = if open {
            let [rhos_a, rhos_b] = self.mol_rho_spin(grid, 1);
            let Dr = ked - (rhos_a.grad_norm2() / rhos_a.val0 + rhos_b.grad_norm2() / rhos_b.val0) / 8.0;
            let D0 = fcs[0] * (rhos_a.val0.powf(alp) + rhos_b.val0.powf(alp));
            (Dr, D0)
        } else {
            let rhos = self.mol_rho_dm(grid, 1);
            let Dr = ked - rhos.grad_norm2() / rhos.val0 / 8.0;
            let D0 = fcs[1] * rhos.val0.powf(alp);
            (Dr, D0)
        };
        ked += 1e-5;
        let elf = 1.0 / (1.0 + (Dr / D0).powi(2));
        elf
    }

    ///  Localized orbital locator (LOL)
    pub fn LOL(&self, grid: &[f64; 3]) -> f64 {
        let fcs = [4.557799872345596, 2.871234000188191]; // Thomas-Fermi constants =  (3/10)*(6*Pi^2)**(2/3) (3/10)*(3*Pi^2)**(2/3)
        let alp = 1.6666666666666; // 5/3
        let ked = self.KED(grid) / 2.0;
        let open = self.mole.is_open();
        let D0 = if open {
            let [rhos_a, rhos_b] = self.mol_rho_spin(grid, 1);
            fcs[0] * (rhos_a.val0.powf(alp) + rhos_b.val0.powf(alp))
        } else {
            let rhos = self.mol_rho_dm(grid, 1);
            fcs[1] * rhos.val0.powf(alp)
        };
        let tau = D0 / ked;
        let lol = tau / (1.0 + tau);
        lol
    }

    /// 计算信息熵密度（Information Entropy Density）
    pub fn infoentro(&self, grid: &[f64; 3]) -> f64 {
        let rho0 = self.mol_rho_dm(grid, 0).val0;
        let ncrg = self.mole.ncrg() as f64;
        let rho = rho0 / ncrg;
        if rho < 1e-100 { 0.0 } else { -rho * f64::ln(rho) }
    }

    /// 计算活化密度（Activation Density）
    pub fn activate(&self, grid: &[f64; 3]) -> f64 {
        let rho = self.mol_rho_dm(grid, 0).val0;
        rho * f64::exp(-1.0 * (rho - 0.4))
    }

    /// 计算电子密度的曲率（Curvature）
    pub fn curvature(&self, grid: &[f64; 3]) -> f64 {
        // 计算电子密度的曲率
        let rhos = self.mol_rho_dm(grid, 2);
        let [dx, dy, dz] = rhos.val1;
        let f = Array2::from_shape_vec((1, 3), vec![dx, dy, dz]).expect("矩阵生成失败");
        let norm2 = (dx.powi(2) + dy.powi(2) + dz.powi(2)).powf(0.50000000); // 梯度模长
        let norm3 = (dx.powi(3) + dy.powi(3) + dz.powi(3)).powf(0.33333333); // 梯度模长
        let xx = rhos.val2[0][0];
        let yy = rhos.val2[1][1];
        let zz = rhos.val2[2][2];
        let tr = xx + yy + zz; // Hessian矩阵的迹
        let h = Array2::from_shape_vec((3, 3), rhos.val2.iter().cloned().flatten().collect()).expect("矩阵生成失败");
        let ft = f.t();
        let mul = f.dot(&h).dot(&ft)[(0, 0)];
        (norm2 * tr - mul) / norm3
    }

    // ===== 静电势（ESP） =====

    /// 计算指定点处由原子核产生的静电势 `Σ_A Z_A/|r-R_A|`
    ///
    /// 使用有效核电荷（`atoms.crgs()`），因此赝势（ECP）体系也自洽。约定与 Multiwfn 一致：
    /// 查询点落在原子核上（含浮点误差范围内的重合）时返回 `1e3`，避免出现无穷大。
    pub fn nuc_esp(
        &self,
        grid: &[f64; 3], // 请求点（要计算的点）的坐标
    ) -> f64 {
        let atoms = self.mole.atoms();
        let crgs = atoms.crgs();
        let xyzs = atoms.xyzs();
        let mut val = 0.0;
        for a in 0..crgs.len() {
            let dist = maths::calc_dist(grid, &xyzs[a]);
            if dist < 1e-6 {
                // 与 Multiwfn 相同：落在原子核上（含浮点误差范围内的重合）时给一个有限大值
                return 1e3;
            }
            val += crgs[a] / dist;
        }
        val
    }

    /// 计算指定点处由电子产生的静电势 `-∫ρ(r')/|r-r'|d³r'`
    ///
    /// 由密度矩阵的解析高斯积分直接得到（见 [`EspCalculator`]），不做实空间格点积分，
    /// 因此没有格点积分误差，也不存在 `1/|r-r'|` 的自作用奇点。
    /// 结果为负值，已包含电子负号，可直接与 [`Calculator::nuc_esp`] 相加。
    pub fn ele_esp(
        &self,
        grid: &[f64; 3], // query position 请求点（要计算的点）的坐标
    ) -> f64 {
        let plan = self.esp_plan.get_or_init(|| {
            let dmat = self.mole.coefs().get_dmat("car", self.mole);
            EspCalculator::new(self.mole, dmat)
        });
        plan.val(grid)
    }

    /// 计算指定点处的总静电势（核项 + 电子项）
    ///
    /// `V(r) = Σ_A Z_A/|r-R_A| - ∫ρ(r')/|r-r'|d³r'`，单位 Hartree。
    /// 原理与推导见 `docs/静电势.md`。
    pub fn tot_esp(&self, grid: &[f64; 3]) -> f64 {
        self.nuc_esp(grid) + self.ele_esp(grid)
    }
}

/// 根据密度矩阵计算电子密度（自由函数版本）
pub fn mol_rho_dm(
    ato_atms: &Vec<usize>,
    dmat_car: &Array2<f64>,
    ato_wfn0: &Vec<f64>,
    ato_wfn1: &Vec<[f64; 3]>,
    ato_wfn2: &Vec<[[f64; 3]; 3]>,
    atms: &Vec<usize>,
    level: usize,
) -> (f64, [f64; 3], [[f64; 3]; 3]) {
    let mut mol_rho0 = 0.0;
    let mut mol_rho1 = [0.0; 3];
    let mut mol_rho2 = [[0.0; 3]; 3];
    let nmat = ato_wfn0.len();

    for i in 0..nmat {
        let atmi = ato_atms[i];
        if !atms.contains(&atmi) {
            continue;
        }
        for j in i..nmat {
            let atmj = ato_atms[j];
            if !atms.contains(&atmj) {
                // 如果不是需要计算的原子则跳过
                continue;
            }

            let d = dmat_car[(i, j)];
            // let d=d*self.dmat_fix[(atmi as usize,atmj as usize)]; // 给密度矩阵加矫正

            if d == 0.0 {
                continue;
            }
            let k = if i == j { 1.0 } else { 2.0 };
            mol_rho0 += ato_wfn0[i] * ato_wfn0[j] * d * k;

            // 一阶导数：dx
            if level >= 1 {
                mol_rho1[0] += (ato_wfn0[i] * ato_wfn1[j][0] + ato_wfn0[j] * ato_wfn1[i][0]) * d * k;
                mol_rho1[1] += (ato_wfn0[i] * ato_wfn1[j][1] + ato_wfn0[j] * ato_wfn1[i][1]) * d * k;
                mol_rho1[2] += (ato_wfn0[i] * ato_wfn1[j][2] + ato_wfn0[j] * ato_wfn1[i][2]) * d * k;
            }

            // 二阶导数：dxx
            if level >= 2 {
                mol_rho2[0][0] += (ato_wfn2[i][0][0] * ato_wfn0[j]
                    + 2.0 * ato_wfn1[i][0] * ato_wfn1[j][0]
                    + ato_wfn0[i] * ato_wfn2[j][0][0])
                    * d
                    * k;
                mol_rho2[1][1] += (ato_wfn2[i][1][1] * ato_wfn0[j]
                    + 2.0 * ato_wfn1[i][1] * ato_wfn1[j][1]
                    + ato_wfn0[i] * ato_wfn2[j][1][1])
                    * d
                    * k;
                mol_rho2[2][2] += (ato_wfn2[i][2][2] * ato_wfn0[j]
                    + 2.0 * ato_wfn1[i][2] * ato_wfn1[j][2]
                    + ato_wfn0[i] * ato_wfn2[j][2][2])
                    * d
                    * k;
                mol_rho2[0][1] += (ato_wfn2[i][0][1] * ato_wfn0[j]
                    + ato_wfn1[i][0] * ato_wfn1[j][1]
                    + ato_wfn1[i][1] * ato_wfn1[j][0]
                    + ato_wfn0[i] * ato_wfn2[j][0][1])
                    * d
                    * k;
                mol_rho2[0][2] += (ato_wfn2[i][0][2] * ato_wfn0[j]
                    + ato_wfn1[i][0] * ato_wfn1[j][2]
                    + ato_wfn1[i][2] * ato_wfn1[j][0]
                    + ato_wfn0[i] * ato_wfn2[j][0][2])
                    * d
                    * k;
                mol_rho2[1][2] += (ato_wfn2[i][1][2] * ato_wfn0[j]
                    + ato_wfn1[i][1] * ato_wfn1[j][2]
                    + ato_wfn1[i][2] * ato_wfn1[j][1]
                    + ato_wfn0[i] * ato_wfn2[j][1][2])
                    * d
                    * k;
            }
        }
    }
    if level >= 2 {
        mol_rho2[1][0] = mol_rho2[0][1];
        mol_rho2[2][0] = mol_rho2[0][2];
        mol_rho2[2][1] = mol_rho2[1][2];
    }
    (mol_rho0, mol_rho1, mol_rho2)
}

/// level=0 的电子密度：只保留显著基函数，且只计算 val0。
fn rho0_core(plan: &RhoPlan<'_>, gx: f64, gy: f64, gz: f64, sig: &mut Vec<usize>, v0: &mut Vec<f64>) -> GridVals {
    sig.clear();
    v0.clear();
    let cars = &plan.atos.cars;
    for k in 0..plan.idxs.len() {
        let [ax, ay, az] = plan.xyzs[k];
        let (dx, dy, dz) = (gx - ax, gy - ay, gz - az);
        if dx * dx + dy * dy + dz * dz > plan.r2cut[k] {
            continue;
        }
        let ato = &cars[plan.idxs[k]];
        let (a0, _, _) = wfnfunc::cgf(&[dx, dy, dz], ato.lmn(), &ato.alps, ato.coes_norm(), 0);
        sig.push(k);
        v0.push(a0);
    }
    let mut rho0 = 0.0;
    for a in 0..sig.len() {
        let i = plan.idxs[sig[a]];
        let va = v0[a];
        if va == 0.0 {
            continue;
        }
        rho0 += va * va * plan.dmat[(i, i)];
        for b in (a + 1)..sig.len() {
            let d = plan.dmat[(i, plan.idxs[sig[b]])];
            if d == 0.0 {
                continue;
            }
            rho0 += va * v0[b] * d * 2.0;
        }
    }
    GridVals::new(rho0, [0.0; 3], [[0.0; 3]; 3])
}

/// level>=1 的电子密度及其一阶、二阶导数：同样只保留显著基函数。
#[rustfmt::skip]
fn rho_core(
    plan: &RhoPlan<'_>,
    gx: f64,
    gy: f64,
    gz: f64,
    level: usize,
    sig: &mut Vec<usize>,
    wfns: &mut Vec<GridVals>,
) -> GridVals {
    sig.clear();
    wfns.clear();
    let cars = &plan.atos.cars;
    for k in 0..plan.idxs.len() {
        let [ax, ay, az] = plan.xyzs[k];
        let (dx, dy, dz) = (gx - ax, gy - ay, gz - az);
        if dx * dx + dy * dy + dz * dz > plan.r2cut[k] {
            continue;
        }
        let ato = &cars[plan.idxs[k]];
        let (a0, a1, a2) = wfnfunc::cgf(&[dx, dy, dz], ato.lmn(), &ato.alps, ato.coes_norm(), level);
        sig.push(k);
        wfns.push(GridVals::new(a0, a1, a2));
    }
    let mut rho0 = 0.0;
    let mut rho1 = [0.0; 3];
    let mut rho2 = [[0.0; 3]; 3];
    for a in 0..sig.len() {
        let i = plan.idxs[sig[a]];
        let wi = wfns[a];
        for b in a..sig.len() {
            let d = plan.dmat[(i, plan.idxs[sig[b]])];
            if d == 0.0 {
                continue;
            }
            let k = if a == b { 1.0 } else { 2.0 };
            let wj = wfns[b];
            rho0 += wi.val0 * wj.val0 * d * k;
            if level >= 1 {
                rho1[0] += (wi.val0 * wj.val1[0] + wj.val0 * wi.val1[0]) * d * k;
                rho1[1] += (wi.val0 * wj.val1[1] + wj.val0 * wi.val1[1]) * d * k;
                rho1[2] += (wi.val0 * wj.val1[2] + wj.val0 * wi.val1[2]) * d * k;
            }
            if level >= 2 {
                rho2[0][0] += (wi.val2[0][0] * wj.val0 + 2.0 * wi.val1[0] * wj.val1[0] + wi.val0 * wj.val2[0][0]) * d * k;
                rho2[1][1] += (wi.val2[1][1] * wj.val0 + 2.0 * wi.val1[1] * wj.val1[1] + wi.val0 * wj.val2[1][1]) * d * k;
                rho2[2][2] += (wi.val2[2][2] * wj.val0 + 2.0 * wi.val1[2] * wj.val1[2] + wi.val0 * wj.val2[2][2]) * d * k;
                rho2[0][1] += (wi.val2[0][1] * wj.val0 + wi.val1[0] * wj.val1[1] + wi.val1[1] * wj.val1[0] + wi.val0 * wj.val2[0][1]) * d * k;
                rho2[0][2] += (wi.val2[0][2] * wj.val0 + wi.val1[0] * wj.val1[2] + wi.val1[2] * wj.val1[0] + wi.val0 * wj.val2[0][2]) * d * k;
                rho2[1][2] += (wi.val2[1][2] * wj.val0 + wi.val1[1] * wj.val1[2] + wi.val1[2] * wj.val1[1] + wi.val0 * wj.val2[1][2]) * d * k;
            }
        }
    }
    if level >= 2 {
        rho2[1][0] = rho2[0][1];
        rho2[2][0] = rho2[0][2];
        rho2[2][1] = rho2[1][2];
    }
    GridVals::new(rho0, rho1, rho2)
}
