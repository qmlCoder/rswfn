use crate::base::Mole;
use crate::props::orbtprop::mocv::MocvDest;
use crate::props::orbtprop::obtmat;
use crate::qmcal::ks;
use crate::reader::LogReader;
use crate::rsint::RsInt;
use crate::rsxcs::XcFunc;
use ndarray::{Array2, Axis};

/// π / σ 能量分解结果（闭壳层，杂化泛函）
///
/// 把完整电子能量按密度拆分 `D = D_π + D_σ` 做加法分解。由于 `J`、`K` 对 `D`
/// 线性，库仑/精确交换部分可精确拆开；`E_xc` 非线性，用 `E_xc[ρ_π]`、`E_xc[ρ_σ]`
/// 并把差值归入交叉项 `E_πσ`，从而保证
///
/// `E_π + E_σ + E_πσ = E_elec`（严格成立）。
#[derive(Debug, Clone, Copy)]
pub struct PiEnergy {
    /// π 自能量 `E_π`：π 自身的动能/核吸引、π-π 库仑/交换与 π 的 XC
    pub e_pi: f64,
    /// σ 自能量 `E_σ`
    pub e_sigma: f64,
    /// π-σ 交叉项 `E_πσ`（含 XC 交叉/残差）
    pub e_cross: f64,
    /// 交叉项中的库仑部分 `Tr(D_π J_σ)`（π-σ 静电排斥，通常为正）
    pub e_cross_coul: f64,
    /// 交叉项中的精确交换部分 `−aₓ/2·Tr(D_π K_σ)`（通常为负）
    pub e_cross_exch: f64,
    /// 交叉项中的 XC 残差 `E_xc[ρ]−E_xc[ρ_π]−E_xc[ρ_σ]`（通常为小负）
    pub e_cross_xc: f64,
    /// 完整电子能量 `E_elec = E_π + E_σ + E_πσ`
    pub e_elec: f64,
    /// 核-核排斥能 `E_nuc`
    pub e_nuc: f64,
    /// 总能量 `E_total = E_elec + E_nuc`（可与 SCF Done 对照）
    pub e_total: f64,
}

impl PiEnergy {
    /// π 对电子能量的可加贡献：`E_π + E_πσ`
    ///
    /// 注意：这是"电子能量"的份额（大负数），不是分子总能量；与总能量可比需再加
    /// 上核排斥能的分摊。
    pub fn pi(&self) -> f64 {
        self.e_pi + self.e_cross
    }
}

/// 能量计算器
pub struct Calculator<'m> {
    pub mole: &'m Mole,
}

impl<'m> Calculator<'m> {
    /// 创建一个新的 Calculator 实例
    pub fn new(mole: &'m Mole) -> Self {
        Calculator { mole }
    }

    /// 从 log 文件读取泛函与基组，返回可用的分子与泛函
    ///
    /// 仅支持由高斯 log/out 文件读入的分子。`Mole` 通常已带基组；若缺失，则按
    /// log 中的基组名用 `rsbse` 加载（此时返回 `Some(克隆后的分子)`）。
    pub(crate) fn prepare_log(&self) -> Result<(Option<Mole>, XcFunc), String> {
        let path = self.mole.path.as_deref().ok_or("分子没有来源文件路径，无法读取泛函与基组")?;
        let lower = path.to_lowercase();
        if !(lower.ends_with(".log") || lower.ends_with(".out")) {
            return Err(format!("仅支持高斯 log/out 文件: {path}"));
        }
        let job = LogReader::new(path).read_gau_job().ok_or("无法从 log 文件中解析作业关键词")?;
        let func_name = job.functional.as_deref().ok_or("log 中未识别到交换相关泛函（可能为 HF，或本库不支持的泛函）")?;
        let func = XcFunc::from_name(func_name).ok_or_else(|| format!("log 中的泛函不受支持: {func_name}"))?;

        if self.mole.is_open() {
            return Err("仅支持闭壳层分子".to_string());
        }
        if self.mole.coefs.is_none() {
            return Err("分子没有轨道系数".to_string());
        }

        let owned_mole = if self.mole.basis.is_none() {
            let basis_name = job.basis_set.as_deref().ok_or("log 中未识别到基组名")?;
            let basis = crate::rsbse::from_name(basis_name, &self.mole.atoms).ok_or_else(|| format!("无法加载基组: {basis_name}"))?;
            let mut mole = self.mole.clone();
            mole.basis = Some(basis);
            Some(mole)
        } else {
            None
        };

        Ok((owned_mole, func))
    }

    /// 由 π 投影系数矩阵对完整电子能量做 π/σ/πσ 分解
    ///
    /// `cmat_pi` 为 π 投影系数（`cmat_form`：POCV 为 `"raw"`，MOCV 为 `"sph"`）。
    fn pi_decompose(mole: &Mole, func: XcFunc, cmat_pi: &Array2<f64>, cmat_form: &str) -> Result<PiEnergy, String> {
        let form = mole.basis().int_form();
        let occ = mole.coefs().obt_idxs(true).clone();

        // 完整密度矩阵 D = 2·C_occ·C_occᵀ
        let cmat_full = mole.coefs().get_cmat(form, mole);
        let coef_occ = cmat_full.select(Axis(1), &occ);
        let dmat = coef_occ.dot(&coef_occ.t()) * 2.0;

        // π 投影密度矩阵 D_π = 2·C_occ^π·(C_occ^π)ᵀ
        let mut mole_pi = mole.clone();
        mole_pi.coefs.as_mut().ok_or("分子没有轨道系数")?.set_cmat(cmat_form, cmat_pi, mole);
        let cmat_pi_form = mole_pi.coefs().get_cmat(form, &mole_pi);
        let coef_pi = cmat_pi_form.select(Axis(1), &occ);
        let dmat_pi = coef_pi.dot(&coef_pi.t()) * 2.0;
        // 补密度（保证 D = D_π + D_σ 精确成立）
        let dmat_sig = &dmat - &dmat_pi;

        // 单/双电子积分（form 形式）
        let int_caler = RsInt::new(mole);
        let [_smat, tmat, vmat] = int_caler.mat1e([1, 1, 1], form);
        let hmat = &tmat + &vmat;

        // J、K 对 D 线性：分别对 D_π、D_σ 求，再相加得总。
        // 双电子部分走壳层块引擎（Schwarz 筛选 + 并行）。
        let (j_pi, k_pi, j_sig, k_sig) = {
            let engine = crate::rsint::Engine::new(mole, form, 1.0e-11, 8 * 1024 * 1024 * 1024);
            let (j_pi, k_pi) = engine.jk(&dmat_pi);
            let (j_sig, k_sig) = engine.jk(&dmat_sig);
            (j_pi, k_pi, j_sig, k_sig)
        };
        let j_tot = &j_pi + &j_sig;
        let k_tot = &k_pi + &k_sig;

        // 交换相关能（数值格点积分）：π、σ、总
        let rks = ks::RKS::new(mole, func);
        let xc_pi = rks.xc_eng(&dmat_pi);
        let xc_sig = rks.xc_eng(&dmat_sig);
        let xc_tot = rks.xc_eng(&dmat);
        let e_nuc = rks.nuc_eng();

        let ax = func.exact_x();
        let tr = |a: &Array2<f64>, b: &Array2<f64>| (a * b).sum();

        let e_pi = tr(&dmat_pi, &hmat) + 0.5 * tr(&dmat_pi, &j_pi) - ax * 0.25 * tr(&dmat_pi, &k_pi) + xc_pi;
        let e_sigma = tr(&dmat_sig, &hmat) + 0.5 * tr(&dmat_sig, &j_sig) - ax * 0.25 * tr(&dmat_sig, &k_sig) + xc_sig;
        let e_cross_coul = tr(&dmat_pi, &j_sig);
        let e_cross_exch = -ax * 0.5 * tr(&dmat_pi, &k_sig);
        let e_cross_xc = xc_tot - xc_pi - xc_sig;
        let e_cross = e_cross_coul + e_cross_exch + e_cross_xc;
        let e_elec = tr(&dmat, &hmat) + 0.5 * tr(&dmat, &j_tot) - ax * 0.25 * tr(&dmat, &k_tot) + xc_tot;

        Ok(PiEnergy {
            e_pi,
            e_sigma,
            e_cross,
            e_cross_coul,
            e_cross_exch,
            e_cross_xc,
            e_elec,
            e_nuc,
            e_total: e_elec + e_nuc,
        })
    }

    /// POCV 方法的 π/σ/πσ 能量分解
    pub fn pi_pocv_parts(&self) -> Result<PiEnergy, String> {
        let (owned_mole, func) = self.prepare_log()?;
        let mole = owned_mole.as_ref().map_or(self.mole, |mole| mole);
        let (_, cmat_pocv) = obtmat::Calculator::new(mole).pi_pocv()?;
        Self::pi_decompose(mole, func, &cmat_pocv, "raw")
    }

    /// MOCV 方法的 π/σ/πσ 能量分解
    pub fn pi_mocv_parts(&self) -> Result<PiEnergy, String> {
        let (owned_mole, func) = self.prepare_log()?;
        let mole = owned_mole.as_ref().map_or(self.mole, |mole| mole);
        let (_stms, _atos, cmat_mocv) = obtmat::Calculator::new(mole).pi_mocv(&MocvDest::PiProperty);
        Self::pi_decompose(mole, func, &cmat_mocv, "sph")
    }

    /// 使用 POCV 方法计算的 π 对电子能量的贡献 `E_π + E_πσ`
    pub fn pi_pocv(&self) -> Result<f64, String> {
        Ok(self.pi_pocv_parts()?.pi())
    }

    /// 使用 MOCV 方法计算的 π 对电子能量的贡献 `E_π + E_πσ`
    pub fn pi_mocv(&self) -> Result<f64, String> {
        Ok(self.pi_mocv_parts()?.pi())
    }
}
