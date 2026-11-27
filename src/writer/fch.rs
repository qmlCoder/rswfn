use ndarray::{Array2, s};
use std::collections::HashMap;

use crate::base::{Atoms, Basis, Coefs, Mole, Shl};
use crate::datas::elems;
use crate::utils;
use crate::writer::Writer;

/// Fch格式文件(.fch)写入器
pub struct FchWriter {
    /// 分子电荷
    pub charge: isize,
    /// 自旋多重度
    pub multip: usize,
    /// 总能量（Hartree），写入文件的 Total Energy 字段
    pub energy: f64,
    /// 原子信息
    pub atoms: Atoms,
    /// 基组信息
    pub basis: Basis,
    /// 波函数信息
    pub coefs: Coefs,
}

impl FchWriter {
    /// 创建FchWriter实例
    pub fn new() -> Self {
        Self {
            charge: 0,
            multip: 1,
            energy: 0.0,
            atoms: Atoms::new(&vec![], &vec![], &vec![]),
            basis: Basis::new(vec![], None),
            coefs: Coefs::new(vec![], vec![], vec![], vec![], vec![], Array2::zeros((0, 0))),
        }
    }

    /// 从Mole构造FchWriter
    pub fn from_mole(mole: &Mole) -> Self {
        Self {
            charge: mole.charge(),
            multip: mole.multip(),
            energy: 0.0,
            atoms: mole.atoms().clone(),
            basis: mole.basis().clone(),
            coefs: mole.coefs().clone(),
        }
    }

    #[rustfmt::skip]
    /// 生成fch文件内容字符串
    pub fn build(&self)->String {
        let mut text = String::new();
        // 按基组实际的壳层形式写入系数：l >= 2 的形式决定系数矩阵是球谐还是笛卡尔
        let form = self.basis.int_form();
        let basis = self.basis.unify_form(form);
        let mole = Mole::new(self.atoms.clone(), Some(basis), Some(self.coefs.clone()));
        let atoms=&self.atoms;
        let natm=atoms.len();
        let nele=mole.coefs().obt_eles(&mole);
        let syms=atoms.syms();
        let xyzs=atoms.xyzs();
        let cmat=mole.coefs().get_cmat(form, &mole);
        let [nato,ncol]=cmat.shape().try_into().expect("获取分子轨道形状失败");
        // 开壳层时系数矩阵为 nato × 2·nato（α 块在前、β 块在后）
        let open=ncol>nato;
        // "Number of independent functions" 指的是每个自旋的分子轨道数
        let nobt=if open{ncol/2}else{ncol};
        // let nbas=self.basis.map_ato().len();
        let nbas=nato;

        text.push_str("generate by rswfn\n");
        text.push_str(&format!("Number of atoms                            I{:>17}\n",natm));
        text.push_str(&format!("Charge                                     I{:>17}\n",self.charge));
        text.push_str(&format!("Multiplicity                               I{:>17}\n",self.multip));
        text.push_str(&format!("Number of electrons                        I{:>17}\n",nele[0] + nele[1]));
        text.push_str(&format!("Number of alpha electrons                  I{:>17}\n",nele[0]));
        text.push_str(&format!("Number of beta electrons                   I{:>17}\n",nele[1]));
        text.push_str(&format!("Number of basis functions                  I{:>17}\n",nbas));
        text.push_str(&format!("Number of independent functions            I{:>17}\n",nobt));
        text.push_str(&format!("Atomic numbers                             I   N={:>12}\n",natm));
        let atms:Vec<isize>=syms.iter().map(|sym|elems::sym2atm(sym) as isize).collect();
        text.push_str(&i2str(&atms));
        text.push_str(&format!("Nuclear charges                            R   N={:>12}\n",natm));
        let crgs:Vec<f64>=atms.iter().map(|atm|*atm as f64).collect();
        text.push_str(&f2str(&crgs));
        text.push_str(&format!("Current cartesian coordinates              R   N={:>12}\n",natm*3));
        let xyz_vec:Vec<f64>=xyzs.iter().flat_map(|xyz|xyz.iter().copied()).collect();
        text.push_str(&f2str(&xyz_vec));

        text.push_str("Translation vectors                        R   N=           0\n");
        // 基组信息
        let (nshl,shl_types,shl_nbass,shl_atoms,pm_alps,pm_coes,sp_coes)=basis_info(&mole);
        let shl_nbass=shl_nbass.iter().map(|v|*v as isize).collect();
        let shl_atoms:Vec<isize>=shl_atoms.iter().map(|v|*v as isize).collect();
        text.push_str(&format!("Shell types                                I   N={:>12}\n",nshl));
        text.push_str(&i2str(&shl_types));
        text.push_str(&format!("Number of primitives per shell             I   N={:>12}\n",nshl));
        text.push_str(&i2str(&shl_nbass));
        text.push_str(&format!("Shell to atom map                          I   N={:>12}\n",nshl));
        text.push_str(&i2str(&shl_atoms));
        text.push_str(&format!("Primitive exponents                        R   N={:>12}\n",pm_alps.len()));
        text.push_str(&f2str(&pm_alps));
        text.push_str(&format!("Contraction coefficients                   R   N={:>12}\n",pm_coes.len()));
        text.push_str(&f2str(&pm_coes));
        text.push_str(&format!("P(S=P) Contraction coefficients            R   N={:>12}\n",sp_coes.len()));
        text.push_str(&f2str(&sp_coes));
        text.push_str(&format!("Coordinates of each shell                  R   N={:>12}\n",nshl*3));
        let xyz_vec:Vec<f64>=shl_atoms.iter().flat_map(|atm|xyzs[*atm as usize-1].iter().copied()).collect();
        text.push_str(&f2str(&xyz_vec));
        text.push_str(&format!("Virial Ratio                               R{:>27.15E}\n",0.0));
        text.push_str(&format!("Total Energy                               R{:>27.15E}\n",self.energy));

        // 分子轨道信息
        let obt_engs=&self.coefs.obt_engs;
        if open{
            // 开壳层：α/β 轨道能分别写入，各 nobt 个
            let engs_a=obt_engs[..nobt].to_vec();
            let engs_b=obt_engs[nobt..2*nobt].to_vec();
            text.push_str(&format!("Alpha Orbital Energies                     R   N={:>12}\n",engs_a.len()));
            text.push_str(&f2str(&engs_a));
            text.push_str(&format!("Beta Orbital Energies                      R   N={:>12}\n",engs_b.len()));
            text.push_str(&f2str(&engs_b));
            let coefs_a=cmat.slice(s![..,..nobt]).t().flatten().to_vec();
            let coefs_b=cmat.slice(s![..,nobt..]).t().flatten().to_vec();
            text.push_str(&format!("Alpha MO coefficients                      R   N={:>12}\n",coefs_a.len()));
            text.push_str(&f2str(&coefs_a));
            text.push_str(&format!("Beta MO coefficients                       R   N={:>12}\n",coefs_b.len()));
            text.push_str(&f2str(&coefs_b));
        }else{
            let coefs=cmat.t().flatten().to_vec();
            text.push_str(&format!("Alpha Orbital Energies                     R   N={:>12}\n",obt_engs.len()));
            text.push_str(&f2str(obt_engs));
            text.push_str(&format!("Alpha MO coefficients                      R   N={:>12}\n",coefs.len()));
            text.push_str(&f2str(&coefs));
        }

        // Total SCF Density（下三角，含对角元）
        let dmat=mole.coefs().get_dmat(form, &mole);
        let mut dvec=Vec::new();
        for i in 0..nato{
            for j in 0..=i{
                dvec.push(dmat[(i,j)]);
            }
        }
        text.push_str(&format!("Total SCF Density                          R   N={:>12}\n",dvec.len()));
        text.push_str(&f2str(&dvec));
        text
    }

    /// 保存到文件
    pub fn save(&self, path: &str) {
        let text = self.build();
        std::fs::write(path, text).expect("Unable to write file");
    }
}

fn i2str(nums: &Vec<isize>) -> String {
    let mut text = String::new();
    for (i, num) in nums.iter().enumerate() {
        text.push_str(&format!("{:>12}", num));
        if (i + 1) % 6 == 0 {
            text.push_str("\n");
        }
    }
    if nums.len() % 6 != 0 {
        text.push_str("\n");
    }
    text
}

fn f2str(nums: &Vec<f64>) -> String {
    let mut text = String::new();
    for (i, num) in nums.iter().enumerate() {
        let fort = utils::fort_form(*num, 8);
        text.push_str(&format!("{:>16}", fort));

        if (i + 1) % 5 == 0 {
            text.push_str("\n");
        }
    }
    if nums.len() % 5 != 0 {
        text.push_str("\n");
    }
    text
}

struct ShlInfo {
    atom: usize,
    alps: Vec<f64>,
    coes: Vec<f64>,
    pm_alps: Vec<f64>,
    pm_coes: Vec<f64>,
    sp_coes: Vec<f64>,
    angs: Vec<usize>,
}

impl ShlInfo {
    pub fn blank() -> Self {
        Self {
            atom: 0,
            alps: Vec::new(),
            coes: Vec::new(),
            pm_alps: Vec::new(),
            pm_coes: Vec::new(),
            sp_coes: Vec::new(),
            angs: Vec::new(),
        }
    }
}

fn basis_info(mole: &Mole) -> (usize, Vec<isize>, Vec<usize>, Vec<usize>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let basis = mole.basis();
    let mut shl_types = Vec::new();
    let mut shl_nbass = Vec::new();
    let mut shl_atoms = Vec::new();
    let mut pm_alps = Vec::new();
    let mut pm_coes = Vec::new();
    let sp_coes = Vec::new();
    let nshl = basis.shls.len();
    for i in 0..nshl {
        let shl = basis.shls.get(i);
        shl_types.push(shl_type(&shl));
        shl_nbass.push(shl.alps.len());
        shl_atoms.push(shl.atm + 1); // Gaussian 文件中原子从 1 开始编号
        pm_alps.extend(shl.alps.clone());
        pm_coes.extend(shl.coes.clone());
    }
    (nshl, shl_types, shl_nbass, shl_atoms, pm_alps, pm_coes, sp_coes)
}

/// 壳层类型：l >= 2 的球谐壳层用负数表示（与 Gaussian 约定一致）
fn shl_type(shl: &Shl) -> isize {
    if shl.form == "sph" && shl.ang > 1 { -(shl.ang as isize) } else { shl.ang as isize }
}

impl Writer for FchWriter {
    fn read_mole(&mut self, mole: &crate::base::Mole) {
        self.atoms = mole.atoms().clone();
        self.basis = mole.basis().clone();
        self.coefs = mole.coefs().clone();
    }
}

fn u2i(nums: &mut Vec<usize>) {
    for num in nums {
        *num as isize;
    }
}
