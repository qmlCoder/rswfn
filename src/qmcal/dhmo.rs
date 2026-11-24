// DHMO方法，大部分工作是徐菁涵做的
use crate::base::{Mole, Vec3};
use crate::writer::Writer;
use crate::{maths, reader, utils, writer};
use ndarray::{Array2, s};
/// 计算sp2杂化原子的有效核电荷数
///
/// 使用经验参数计算原子的有效核电荷数,目前支持B、C、N、O、F、S、Cl元素。
/// 其余元素直接返回核电荷数。
pub fn effective_charge(mole: &Mole, atm: usize) -> f64 {
    let atom = mole.atoms().get(atm);
    let sym = atom.sym.as_str();
    let nebs = mole.atoms().nebs(atm);
    match (sym, nebs.len()) {
        ("B", _) => 4.044,  // 通过苯硼酸确定，给0个电子
        ("C", _) => 4.085,  // 通过联苯确定，给1个电子
        ("N", 2) => 4.153,  // 说明这个N有1个电子参与共轭，通过甲亚胺确定他的有效核电荷
        ("N", 3) => 4.259,  // 说明这个N有2个电子参与共轭，通过苯胺来确定他的有效核电荷
        ("O", 1) => 4.213,  // 说明这个O有1个电子参与共轭，通过甲醛确定他的有效核电荷
        ("O", 2) => 4.306,  // 说明这个O有2个电子参与共轭，通过苯酚确定他的有效核电荷
        ("F", _) => 4.396,  // 通过氟乙烯确定, 给2个电子
        ("S", _) => 6.556,  // 通过噻吩确定
        ("Cl", _) => 6.791, // 通过氯乙烯确定, 给2个电子
        _ => atom.crg,      // 没测，直接就是核电荷
    }
}

/// 计算两个原子间主量子数最高的Slater型p轨道的重叠积分
///
/// 目前可计算:
/// 1. 2pπ-2pπ (同核和异核)
/// 2. 3pπ-2pπ
/// 其他情况直接返回经验值0.14450
///
/// # 参数
/// - `mole`: 分子结构
/// - `ai`: 第一个原子的索引
/// - `aj`: 第二个原子的索引
pub fn overlap_integral(mole: &Mole, ai: usize, aj: usize) -> f64 {
    let atom_i = mole.atoms().get(ai);
    let atom_j = mole.atoms().get(aj);
    let row_i = atom_i.irow();
    let row_j = atom_j.irow();
    let miu_i = effective_charge(mole, ai) / row_i as f64;
    let miu_j = effective_charge(mole, aj) / row_j as f64;
    let dist = mole.atoms().dist(ai, aj);
    let p = (miu_i * miu_i) * dist / 2.0;
    let t = f64::abs((miu_i - miu_j) / (miu_i + miu_j));

    let A = |k| -> f64 {
        let mut sum_A = 0.0;
        for miu in 1..k + 2 {
            sum_A += maths::fac(k) as f64 / (p.powi(miu as i32) * maths::fac(k - miu + 1) as f64);
        }
        f64::exp(-p) * sum_A
    };

    let B = |k: isize| -> f64 {
        let mut sum_B1 = 0.0;
        let mut sum_B2 = 0.0;
        for miu in 1..k + 2 {
            sum_B1 += maths::fac(k) as f64 / ((p * t).powi(miu as i32) * maths::fac(k - miu + 1) as f64);
            sum_B2 += -1.0_f64.powi((k - miu) as i32) * maths::fac(k) as f64
                / ((p * t).powi(miu as i32) * maths::fac(k - miu + 1) as f64);
        }
        -f64::exp(-p * t) * sum_B1 - f64::exp(p * t) * sum_B2
    };
    if row_i == 2 && row_j == 2 {
        if t != 0.0 {
            return p.powi(5) / 32.0
                * (1.0 - t.powi(2)).powf(2.5)
                * (A(4) * (B(0) - B(2)) + A(2) * (B(4) - B(0)) + A(0) * (B(2) - B(4)));
        } else {
            return p.powi(5) / 120.0 * (5.0 * A(4) - 6.0 * A(2) + A(0));
        }
    } else if (row_i == 2 && row_j == 3) || (row_i == 3 && row_j == 2) {
        return p.powi(6) / (32.0 * f64::sqrt(30.0)) * (1.0 + t).powf(2.5) * (1.0 - t).powf(3.5) * // 换行
        (A(5) * (B(0) - B(2)) + A(4) * (B(3) - B(1)) + A(3) * (B(4) - B(0)) + A(2) * (B(1) - B(5)) + A(1) * (B(2) - B(4)) + A(0) * (B(5) - B(3)));
    } else {
        return 0.14450;
    }
}

/// 计算原子的自由价
///
/// 根据原子类型、成键环境以及当前键级总和计算自由价。
/// 自由价 = 最大成键能力 - 配位数 - 总键级
pub fn free_valence(mole: &Mole, atm: usize, sum_order: f64) -> f64 {
    let atoms = mole.atoms();
    let atom = atoms.get(atm);
    let sym = &atom.sym;
    let nebs = atoms.nebs(atm);

    let mut nnbs = 0; //非sp3杂化的重原子数量
    for idx in nebs.iter() {
        if atoms.get(*idx).sym == "H" && atoms.nebs(*idx).len() < 4 {
            nnbs += 1;
        }
    }
    // let mut maxfv = 0.0;
    let maxfv = match (sym.as_str(), nnbs, nebs.len()) {
        ("B", 3, _) => 4.7172, // B(CH2)3 -
        ("C", 1, _) => 4.0000, // 乙烯
        ("C", 2, _) => 4.3333, // 苯
        ("C", 3, _) => 4.7321, // C(CH2)3
        ("N", 1, 2) => 2.9784, // 甲亚胺 提供1个电子的N
        ("N", 2, 2) => 3.3172, // 吡啶
        ("N", 1, 3) => 3.4906, // 乙酰胺 提供2个电子的N
        ("N", 2, 3) => 3.9980, // 吡咯
        ("O", _, 1) => 1.9443, // 甲醛   提供1个电子的O
        ("O", 1, 2) => 2.3372, // 苯酚
        ("O", 2, 2) => 2.9090, // 呋喃
        ("S", 2, _) => 2.8788, // 噻吩
        _ => 0.0,
    };
    maxfv - nebs.len() as f64 - sum_order
}

/// HMO分析,计算键级和自由价
///
/// 基于Hückel分子轨道理论,根据键连矩阵、系数矩阵和电子数计算:
/// - 原子间的键级矩阵
/// - 每个原子的自由价
/// - 键级的标准差
///
/// # 参数
/// - `mole`: 分子结构
/// - `BM`: 原子轨道邻接矩阵
/// - `CM`: 轨道系数矩阵
/// - `ato_atms`: 原子轨道对应的原子索引列表
/// - `nele`: 总电子数
pub fn HMO_Analysis(
    mole: &Mole,
    BM: &Array2<f64>,
    CM: &Array2<f64>,
    ato_atms: &Vec<usize>,
    nele: f64,
) -> (Array2<f64>, Vec<f64>, f64) {
    let atoms = mole.atoms();
    let natm = atoms.len();
    let nato = ato_atms.len();
    let mut omat = Array2::<f64>::zeros((natm, natm));
    let mut fvs = vec![0.0; natm];
    let nobt = CM.ncols();
    let mut obt_occ = vec![]; // 每个分子轨道的电子占据数
    let mut nele = nele;
    for _ in 0..nobt {
        if nele - 2.0 > 0.0 {
            obt_occ.push(2.0);
        } else {
            if nele > 0.0 {
                obt_occ.push(nele);
            } else {
                obt_occ.push(0.0);
            }
        }
        nele -= 2.0;
    }
    // let occ_num = (nele / 2.0).ceil() as usize;
    // log::info!("占据分子轨道数量:{}", occ_num);
    // let CMp = CM.slice(s![.., ..occ_num]).to_owned();
    // let mut orders = vec![];
    for atoi in 0..nato {
        let atmi = ato_atms[atoi];
        let mut sum_order = 0.0;
        for atoj in 0..nato {
            let atmj = ato_atms[atoj];
            if BM[(atoi, atoj)] == 0.0 {
                continue;
            }
            // let order = CMp.row(atoi).dot(&CMp.row(atoj)) * 2.0;
            let mut order = 0.0;
            for obt in 0..nobt {
                if obt_occ[obt] == 0.0 {
                    break;
                }
                order += CM[(atoi, obt)] * CM[(atoj, obt)] * obt_occ[obt];
            }
            omat[(atmi, atmj)] += order;
            sum_order += order;
        }
        fvs[atmi] = free_valence(mole, atmi, sum_order); // 每个原子的自由价
    }
    let mut heavy_order = vec![]; // 不包含H的键的键级
    for ai in 0..natm - 1 {
        for aj in ai + 1..natm {
            if !mole.atoms().is_bond(ai, aj) {
                continue;
            }
            if atoms.get(ai).sym != "H" || atoms.get(aj).sym != "H" {
                heavy_order.push(omat[(ai, aj)]);
            }
        }
    }
    println!("计算标准差:{:?}", heavy_order);
    let sd = maths::vec_psd(&heavy_order);
    (omat, fvs, sd)
}

#[derive(Clone)]
/// DHMO的计算结果
pub struct DmoRes {
    /// 原子轨道列表,每个元素包含(原子索引,壳层索引,轨道符号,轨道方向向量)
    pub atos: Vec<(usize, usize, String, [f64; 3])>,
    /// 键连矩阵
    pub BM: Array2<f64>,
    /// 轨道系数矩阵
    pub CM: Array2<f64>,
    /// 分子轨道的能量
    pub es: Vec<f64>,
}

/// 执行dHMO计算
///
/// 读取.dmi输入文件,构建哈密顿矩阵并求解广义本征值问题,
/// 得到轨道系数和能量,返回分子结构和dHMO计算结果。
pub fn dhmo_calc(path: String) -> (Box<Mole>, DmoRes) {
    // 计算原子的局部坐标系

    let dmi_reader = reader::DmoReader::new(&path);
    let mole = Mole::from_read(&dmi_reader);
    let atoms = mole.atoms();
    let atos = dmi_reader.read_atos().clone();
    let nato = atos.len(); // 原子轨道数量
    let mut BM = Array2::<f64>::zeros((nato, nato)); // 键连矩阵
    let s135 = 0.14450;
    for i in 0..nato {
        let ai = atos[i].0;
        let irow = mole.atoms().get(ai).irow() as f64; // 该元素在元素周期表中的行
        let diri = Vec3::from_arr(atos[i].3);
        BM[(i, i)] = effective_charge(&mole, ai).powi(2) / (2.0 * irow.powi(2));
        for j in i..nato {
            let aj = atos[j].0;
            if ai >= aj {
                continue;
            }
            if !atoms.is_bond(ai, aj) {
                continue;
            }
            let dirj = Vec3::from_arr(atos[j].3);

            // 角度校正因子
            let bond_vec: Vec3 = atoms.vec(ai, aj).normalize();
            let bond_len = bond_vec.length();
            let cos_ai = f64::sqrt(1.0 - diri.dot(&bond_vec) / diri.length() * bond_len).powi(2);
            let cos_aj = f64::sqrt(1.0 - dirj.dot(&bond_vec) / dirj.length() * bond_len).powi(2);
            let cross_i = diri.cross(&bond_vec);
            let cross_j = dirj.cross(&bond_vec);
            let cos_aij = cross_i.dot(&cross_j) / (cross_i.length() * cross_j.length());
            if cos_ai.is_nan() || cos_aj.is_nan() || cos_aij.is_nan() {
                log::error!("角度矫正包含NAN,{:>10.4}{:>10.4}{:>10.4}", cos_ai, cos_aj, cos_aij);
                panic!("角度矫正包含NAN");
            }
            let fa = cos_ai * cos_aj * cos_aij;
            // 距离校正因子
            let sij = overlap_integral(&mole, ai, aj);
            let fd = sij / s135;
            // 元素矫正因子
            let fe =
                (effective_charge(&mole, ai).powi(2) + effective_charge(&mole, aj).powi(2)) / (2.0 * 4.085_f64.powi(2));
            if fa.is_nan() || fd.is_nan() || fe.is_nan() {
                log::error!("校正系数包含NaN,fa={:>10.4},fd={:>10.4},fe={:>10.4}", fa, fd, fe);
                panic!("校正系数包含NaN");
            }
            BM[(i, j)] = fa * fd * fe * 0.13906;
            BM[(j, i)] = fa * fd * fe * 0.13906;
        }
    }
    let (engs, cmat) = maths::eigh(&BM);

    let mut engs = engs.flatten().to_vec();
    let cmat = cmat.slice_move(s![..,..;-1]);
    engs.reverse();
    let res = DmoRes {
        atos,
        BM,
        CM: cmat,
        es: engs,
    };
    (Box::new(mole), res)
}

/// 计算并打印结果
pub fn output(dmi_path: &str, dmo_path: Option<String>) {
    let dmo_path = match dmo_path {
        Some(dmo_path) => dmo_path.to_string(),
        None => dmi_path.replace(".dmi", ".dmo"),
    };
    let (mole, dmo_res) = dhmo_calc(dmi_path.to_string());
    let mut writer = writer::DmoWriter::blank();
    writer.read_mole(&mole);
    writer.atos = dmo_res.atos.clone();
    writer.BM = dmo_res.BM.clone();
    writer.CM = dmo_res.CM.clone();
    writer.es = dmo_res.es.clone();
    writer.save(&dmo_path);
}
