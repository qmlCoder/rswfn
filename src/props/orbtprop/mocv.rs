#![allow(non_snake_case)]
use crate::base::Mole;
use crate::datas::FACS;
use crate::datas::btrans;
use crate::maths;
use crate::utils;
/// 局部坐标系与全局坐标系的变换矩阵
/// C^l=(T^T)^{-1} @ C^g
/// C^g=T^T @ C^l
use ndarray::Array2;
use std::collections::HashMap;

/// 使用 MOCV 方法将全局坐标系系数映射到局部坐标系再映射回来
pub fn mocv(L: &Array2<f64>, rcoefs: &Array2<f64>, keeps: &Array2<f64>) -> Array2<f64> {
    let nrow = rcoefs.nrows();
    let Mo = if [1, 3, 5, 7, 9, 11].contains(&nrow) {
        get_tmat_sph(nrow, L)
    } else {
        // mocv::get_tmat_car(nrow, L)
        panic!("只能使用球谐轨道系数")
    };
    // let Mi = mocv::get_tmat_sph(nrow, L);
    let Mt = Mo.t().to_owned();
    let Mi = maths::inv(&Mt);
    let mut tcoefs = Mi.dot(rcoefs); //全局系数映射到局部坐标系
    tcoefs = &tcoefs * keeps; // 保留某个维度
    let fcoefs = Mt.dot(&tcoefs); // 再映射回全局坐标系
    fcoefs
}

/// dest：计算目标，可以为 0:popul 或 1:order
/// 根据计算目标选择用于 MOCV 的原子轨道索引
pub fn mocv_atos(mole: &Mole, dest: &MocvDest) -> Vec<usize> {
    let mut atos = vec![];
    let ato_syms = mole.basis().atos().syms("sph");
    let ato_atms = mole.basis().atos().atms("sph");

    let mut atm_syms = HashMap::new(); // 每个原子保留的符号
    let natm = mole.atoms.len();
    for i in 0..natm {
        let geom = if mole.atoms.is_line(i) { 1 } else { 0 };
        let syms = mocv_syms(dest, geom);
        atm_syms.insert(i, syms);
    }
    let nato = ato_syms.len();
    let ato_locs = mole.basis().atos().locs("sph");
    for i in 0..nato {
        let atm = ato_atms[i];
        let sym = &ato_syms[i];
        let keeps = atm_syms.get(&atm).unwrap();
        if keeps.contains(sym) {
            atos.push(i);
        }
    }
    atos
}

pub enum OrderDType {
    SI, // σ键级
    PI, // π键级
    DE, // δ键级
    PH, // φ键级
}

pub enum MocvDest {
    PiProperty,             // pi电子性质
    OrderDecom(OrderDType), // 键级分解
    RoundOrder,             // 环轴键级
}

impl MocvDest {
    pub fn from_str(s: &str) -> Option<Self> {
        match s {
            "PiProperty" => Some(MocvDest::PiProperty),
            "OrderDecom-SI" => Some(MocvDest::OrderDecom(OrderDType::SI)),
            "OrderDecom-PI" => Some(MocvDest::OrderDecom(OrderDType::PI)),
            "OrderDecom-DE" => Some(MocvDest::OrderDecom(OrderDType::DE)),
            "OrderDecom-PH" => Some(MocvDest::OrderDecom(OrderDType::PH)),
            "RoundOrder" => Some(MocvDest::RoundOrder),
            _ => None,
        }
    }
}

// 计算目标：广义pi性质-平面/直线、键级分解、环轴键级
/// dest：计算目标，可以为 MocvDest 枚举值
/// geom：几何类型，可以为 0:plan 或 1:line
// #[cached]
/// 根据计算目标和几何类型返回要保留的轨道对称性
pub fn mocv_syms(dest: &MocvDest, geom: usize) -> Vec<String> {
    let syms = match (dest, geom) {
        (MocvDest::PiProperty, 0) => {
            vec![
                "PX", // P轨道
                "D+1", "D-2", // 球谐D轨道
                "F+1", "F-2", "F+3", // 球谐F轨道
                "G+1", "G-2", "G+3", "G-4",
            ] // 平面分子电子布局
        }
        (MocvDest::PiProperty, 1) => {
            vec![
                "PX", "PY", // P轨道
                "D+1", "D-1", "D+2", "D-2", // 球谐D轨道
                "F+1", "F-1", "F+2", "F-2", "F+3", "F-3", // 球谐F轨道
                "G+1", "G-1", "G+2", "G-2", "G+3", "G-3", "G+4", "G-4", // 球谐G轨道
            ] // 线性分子电子布居
        }
        (MocvDest::OrderDecom(OrderDType::SI), _) => {
            vec!["S", "PZ", "D 0", "F 0", "G 0"]
        }
        (MocvDest::OrderDecom(OrderDType::PI), _) => {
            vec!["PX", "PY", "D+1", "D-1", "F+1", "F-1", "G+1", "G-1"]
        }
        (MocvDest::OrderDecom(OrderDType::DE), _) => {
            vec!["D+2", "D-2", "F+2", "F-2", "G+2", "G-2"]
        }
        (MocvDest::OrderDecom(OrderDType::PH), _) => {
            vec!["F+3", "F-3", "G+3", "G-3"]
        }
        (MocvDest::RoundOrder, _) => {
            vec![
                "PX", // P轨道
                "D+1", "D-2", // 球谐D轨道
                "F+1", "F-2", "F+3", // 球谐F轨道
                "G+1", "G-2", "G+3", "G-4",
            ] // 环轴键级
        }
        _ => vec![],
    };
    syms.into_iter().map(|sym| sym.to_string()).collect()
}

/// 根据壳层维度自动选择合适的变换矩阵
pub fn get_tmat(num: usize, lmat: &Array2<f64>) -> Array2<f64> {
    if [1, 3, 5, 7, 9].contains(&num) {
        return get_tmat_sph(num, lmat);
    }
    if [1, 3, 6, 10, 15].contains(&num) {
        return get_tmat_car(num, lmat);
    }
    panic!("不支持的类型:{}", num)
}

/// 获取球谐函数形式的变换矩阵
///
/// 球谐基与笛卡尔基通过 U (DMAT/FMAT/GMAT) 相关联：φ = χ · U，且 U 的列相对
/// 笛卡尔角度重叠 S 正交归一（Uᵀ S U = I），而非相对欧氏内积正交（Uᵀ U ≠ I）。
/// 设笛卡尔旋转矩阵 R_car = get_tmat_car(..)ᵀ（满足 C^g = R_car C^l），则
/// 球谐基下的旋转矩阵为
///
///   R_sph = Uᵀ S R_car U
///
/// 对应本函数的返回（T = R_sphᵀ，满足 C^g = Tᵀ C^l）：
///
///   T_sph = R_sphᵀ = Uᵀ · get_tmat_car(..) · S · U
///
/// 当无旋转（R_car = I）时，T_sph = Uᵀ S U = I，与局部系==全局系的物理意义一致。
pub fn get_tmat_sph(num: usize, lmat: &Array2<f64>) -> Array2<f64> {
    match num {
        01 => tmat_s(lmat),
        03 => tmat_p(lmat),
        05 => btrans::DMAT
            .t()
            .dot(&get_tmat_car(06, lmat))
            .dot(&sang_d())
            .dot(&*btrans::DMAT),
        07 => btrans::FMAT
            .t()
            .dot(&get_tmat_car(10, lmat))
            .dot(&sang_f())
            .dot(&*btrans::FMAT),
        09 => btrans::GMAT
            .t()
            .dot(&get_tmat_car(15, lmat))
            .dot(&sang_g())
            .dot(&*btrans::GMAT),
        _ => panic!("不支持的类型:{}", num),
    }
}

/// 获取变换矩阵,输入的矩阵是列为基向量，要转为行为基向量
pub fn get_tmat_car(num: usize, lmat: &Array2<f64>) -> Array2<f64> {
    match num {
        01 => tmat_s(lmat),
        03 => tmat_p(lmat),
        06 => tmat_d(lmat) * nmat_d(),
        10 => tmat_f(lmat) * nmat_f(),
        15 => tmat_g(lmat) * nmat_g(),
        _ => panic!("不支持的类型:{}", num),
    }
}

fn get_ts(L: &Array2<f64>) -> [f64; 9] {
    let L11 = L[(0, 0)];
    let L12 = L[(0, 1)];
    let L13 = L[(0, 2)];
    let L21 = L[(1, 0)];
    let L22 = L[(1, 1)];
    let L23 = L[(1, 2)];
    let L31 = L[(2, 0)];
    let L32 = L[(2, 1)];
    let L33 = L[(2, 2)];
    return [L11, L12, L13, L21, L22, L23, L31, L32, L33];
}

/// 获取 s 轨道的变换矩阵（始终返回1x1单位矩阵）
pub fn tmat_s(L: &Array2<f64>) -> Array2<f64> {
    Array2::ones((1, 1))
}

/// 获取 p 轨道的变换矩阵
pub fn tmat_p(L: &Array2<f64>) -> Array2<f64> {
    let mut MP: Array2<f64> = Array2::zeros((3, 3));
    let [L11, L12, L13, L21, L22, L23, L31, L32, L33] = get_ts(L);
    MP[(0, 0)] = L11;
    MP[(0, 1)] = L21;
    MP[(0, 2)] = L31;
    MP[(1, 0)] = L12;
    MP[(1, 1)] = L22;
    MP[(1, 2)] = L32;
    MP[(2, 0)] = L13;
    MP[(2, 1)] = L23;
    MP[(2, 2)] = L33;
    MP
}

/// 获取 d 轨道的变换矩阵
pub fn tmat_d(L: &Array2<f64>) -> Array2<f64> {
    let [L11, L12, L13, L21, L22, L23, L31, L32, L33] = get_ts(L);
    let mut MD: Array2<f64> = Array2::zeros((6, 6));
    MD[(0, 0)] = L11.powi(2);
    MD[(0, 1)] = L21.powi(2);
    MD[(0, 2)] = L31.powi(2);
    MD[(0, 3)] = 2.0 * L11 * L21;
    MD[(0, 4)] = 2.0 * L11 * L31;
    MD[(0, 5)] = 2.0 * L21 * L31;
    MD[(1, 0)] = L12.powi(2);
    MD[(1, 1)] = L22.powi(2);
    MD[(1, 2)] = L32.powi(2);
    MD[(1, 3)] = 2.0 * L12 * L22;
    MD[(1, 4)] = 2.0 * L12 * L32;
    MD[(1, 5)] = 2.0 * L22 * L32;
    MD[(2, 0)] = L13.powi(2);
    MD[(2, 1)] = L23.powi(2);
    MD[(2, 2)] = L33.powi(2);
    MD[(2, 3)] = 2.0 * L13 * L23;
    MD[(2, 4)] = 2.0 * L13 * L33;
    MD[(2, 5)] = 2.0 * L23 * L33;
    MD[(3, 0)] = L11 * L12;
    MD[(3, 1)] = L21 * L22;
    MD[(3, 2)] = L31 * L32;
    MD[(3, 3)] = L11 * L22 + L12 * L21;
    MD[(3, 4)] = L11 * L32 + L12 * L31;
    MD[(3, 5)] = L21 * L32 + L22 * L31;
    MD[(4, 0)] = L11 * L13;
    MD[(4, 1)] = L21 * L23;
    MD[(4, 2)] = L31 * L33;
    MD[(4, 3)] = L11 * L23 + L13 * L21;
    MD[(4, 4)] = L11 * L33 + L13 * L31;
    MD[(4, 5)] = L21 * L33 + L23 * L31;
    MD[(5, 0)] = L12 * L13;
    MD[(5, 1)] = L22 * L23;
    MD[(5, 2)] = L32 * L33;
    MD[(5, 3)] = L12 * L23 + L13 * L22;
    MD[(5, 4)] = L12 * L33 + L13 * L32;
    MD[(5, 5)] = L22 * L33 + L23 * L32;
    MD
}

// #[rustfmt::skip]
/// 获取 f 轨道的变换矩阵
pub fn tmat_f(L: &Array2<f64>) -> Array2<f64> {
    let [L11, L12, L13, L21, L22, L23, L31, L32, L33] = get_ts(L);
    let mut MF: Array2<f64> = Array2::zeros((10, 10));
    MF[(0, 0)] = L11.powi(3);
    MF[(0, 1)] = L21.powi(3);
    MF[(0, 2)] = L31.powi(3);
    MF[(0, 3)] = 3.0 * L11 * L21.powi(2);
    MF[(0, 4)] = 3.0 * L11.powi(2) * L21;
    MF[(0, 5)] = 3.0 * L11.powi(2) * L31;
    MF[(0, 6)] = 3.0 * L11 * L31.powi(2);
    MF[(0, 7)] = 3.0 * L21 * L31.powi(2);
    MF[(0, 8)] = 3.0 * L21.powi(2) * L31;
    MF[(0, 9)] = 6.0 * L11 * L21 * L31;
    MF[(1, 0)] = L12.powi(3);
    MF[(1, 1)] = L22.powi(3);
    MF[(1, 2)] = L32.powi(3);
    MF[(1, 3)] = 3.0 * L12 * L22.powi(2);
    MF[(1, 4)] = 3.0 * L12.powi(2) * L22;
    MF[(1, 5)] = 3.0 * L12.powi(2) * L32;
    MF[(1, 6)] = 3.0 * L12 * L32.powi(2);
    MF[(1, 7)] = 3.0 * L22 * L32.powi(2);
    MF[(1, 8)] = 3.0 * L22.powi(2) * L32;
    MF[(1, 9)] = 6.0 * L12 * L22 * L32;
    MF[(2, 0)] = L13.powi(3);
    MF[(2, 1)] = L23.powi(3);
    MF[(2, 2)] = L33.powi(3);
    MF[(2, 3)] = 3.0 * L13 * L23.powi(2);
    MF[(2, 4)] = 3.0 * L13.powi(2) * L23;
    MF[(2, 5)] = 3.0 * L13.powi(2) * L33;
    MF[(2, 6)] = 3.0 * L13 * L33.powi(2);
    MF[(2, 7)] = 3.0 * L23 * L33.powi(2);
    MF[(2, 8)] = 3.0 * L23.powi(2) * L33;
    MF[(2, 9)] = 6.0 * L13 * L23 * L33;
    MF[(3, 0)] = L11 * L12.powi(2);
    MF[(3, 1)] = L21 * L22.powi(2);
    MF[(3, 2)] = L31 * L32.powi(2);
    MF[(3, 3)] = L11 * L22.powi(2) + 2.0 * L12 * L21 * L22;
    MF[(3, 4)] = 2.0 * L11 * L12 * L22 + L12.powi(2) * L21;
    MF[(3, 5)] = 2.0 * L11 * L12 * L32 + L12.powi(2) * L31;
    MF[(3, 6)] = L11 * L32.powi(2) + 2.0 * L12 * L31 * L32;
    MF[(3, 7)] = L21 * L32.powi(2) + 2.0 * L22 * L31 * L32;
    MF[(3, 8)] = 2.0 * L21 * L22 * L32 + L22.powi(2) * L31;
    MF[(3, 9)] = 2.0 * L11 * L22 * L32 + 2.0 * L12 * L21 * L32 + 2.0 * L12 * L22 * L31;
    MF[(4, 0)] = L11.powi(2) * L12;
    MF[(4, 1)] = L21.powi(2) * L22;
    MF[(4, 2)] = L31.powi(2) * L32;
    MF[(4, 3)] = 2.0 * L11 * L21 * L22 + L12 * L21.powi(2);
    MF[(4, 4)] = L11.powi(2) * L22 + 2.0 * L11 * L12 * L21;
    MF[(4, 5)] = L11.powi(2) * L32 + 2.0 * L11 * L12 * L31;
    MF[(4, 6)] = 2.0 * L11 * L31 * L32 + L12 * L31.powi(2);
    MF[(4, 7)] = 2.0 * L21 * L31 * L32 + L22 * L31.powi(2);
    MF[(4, 8)] = L21.powi(2) * L32 + 2.0 * L21 * L22 * L31;
    MF[(4, 9)] = 2.0 * L11 * L21 * L32 + 2.0 * L11 * L22 * L31 + 2.0 * L12 * L21 * L31;
    MF[(5, 0)] = L11.powi(2) * L13;
    MF[(5, 1)] = L21.powi(2) * L23;
    MF[(5, 2)] = L31.powi(2) * L33;
    MF[(5, 3)] = 2.0 * L11 * L21 * L23 + L13 * L21.powi(2);
    MF[(5, 4)] = L11.powi(2) * L23 + 2.0 * L11 * L13 * L21;
    MF[(5, 5)] = L11.powi(2) * L33 + 2.0 * L11 * L13 * L31;
    MF[(5, 6)] = 2.0 * L11 * L31 * L33 + L13 * L31.powi(2);
    MF[(5, 7)] = 2.0 * L21 * L31 * L33 + L23 * L31.powi(2);
    MF[(5, 8)] = L21.powi(2) * L33 + 2.0 * L21 * L23 * L31;
    MF[(5, 9)] = 2.0 * L11 * L21 * L33 + 2.0 * L11 * L23 * L31 + 2.0 * L13 * L21 * L31;
    MF[(6, 0)] = L11 * L13.powi(2);
    MF[(6, 1)] = L21 * L23.powi(2);
    MF[(6, 2)] = L31 * L33.powi(2);
    MF[(6, 3)] = L11 * L23.powi(2) + 2.0 * L13 * L21 * L23;
    MF[(6, 4)] = 2.0 * L11 * L13 * L23 + L13.powi(2) * L21;
    MF[(6, 5)] = 2.0 * L11 * L13 * L33 + L13.powi(2) * L31;
    MF[(6, 6)] = L11 * L33.powi(2) + 2.0 * L13 * L31 * L33;
    MF[(6, 7)] = L21 * L33.powi(2) + 2.0 * L23 * L31 * L33;
    MF[(6, 8)] = 2.0 * L21 * L23 * L33 + L23.powi(2) * L31;
    MF[(6, 9)] = 2.0 * L11 * L23 * L33 + 2.0 * L13 * L21 * L33 + 2.0 * L13 * L23 * L31;
    MF[(7, 0)] = L12 * L13.powi(2);
    MF[(7, 1)] = L22 * L23.powi(2);
    MF[(7, 2)] = L32 * L33.powi(2);
    MF[(7, 3)] = L12 * L23.powi(2) + 2.0 * L13 * L22 * L23;
    MF[(7, 4)] = 2.0 * L12 * L13 * L23 + L13.powi(2) * L22;
    MF[(7, 5)] = 2.0 * L12 * L13 * L33 + L13.powi(2) * L32;
    MF[(7, 6)] = L12 * L33.powi(2) + 2.0 * L13 * L32 * L33;
    MF[(7, 7)] = L22 * L33.powi(2) + 2.0 * L23 * L32 * L33;
    MF[(7, 8)] = 2.0 * L22 * L23 * L33 + L23.powi(2) * L32;
    MF[(7, 9)] = 2.0 * L12 * L23 * L33 + 2.0 * L13 * L22 * L33 + 2.0 * L13 * L23 * L32;
    MF[(8, 0)] = L12.powi(2) * L13;
    MF[(8, 1)] = L22.powi(2) * L23;
    MF[(8, 2)] = L32.powi(2) * L33;
    MF[(8, 3)] = 2.0 * L12 * L22 * L23 + L13 * L22.powi(2);
    MF[(8, 4)] = L12.powi(2) * L23 + 2.0 * L12 * L13 * L22;
    MF[(8, 5)] = L12.powi(2) * L33 + 2.0 * L12 * L13 * L32;
    MF[(8, 6)] = 2.0 * L12 * L32 * L33 + L13 * L32.powi(2);
    MF[(8, 7)] = 2.0 * L22 * L32 * L33 + L23 * L32.powi(2);
    MF[(8, 8)] = L22.powi(2) * L33 + 2.0 * L22 * L23 * L32;
    MF[(8, 9)] = 2.0 * L12 * L22 * L33 + 2.0 * L12 * L23 * L32 + 2.0 * L13 * L22 * L32;
    MF[(9, 0)] = L11 * L12 * L13;
    MF[(9, 1)] = L21 * L22 * L23;
    MF[(9, 2)] = L31 * L32 * L33;
    MF[(9, 3)] = L11 * L22 * L23 + L12 * L21 * L23 + L13 * L21 * L22;
    MF[(9, 4)] = L11 * L12 * L23 + L11 * L13 * L22 + L12 * L13 * L21;
    MF[(9, 5)] = L11 * L12 * L33 + L11 * L13 * L32 + L12 * L13 * L31;
    MF[(9, 6)] = L11 * L32 * L33 + L12 * L31 * L33 + L13 * L31 * L32;
    MF[(9, 7)] = L21 * L32 * L33 + L22 * L31 * L33 + L23 * L31 * L32;
    MF[(9, 8)] = L21 * L22 * L33 + L21 * L23 * L32 + L22 * L23 * L31;
    MF[(9, 9)] =
        L11 * L22 * L33 + L11 * L23 * L32 + L12 * L21 * L33 + L12 * L23 * L31 + L13 * L21 * L32 + L13 * L22 * L31;
    MF
}

#[rustfmt::skip]
/// 获取 g 轨道的变换矩阵
pub fn tmat_g(L: &Array2<f64>) -> Array2<f64> {
    let [L11, L12, L13, L21, L22, L23, L31, L32, L33] = get_ts(L);
    let mut MG: Array2<f64> = Array2::zeros((15, 15));
    MG[(0,0)]=L33.powi(4);
    MG[(0,1)]=4.0*L23*L33.powi(3);
    MG[(0,2)]=6.0*L23.powi(2)*L33.powi(2);
    MG[(0,3)]=4.0*L23.powi(3)*L33;
    MG[(0,4)]=L23.powi(4);
    MG[(0,5)]=4.0*L13*L33.powi(3);
    MG[(0,6)]=12.0*L13*L23*L33.powi(2);
    MG[(0,7)]=12.0*L13*L23.powi(2)*L33;
    MG[(0,8)]=4.0*L13*L23.powi(3);
    MG[(0,9)]=6.0*L13.powi(2)*L33.powi(2);
    MG[(0,10)]=12.0*L13.powi(2)*L23*L33;
    MG[(0,11)]=6.0*L13.powi(2)*L23.powi(2);
    MG[(0,12)]=4.0*L13.powi(3)*L33;
    MG[(0,13)]=4.0*L13.powi(3)*L23;
    MG[(0,14)]=L13.powi(4);
    MG[(1,0)]=L32*L33.powi(3);
    MG[(1,1)]=L22*L33.powi(3) + 3.0*L23*L32*L33.powi(2);
    MG[(1,2)]=3.0*L22*L23*L33.powi(2) + 3.0*L23.powi(2)*L32*L33;
    MG[(1,3)]=3.0*L22*L23.powi(2)*L33 + L23.powi(3)*L32;
    MG[(1,4)]=L22*L23.powi(3);
    MG[(1,5)]=L12*L33.powi(3) + 3.0*L13*L32*L33.powi(2);
    MG[(1,6)]=3.0*L12*L23*L33.powi(2) + 3.0*L13*L22*L33.powi(2) + 6.0*L13*L23*L32*L33;
    MG[(1,7)]=3.0*L12*L23.powi(2)*L33 + 6.0*L13*L22*L23*L33 + 3.0*L13*L23.powi(2)*L32;
    MG[(1,8)]=L12*L23.powi(3) + 3.0*L13*L22*L23.powi(2);
    MG[(1,9)]=3.0*L12*L13*L33.powi(2) + 3.0*L13.powi(2)*L32*L33;
    MG[(1,10)]=6.0*L12*L13*L23*L33 + 3.0*L13.powi(2)*L22*L33 + 3.0*L13.powi(2)*L23*L32;
    MG[(1,11)]=3.0*L12*L13*L23.powi(2) + 3.0*L13.powi(2)*L22*L23;
    MG[(1,12)]=3.0*L12*L13.powi(2)*L33 + L13.powi(3)*L32;
    MG[(1,13)]=3.0*L12*L13.powi(2)*L23 + L13.powi(3)*L22;
    MG[(1,14)]=L12*L13.powi(3);
    MG[(2,0)]=L32.powi(2)*L33.powi(2);
    MG[(2,1)]=2.0*L22*L32*L33.powi(2) + 2.0*L23*L32.powi(2)*L33;
    MG[(2,2)]=L22.powi(2)*L33.powi(2) + 4.0*L22*L23*L32*L33 + L23.powi(2)*L32.powi(2);
    MG[(2,3)]=2.0*L22.powi(2)*L23*L33 + 2.0*L22*L23.powi(2)*L32;
    MG[(2,4)]=L22.powi(2)*L23.powi(2);
    MG[(2,5)]=2.0*L12*L32*L33.powi(2) + 2.0*L13*L32.powi(2)*L33;
    MG[(2,6)]=2.0*L12*L22*L33.powi(2) + 4.0*L12*L23*L32*L33 + 4.0*L13*L22*L32*L33 + 2.0*L13*L23*L32.powi(2);
    MG[(2,7)]=4.0*L12*L22*L23*L33 + 2.0*L12*L23.powi(2)*L32 + 2.0*L13*L22.powi(2)*L33 + 4.0*L13*L22*L23*L32;
    MG[(2,8)]=2.0*L12*L22*L23.powi(2) + 2.0*L13*L22.powi(2)*L23;
    MG[(2,9)]=L12.powi(2)*L33.powi(2) + 4.0*L12*L13*L32*L33 + L13.powi(2)*L32.powi(2);
    MG[(2,10)]=2.0*L12.powi(2)*L23*L33 + 4.0*L12*L13*L22*L33 + 4.0*L12*L13*L23*L32 + 2.0*L13.powi(2)*L22*L32;
    MG[(2,11)]=L12.powi(2)*L23.powi(2) + 4.0*L12*L13*L22*L23 + L13.powi(2)*L22.powi(2);
    MG[(2,12)]=2.0*L12.powi(2)*L13*L33 + 2.0*L12*L13.powi(2)*L32;
    MG[(2,13)]=2.0*L12.powi(2)*L13*L23 + 2.0*L12*L13.powi(2)*L22;
    MG[(2,14)]=L12.powi(2)*L13.powi(2);
    MG[(3,0)]=L32.powi(3)*L33;
    MG[(3,1)]=3.0*L22*L32.powi(2)*L33 + L23*L32.powi(3);
    MG[(3,2)]=3.0*L22.powi(2)*L32*L33 + 3.0*L22*L23*L32.powi(2);
    MG[(3,3)]=L22.powi(3)*L33 + 3.0*L22.powi(2)*L23*L32;
    MG[(3,4)]=L22.powi(3)*L23;
    MG[(3,5)]=3.0*L12*L32.powi(2)*L33 + L13*L32.powi(3);
    MG[(3,6)]=6.0*L12*L22*L32*L33 + 3.0*L12*L23*L32.powi(2) + 3.0*L13*L22*L32.powi(2);
    MG[(3,7)]=3.0*L12*L22.powi(2)*L33 + 6.0*L12*L22*L23*L32 + 3.0*L13*L22.powi(2)*L32;
    MG[(3,8)]=3.0*L12*L22.powi(2)*L23 + L13*L22.powi(3);
    MG[(3,9)]=3.0*L12.powi(2)*L32*L33 + 3.0*L12*L13*L32.powi(2);
    MG[(3,10)]=3.0*L12.powi(2)*L22*L33 + 3.0*L12.powi(2)*L23*L32 + 6.0*L12*L13*L22*L32;
    MG[(3,11)]=3.0*L12.powi(2)*L22*L23 + 3.0*L12*L13*L22.powi(2);
    MG[(3,12)]=L12.powi(3)*L33 + 3.0*L12.powi(2)*L13*L32;
    MG[(3,13)]=L12.powi(3)*L23 + 3.0*L12.powi(2)*L13*L22;
    MG[(3,14)]=L12.powi(3)*L13;
    MG[(4,0)]=L32.powi(4);
    MG[(4,1)]=4.0*L22*L32.powi(3);
    MG[(4,2)]=6.0*L22.powi(2)*L32.powi(2);
    MG[(4,3)]=4.0*L22.powi(3)*L32;
    MG[(4,4)]=L22.powi(4);
    MG[(4,5)]=4.0*L12*L32.powi(3);
    MG[(4,6)]=12.0*L12*L22*L32.powi(2);
    MG[(4,7)]=12.0*L12*L22.powi(2)*L32;
    MG[(4,8)]=4.0*L12*L22.powi(3);
    MG[(4,9)]=6.0*L12.powi(2)*L32.powi(2);
    MG[(4,10)]=12.0*L12.powi(2)*L22*L32;
    MG[(4,11)]=6.0*L12.powi(2)*L22.powi(2);
    MG[(4,12)]=4.0*L12.powi(3)*L32;
    MG[(4,13)]=4.0*L12.powi(3)*L22;
    MG[(4,14)]=L12.powi(4);
    MG[(5,0)]=L31*L33.powi(3);
    MG[(5,1)]=L21*L33.powi(3) + 3.0*L23*L31*L33.powi(2);
    MG[(5,2)]=3.0*L21*L23*L33.powi(2) + 3.0*L23.powi(2)*L31*L33;
    MG[(5,3)]=3.0*L21*L23.powi(2)*L33 + L23.powi(3)*L31;
    MG[(5,4)]=L21*L23.powi(3);
    MG[(5,5)]=L11*L33.powi(3) + 3.0*L13*L31*L33.powi(2);
    MG[(5,6)]=3.0*L11*L23*L33.powi(2) + 3.0*L13*L21*L33.powi(2) + 6.0*L13*L23*L31*L33;
    MG[(5,7)]=3.0*L11*L23.powi(2)*L33 + 6.0*L13*L21*L23*L33 + 3.0*L13*L23.powi(2)*L31;
    MG[(5,8)]=L11*L23.powi(3) + 3.0*L13*L21*L23.powi(2);
    MG[(5,9)]=3.0*L11*L13*L33.powi(2) + 3.0*L13.powi(2)*L31*L33;
    MG[(5,10)]=6.0*L11*L13*L23*L33 + 3.0*L13.powi(2)*L21*L33 + 3.0*L13.powi(2)*L23*L31;
    MG[(5,11)]=3.0*L11*L13*L23.powi(2) + 3.0*L13.powi(2)*L21*L23;
    MG[(5,12)]=3.0*L11*L13.powi(2)*L33 + L13.powi(3)*L31;
    MG[(5,13)]=3.0*L11*L13.powi(2)*L23 + L13.powi(3)*L21;
    MG[(5,14)]=L11*L13.powi(3);
    MG[(6,0)]=L31*L32*L33.powi(2);
    MG[(6,1)]=L21*L32*L33.powi(2) + L22*L31*L33.powi(2) + 2.0*L23*L31*L32*L33;
    MG[(6,2)]=L21*L22*L33.powi(2) + 2.0*L21*L23*L32*L33 + 2.0*L22*L23*L31*L33 + L23.powi(2)*L31*L32;
    MG[(6,3)]=2.0*L21*L22*L23*L33 + L21*L23.powi(2)*L32 + L22*L23.powi(2)*L31;
    MG[(6,4)]=L21*L22*L23.powi(2);
    MG[(6,5)]=L11*L32*L33.powi(2) + L12*L31*L33.powi(2) + 2.0*L13*L31*L32*L33;
    MG[(6,6)]=L11*L22*L33.powi(2) + 2.0*L11*L23*L32*L33 + L12*L21*L33.powi(2) + 2.0*L12*L23*L31*L33 + 2.0*L13*L21*L32*L33 + 2.0*L13*L22*L31*L33 + 2.0*L13*L23*L31*L32;
    MG[(6,7)]=2.0*L11*L22*L23*L33 + L11*L23.powi(2)*L32 + 2.0*L12*L21*L23*L33 + L12*L23.powi(2)*L31 + 2.0*L13*L21*L22*L33 + 2.0*L13*L21*L23*L32 + 2.0*L13*L22*L23*L31;
    MG[(6,8)]=L11*L22*L23.powi(2) + L12*L21*L23.powi(2) + 2.0*L13*L21*L22*L23;
    MG[(6,9)]=L11*L12*L33.powi(2) + 2.0*L11*L13*L32*L33 + 2.0*L12*L13*L31*L33 + L13.powi(2)*L31*L32;
    MG[(6,10)]=2.0*L11*L12*L23*L33 + 2.0*L11*L13*L22*L33 + 2.0*L11*L13*L23*L32 + 2.0*L12*L13*L21*L33 + 2.0*L12*L13*L23*L31 + L13.powi(2)*L21*L32 + L13.powi(2)*L22*L31;
    MG[(6,11)]=L11*L12*L23.powi(2) + 2.0*L11*L13*L22*L23 + 2.0*L12*L13*L21*L23 + L13.powi(2)*L21*L22;
    MG[(6,12)]=2.0*L11*L12*L13*L33 + L11*L13.powi(2)*L32 + L12*L13.powi(2)*L31;
    MG[(6,13)]=2.0*L11*L12*L13*L23 + L11*L13.powi(2)*L22 + L12*L13.powi(2)*L21;
    MG[(6,14)]=L11*L12*L13.powi(2);
    MG[(7,0)]=L31*L32.powi(2)*L33;
    MG[(7,1)]=L21*L32.powi(2)*L33 + 2.0*L22*L31*L32*L33 + L23*L31*L32.powi(2);
    MG[(7,2)]=2.0*L21*L22*L32*L33 + L21*L23*L32.powi(2) + L22.powi(2)*L31*L33 + 2.0*L22*L23*L31*L32;
    MG[(7,3)]=L21*L22.powi(2)*L33 + 2.0*L21*L22*L23*L32 + L22.powi(2)*L23*L31;
    MG[(7,4)]=L21*L22.powi(2)*L23;
    MG[(7,5)]=L11*L32.powi(2)*L33 + 2.0*L12*L31*L32*L33 + L13*L31*L32.powi(2);
    MG[(7,6)]=2.0*L11*L22*L32*L33 + L11*L23*L32.powi(2) + 2.0*L12*L21*L32*L33 + 2.0*L12*L22*L31*L33 + 2.0*L12*L23*L31*L32 + L13*L21*L32.powi(2) + 2.0*L13*L22*L31*L32;
    MG[(7,7)]=L11*L22.powi(2)*L33 + 2.0*L11*L22*L23*L32 + 2.0*L12*L21*L22*L33 + 2.0*L12*L21*L23*L32 + 2.0*L12*L22*L23*L31 + 2.0*L13*L21*L22*L32 + L13*L22.powi(2)*L31;
    MG[(7,8)]=L11*L22.powi(2)*L23 + 2.0*L12*L21*L22*L23 + L13*L21*L22.powi(2);
    MG[(7,9)]=2.0*L11*L12*L32*L33 + L11*L13*L32.powi(2) + L12.powi(2)*L31*L33 + 2.0*L12*L13*L31*L32;
    MG[(7,10)]=2.0*L11*L12*L22*L33 + 2.0*L11*L12*L23*L32 + 2.0*L11*L13*L22*L32 + L12.powi(2)*L21*L33 + L12.powi(2)*L23*L31 + 2.0*L12*L13*L21*L32 + 2.0*L12*L13*L22*L31;
    MG[(7,11)]=2.0*L11*L12*L22*L23 + L11*L13*L22.powi(2) + L12.powi(2)*L21*L23 + 2.0*L12*L13*L21*L22;
    MG[(7,12)]=L11*L12.powi(2)*L33 + 2.0*L11*L12*L13*L32 + L12.powi(2)*L13*L31;
    MG[(7,13)]=L11*L12.powi(2)*L23 + 2.0*L11*L12*L13*L22 + L12.powi(2)*L13*L21;
    MG[(7,14)]=L11*L12.powi(2)*L13;
    MG[(8,0)]=L31*L32.powi(3);
    MG[(8,1)]=L21*L32.powi(3) + 3.0*L22*L31*L32.powi(2);
    MG[(8,2)]=3.0*L21*L22*L32.powi(2) + 3.0*L22.powi(2)*L31*L32;
    MG[(8,3)]=3.0*L21*L22.powi(2)*L32 + L22.powi(3)*L31;
    MG[(8,4)]=L21*L22.powi(3);
    MG[(8,5)]=L11*L32.powi(3) + 3.0*L12*L31*L32.powi(2);
    MG[(8,6)]=3.0*L11*L22*L32.powi(2) + 3.0*L12*L21*L32.powi(2) + 6.0*L12*L22*L31*L32;
    MG[(8,7)]=3.0*L11*L22.powi(2)*L32 + 6.0*L12*L21*L22*L32 + 3.0*L12*L22.powi(2)*L31;
    MG[(8,8)]=L11*L22.powi(3) + 3.0*L12*L21*L22.powi(2);
    MG[(8,9)]=3.0*L11*L12*L32.powi(2) + 3.0*L12.powi(2)*L31*L32;
    MG[(8,10)]=6.0*L11*L12*L22*L32 + 3.0*L12.powi(2)*L21*L32 + 3.0*L12.powi(2)*L22*L31;
    MG[(8,11)]=3.0*L11*L12*L22.powi(2) + 3.0*L12.powi(2)*L21*L22;
    MG[(8,12)]=3.0*L11*L12.powi(2)*L32 + L12.powi(3)*L31;
    MG[(8,13)]=3.0*L11*L12.powi(2)*L22 + L12.powi(3)*L21;
    MG[(8,14)]=L11*L12.powi(3);
    MG[(9,0)]=L31.powi(2)*L33.powi(2);
    MG[(9,1)]=2.0*L21*L31*L33.powi(2) + 2.0*L23*L31.powi(2)*L33;
    MG[(9,2)]=L21.powi(2)*L33.powi(2) + 4.0*L21*L23*L31*L33 + L23.powi(2)*L31.powi(2);
    MG[(9,3)]=2.0*L21.powi(2)*L23*L33 + 2.0*L21*L23.powi(2)*L31;
    MG[(9,4)]=L21.powi(2)*L23.powi(2);
    MG[(9,5)]=2.0*L11*L31*L33.powi(2) + 2.0*L13*L31.powi(2)*L33;
    MG[(9,6)]=2.0*L11*L21*L33.powi(2) + 4.0*L11*L23*L31*L33 + 4.0*L13*L21*L31*L33 + 2.0*L13*L23*L31.powi(2);
    MG[(9,7)]=4.0*L11*L21*L23*L33 + 2.0*L11*L23.powi(2)*L31 + 2.0*L13*L21.powi(2)*L33 + 4.0*L13*L21*L23*L31;
    MG[(9,8)]=2.0*L11*L21*L23.powi(2) + 2.0*L13*L21.powi(2)*L23;
    MG[(9,9)]=L11.powi(2)*L33.powi(2) + 4.0*L11*L13*L31*L33 + L13.powi(2)*L31.powi(2);
    MG[(9,10)]=2.0*L11.powi(2)*L23*L33 + 4.0*L11*L13*L21*L33 + 4.0*L11*L13*L23*L31 + 2.0*L13.powi(2)*L21*L31;
    MG[(9,11)]=L11.powi(2)*L23.powi(2) + 4.0*L11*L13*L21*L23 + L13.powi(2)*L21.powi(2);
    MG[(9,12)]=2.0*L11.powi(2)*L13*L33 + 2.0*L11*L13.powi(2)*L31;
    MG[(9,13)]=2.0*L11.powi(2)*L13*L23 + 2.0*L11*L13.powi(2)*L21;
    MG[(9,14)]=L11.powi(2)*L13.powi(2);
    MG[(10,0)]=L31.powi(2)*L32*L33;
    MG[(10,1)]=2.0*L21*L31*L32*L33 + L22*L31.powi(2)*L33 + L23*L31.powi(2)*L32;
    MG[(10,2)]=L21.powi(2)*L32*L33 + 2.0*L21*L22*L31*L33 + 2.0*L21*L23*L31*L32 + L22*L23*L31.powi(2);
    MG[(10,3)]=L21.powi(2)*L22*L33 + L21.powi(2)*L23*L32 + 2.0*L21*L22*L23*L31;
    MG[(10,4)]=L21.powi(2)*L22*L23;
    MG[(10,5)]=2.0*L11*L31*L32*L33 + L12*L31.powi(2)*L33 + L13*L31.powi(2)*L32;
    MG[(10,6)]=2.0*L11*L21*L32*L33 + 2.0*L11*L22*L31*L33 + 2.0*L11*L23*L31*L32 + 2.0*L12*L21*L31*L33 + L12*L23*L31.powi(2) + 2.0*L13*L21*L31*L32 + L13*L22*L31.powi(2);
    MG[(10,7)]=2.0*L11*L21*L22*L33 + 2.0*L11*L21*L23*L32 + 2.0*L11*L22*L23*L31 + L12*L21.powi(2)*L33 + 2.0*L12*L21*L23*L31 + L13*L21.powi(2)*L32 + 2.0*L13*L21*L22*L31;
    MG[(10,8)]=2.0*L11*L21*L22*L23 + L12*L21.powi(2)*L23 + L13*L21.powi(2)*L22;
    MG[(10,9)]=L11.powi(2)*L32*L33 + 2.0*L11*L12*L31*L33 + 2.0*L11*L13*L31*L32 + L12*L13*L31.powi(2);
    MG[(10,10)]=L11.powi(2)*L22*L33 + L11.powi(2)*L23*L32 + 2.0*L11*L12*L21*L33 + 2.0*L11*L12*L23*L31 + 2.0*L11*L13*L21*L32 + 2.0*L11*L13*L22*L31 + 2.0*L12*L13*L21*L31;
    MG[(10,11)]=L11.powi(2)*L22*L23 + 2.0*L11*L12*L21*L23 + 2.0*L11*L13*L21*L22 + L12*L13*L21.powi(2);
    MG[(10,12)]=L11.powi(2)*L12*L33 + L11.powi(2)*L13*L32 + 2.0*L11*L12*L13*L31;
    MG[(10,13)]=L11.powi(2)*L12*L23 + L11.powi(2)*L13*L22 + 2.0*L11*L12*L13*L21;
    MG[(10,14)]=L11.powi(2)*L12*L13;
    MG[(11,0)]=L31.powi(2)*L32.powi(2);
    MG[(11,1)]=2.0*L21*L31*L32.powi(2) + 2.0*L22*L31.powi(2)*L32;
    MG[(11,2)]=L21.powi(2)*L32.powi(2) + 4.0*L21*L22*L31*L32 + L22.powi(2)*L31.powi(2);
    MG[(11,3)]=2.0*L21.powi(2)*L22*L32 + 2.0*L21*L22.powi(2)*L31;
    MG[(11,4)]=L21.powi(2)*L22.powi(2);
    MG[(11,5)]=2.0*L11*L31*L32.powi(2) + 2.0*L12*L31.powi(2)*L32;
    MG[(11,6)]=2.0*L11*L21*L32.powi(2) + 4.0*L11*L22*L31*L32 + 4.0*L12*L21*L31*L32 + 2.0*L12*L22*L31.powi(2);
    MG[(11,7)]=4.0*L11*L21*L22*L32 + 2.0*L11*L22.powi(2)*L31 + 2.0*L12*L21.powi(2)*L32 + 4.0*L12*L21*L22*L31;
    MG[(11,8)]=2.0*L11*L21*L22.powi(2) + 2.0*L12*L21.powi(2)*L22;
    MG[(11,9)]=L11.powi(2)*L32.powi(2) + 4.0*L11*L12*L31*L32 + L12.powi(2)*L31.powi(2);
    MG[(11,10)]=2.0*L11.powi(2)*L22*L32 + 4.0*L11*L12*L21*L32 + 4.0*L11*L12*L22*L31 + 2.0*L12.powi(2)*L21*L31;
    MG[(11,11)]=L11.powi(2)*L22.powi(2) + 4.0*L11*L12*L21*L22 + L12.powi(2)*L21.powi(2);
    MG[(11,12)]=2.0*L11.powi(2)*L12*L32 + 2.0*L11*L12.powi(2)*L31;
    MG[(11,13)]=2.0*L11.powi(2)*L12*L22 + 2.0*L11*L12.powi(2)*L21;
    MG[(11,14)]=L11.powi(2)*L12.powi(2);
    MG[(12,0)]=L31.powi(3)*L33;
    MG[(12,1)]=3.0*L21*L31.powi(2)*L33 + L23*L31.powi(3);
    MG[(12,2)]=3.0*L21.powi(2)*L31*L33 + 3.0*L21*L23*L31.powi(2);
    MG[(12,3)]=L21.powi(3)*L33 + 3.0*L21.powi(2)*L23*L31;
    MG[(12,4)]=L21.powi(3)*L23;
    MG[(12,5)]=3.0*L11*L31.powi(2)*L33 + L13*L31.powi(3);
    MG[(12,6)]=6.0*L11*L21*L31*L33 + 3.0*L11*L23*L31.powi(2) + 3.0*L13*L21*L31.powi(2);
    MG[(12,7)]=3.0*L11*L21.powi(2)*L33 + 6.0*L11*L21*L23*L31 + 3.0*L13*L21.powi(2)*L31;
    MG[(12,8)]=3.0*L11*L21.powi(2)*L23 + L13*L21.powi(3);
    MG[(12,9)]=3.0*L11.powi(2)*L31*L33 + 3.0*L11*L13*L31.powi(2);
    MG[(12,10)]=3.0*L11.powi(2)*L21*L33 + 3.0*L11.powi(2)*L23*L31 + 6.0*L11*L13*L21*L31;
    MG[(12,11)]=3.0*L11.powi(2)*L21*L23 + 3.0*L11*L13*L21.powi(2);
    MG[(12,12)]=L11.powi(3)*L33 + 3.0*L11.powi(2)*L13*L31;
    MG[(12,13)]=L11.powi(3)*L23 + 3.0*L11.powi(2)*L13*L21;
    MG[(12,14)]=L11.powi(3)*L13;
    MG[(13,0)]=L31.powi(3)*L32;
    MG[(13,1)]=3.0*L21*L31.powi(2)*L32 + L22*L31.powi(3);
    MG[(13,2)]=3.0*L21.powi(2)*L31*L32 + 3.0*L21*L22*L31.powi(2);
    MG[(13,3)]=L21.powi(3)*L32 + 3.0*L21.powi(2)*L22*L31;
    MG[(13,4)]=L21.powi(3)*L22;
    MG[(13,5)]=3.0*L11*L31.powi(2)*L32 + L12*L31.powi(3);
    MG[(13,6)]=6.0*L11*L21*L31*L32 + 3.0*L11*L22*L31.powi(2) + 3.0*L12*L21*L31.powi(2);
    MG[(13,7)]=3.0*L11*L21.powi(2)*L32 + 6.0*L11*L21*L22*L31 + 3.0*L12*L21.powi(2)*L31;
    MG[(13,8)]=3.0*L11*L21.powi(2)*L22 + L12*L21.powi(3);
    MG[(13,9)]=3.0*L11.powi(2)*L31*L32 + 3.0*L11*L12*L31.powi(2);
    MG[(13,10)]=3.0*L11.powi(2)*L21*L32 + 3.0*L11.powi(2)*L22*L31 + 6.0*L11*L12*L21*L31;
    MG[(13,11)]=3.0*L11.powi(2)*L21*L22 + 3.0*L11*L12*L21.powi(2);
    MG[(13,12)]=L11.powi(3)*L32 + 3.0*L11.powi(2)*L12*L31;
    MG[(13,13)]=L11.powi(3)*L22 + 3.0*L11.powi(2)*L12*L21;
    MG[(13,14)]=L11.powi(3)*L12;
    MG[(14,0)]=L31.powi(4);
    MG[(14,1)]=4.0*L21*L31.powi(3);
    MG[(14,2)]=6.0*L21.powi(2)*L31.powi(2);
    MG[(14,3)]=4.0*L21.powi(3)*L31;
    MG[(14,4)]=L21.powi(4);
    MG[(14,5)]=4.0*L11*L31.powi(3);
    MG[(14,6)]=12.0*L11*L21*L31.powi(2);
    MG[(14,7)]=12.0*L11*L21.powi(2)*L31;
    MG[(14,8)]=4.0*L11*L21.powi(3);
    MG[(14,9)]=6.0*L11.powi(2)*L31.powi(2);
    MG[(14,10)]=12.0*L11.powi(2)*L21*L31;
    MG[(14,11)]=6.0*L11.powi(2)*L21.powi(2);
    MG[(14,12)]=4.0*L11.powi(3)*L31;
    MG[(14,13)]=4.0*L11.powi(3)*L21;
    MG[(14,14)]=L11.powi(4);
    MG
}

pub fn get_nmat(lmns: &Vec<[usize; 3]>) -> Array2<f64> {
    let nval = |lmn: &[usize; 3]| {
        let [l, m, n] = lmn;
        let fac = FACS[*l] * FACS[*m] * FACS[*n];
        1.0 / fac.sqrt()
    };
    let size = lmns.len();
    let mut nmat = Array2::zeros([size, size]);
    for (i, lmn_i) in lmns.iter().enumerate() {
        let val_i = nval(lmn_i);
        for (j, lmn_j) in lmns.iter().enumerate() {
            let val_j = nval(lmn_j);
            nmat[(i, j)] = val_i / val_j;
        }
    }
    nmat
}

#[rustfmt::skip]
/// 获取 d 轨道的归一化矩阵
pub fn nmat_d()->Array2<f64>{
    let lmns=btrans::CAR_D_LMNS.to_vec();
    get_nmat(&lmns)
}

#[rustfmt::skip]
/// 获取 f 轨道的归一化矩阵
pub fn nmat_f()->Array2<f64>{
    let lmns=btrans::CAR_F_LMNS.to_vec();
    get_nmat(&lmns)
}

#[rustfmt::skip]
/// 获取 g 轨道的归一化矩阵
pub fn nmat_g()->Array2<f64>{
    let lmns=btrans::CAR_G_LMNS.to_vec();
    get_nmat(&lmns)
}

/// 计算归一化笛卡尔基函数在同一壳层内的角度重叠矩阵 S
///
/// 对同一角动量 l 的壳层，所有函数共用指数 α，归一化后的笛卡尔函数
/// χ_i = N_i x^{li} y^{mi} z^{ni} e^{-αr²} 之间的重叠
///
///   S[i][j] = (li+lj-1)!!(mi+mj-1)!!(ni+nj-1)!! / sqrt(fac_i·fac_j)
///
/// 只与角量子数有关，与 α 无关，是纯角度量。
fn get_sang(lmns: &[[usize; 3]]) -> Array2<f64> {
    let n = lmns.len();
    let mut s = Array2::<f64>::zeros((n, n));
    for i in 0..n {
        let [li, mi, ni] = lmns[i];
        let fac_i = FACS[li] * FACS[mi] * FACS[ni];
        for j in 0..n {
            let [lj, mj, nj] = lmns[j];
            let lx = li + lj;
            let my = mi + mj;
            let nz = ni + nj;
            // 任意一个方向为奇次幂时，高斯积分恒为零
            if lx % 2 == 1 || my % 2 == 1 || nz % 2 == 1 {
                continue;
            }
            let num = FACS[lx / 2] * FACS[my / 2] * FACS[nz / 2];
            let fac_j = FACS[lj] * FACS[mj] * FACS[nj];
            s[(i, j)] = num / (fac_i * fac_j).sqrt();
        }
    }
    s
}

/// 归一化笛卡尔 d 基的角度重叠矩阵
fn sang_d() -> Array2<f64> {
    get_sang(&btrans::CAR_D_LMNS)
}

/// 归一化笛卡尔 f 基的角度重叠矩阵
fn sang_f() -> Array2<f64> {
    get_sang(&btrans::CAR_F_LMNS)
}

/// 归一化笛卡尔 g 基的角度重叠矩阵
fn sang_g() -> Array2<f64> {
    get_sang(&btrans::CAR_G_LMNS)
}
