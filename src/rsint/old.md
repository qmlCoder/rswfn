```rs
//计算高斯函数积分：⾼斯-埃尔⽶特求积法
pub fn gtf_integ(alps: &[f64; 2], xyzs: &[[f64; 3]; 2], lmns: &[[usize; 3]; 2]) -> f64 {
    let [e1, e2] = alps;
    let [x1, y1, z1] = &xyzs[0];
    let [x2, y2, z2] = &xyzs[1];
    let alp = e1 + e2;
    let px = (e1 * x1 + e2 * x2) / alp;
    let py = (e1 * y1 + e2 * y2) / alp;
    let pz = (e1 * z1 + e2 * z2) / alp;
    let [l1, m1, n1] = lmns[0];
    let [l2, m2, n2] = lmns[1];
    let dx = x2 - x1;
    let dy = y2 - y1;
    let dz = z2 - z1;
    let expx = -e1 * e2 * (dx.powi(2) + dy.powi(2) + dz.powi(2)) / alp;
    let expv = f64::exp(expx);
    let nx = ((l1 + l2 + 1) as f64 / 2.0).ceil() as usize;
    let ny = ((m1 + m2 + 1) as f64 / 2.0).ceil() as usize;
    let nz = ((n1 + n2 + 1) as f64 / 2.0).ceil() as usize;
    let sqrv = alp.sqrt();
    let mut sx = 0.0;
    for i in 0..nx {
        let tmpv = RHM[(nx - 1, i)] / sqrv + px;
        sx += WHM[(nx - 1, i)] * (tmpv - x1).powi(l1 as i32) * (tmpv - x2).powi(l2 as i32);
    }
    sx = sx / sqrv;

    let mut sy = 0.0;
    for i in 0..ny {
        let tmpv = RHM[(ny - 1, i)] / sqrv + py;
        sy += WHM[(ny - 1, i)] * (tmpv - y1).powi(m1 as i32) * (tmpv - y2).powi(m2 as i32);
    }
    sy = sy / sqrv;

    let mut sz = 0.0;
    for i in 0..nz {
        let tmpv = RHM[(nz - 1, i)] / sqrv + pz;
        sz += WHM[(nz - 1, i)] * (tmpv - z1).powi(n1 as i32) * (tmpv - z2).powi(n2 as i32);
    }
    sz = sz / sqrv;
    let val = sx * sy * sz * expv;
    val
}

pub static CAR_LMNS: LazyLock<Vec<Vec<[usize; 3]>>> = LazyLock::new(|| vec![btrans::MIX_S_LMNS.clone().to_vec(), btrans::MIX_P_LMNS.clone().to_vec(), btrans::CAR_D_LMNS.clone().to_vec(), btrans::CAR_F_LMNS.clone().to_vec(), btrans::CAR_G_LMNS.clone().to_vec(), btrans::CAR_H_LMNS.clone().to_vec()]);

pub fn lmn2vec(lmn: &[usize; 3]) -> (Vec<f64>, Vec<[usize; 3]>) {
    let mut vec;
    for lmns in CAR_LMNS.iter() {
        vec = vec![0.0; lmns.len()];
        for (j, each) in lmns.iter().enumerate() {
            if each == lmn {
                vec[j] = 1.0;
                return (vec, lmns.clone());
            }
        }
    }
    panic!("lmn2vec error");
}

pub fn sym2vec(sym: &str) -> (Vec<f64>, Vec<[usize; 3]>) {
    if btrans::MIX_S_SYMS.contains(&sym) {
        let vec = vec![1.0];
        return (vec, btrans::MIX_S_LMNS.clone().to_vec());
    }
    if btrans::MIX_P_SYMS.contains(&sym) {
        let idx = btrans::MIX_P_SYMS.iter().position(|x| *x == sym).unwrap();
        let mut vec = vec![0.0; 3];
        vec[idx] = 1.0;
        return (vec, btrans::MIX_P_LMNS.clone().to_vec());
    }
    if btrans::CAR_D_SYMS.contains(&sym) {
        let idx = btrans::CAR_D_SYMS.iter().position(|x| *x == sym).unwrap();
        let mut vec = vec![0.0; 6];
        vec[idx] = 1.0;
        return (vec, btrans::CAR_D_LMNS.clone().to_vec());
    }
    if btrans::CAR_F_SYMS.contains(&sym) {
        let idx = btrans::CAR_F_SYMS.iter().position(|x| *x == sym).unwrap();
        let mut vec = vec![0.0; 10];
        vec[idx] = 1.0;
        return (vec, btrans::CAR_F_LMNS.clone().to_vec());
    }
    if btrans::SPH_D_SYMS.contains(&sym) {
        let idx = btrans::SPH_D_SYMS.iter().position(|x| *x == sym).unwrap();
        let mut vec = vec![0.0; 5];
        vec[idx] = 1.0;
        let vec = Array2::from_shape_vec((vec.len(), 1), vec).expect("生成向量失败");
        let vec = btrans::DMAT.dot(&vec).flatten().to_vec();
        return (vec, btrans::CAR_D_LMNS.clone().to_vec());
    }
    if btrans::SPH_F_SYMS.contains(&sym) {
        let idx = btrans::SPH_F_SYMS.iter().position(|x| *x == sym).unwrap();
        let mut vec = vec![0.0; 7];
        vec[idx] = 1.0;
        let vec = Array2::from_shape_vec((vec.len(), 1), vec).expect("生成向量失败");
        let vec = btrans::DMAT.dot(&vec).flatten().to_vec();
        return (vec, btrans::CAR_F_LMNS.clone().to_vec());
    }
    panic!("sym2vec error:{}", sym);
}

pub fn local_gtf_integ(
    alps: &[f64; 2],
    xyzs: &[[f64; 3]; 2],
    syms: &[&str; 2],
    stms: &[Array2<f64>; 2], // 两个局部坐标系，每一列代表一个单位向量
) -> f64 {
    let sym_i = &syms[0];
    let sym_j = &syms[1];
    let stm_i = &stms[0];
    let stm_j = &stms[1];
    let (vec_i, lmns_i) = sym2vec(&sym_i);
    let (vec_j, lmns_j) = sym2vec(&sym_j);
    let tmat_i = match vec_i.len() {
        1 => mocv::tmat_s(&stm_i),
        3 => mocv::tmat_p(&stm_i),
        6 => mocv::tmat_d(&stm_i),
        10 => mocv::tmat_f(&stm_i),
        _ => panic!("lmn2vec error:{:?}", vec_i),
    };
    let tmat_j = match vec_j.len() {
        1 => mocv::tmat_s(&stm_i),
        3 => mocv::tmat_p(&stm_j),
        6 => mocv::tmat_d(&stm_j),
        10 => mocv::tmat_f(&stm_j),
        _ => panic!("lmn2vec error:{:?}", vec_j),
    };
    let vec_i = Array2::from_shape_vec((vec_i.len(), 1), vec_i.to_vec()).expect("生成向量失败");
    let vec_j = Array2::from_shape_vec((vec_j.len(), 1), vec_j.to_vec()).expect("生成向量失败");
    let coes_i = tmat_i.t().dot(&vec_i).flatten().to_vec();
    let coes_j = tmat_j.t().dot(&vec_j).flatten().to_vec();
    let len_i = vec_i.len();
    let len_j = vec_j.len();
    let mut val = 0.0;
    for i in 0..len_i {
        for j in 0..len_j {
            val += gtf_integ(alps, xyzs, &[lmns_i[i], lmns_j[j]]) * coes_i[i] * coes_j[j];
        }
    }
    val
}

// 重叠矩阵，全都是正交基函数
pub fn calc_smat(
    alps: &Vec<Vec<f64>>,   // 每个原子轨道的 指数列表
    coes: &Vec<Vec<f64>>,   // 每个原子轨道的 系数列表
    lmns: &Vec<[usize; 3]>, // 每个原子轨道 的角动量分量
    xyzs: &Vec<[f64; 3]>,   // 每个原子轨道的坐标
) -> Array2<f64> {
    let nato = xyzs.len(); // 原子轨道数量
    let mut smat = Array2::zeros((nato, nato)); // 重叠矩阵
    let mut coes = coes.clone();
    for i in 0..coes.len() {
        let [l, m, n] = lmns[i];
        for j in 0..coes[i].len() {
            let fac = FACS[l as usize] * FACS[m as usize] * FACS[n as usize];
            let ang = l + m + n;
            let alp = alps[i][j];
            let nm = (2.0 * alp / PI).powf(0.75) * f64::sqrt((4.0 * alp).powi(ang as i32) / fac);
            coes[i][j] = nm * coes[i][j];
        }
    }
    for i in 0..nato {
        let nbas_i = alps[i].len();
        for j in 0..nato {
            let nbas_j = alps[j].len();

            for k in 0..nbas_i {
                for l in 0..nbas_j {
                    let sval = gtf_integ(&[alps[i][k], alps[j][l]], &[xyzs[i], xyzs[j]], &[lmns[i], lmns[j]]);
                    smat[(i, j)] += coes[i][k] * coes[j][l] * sval;
                    assert_ne!(coes[i][k], 0.0);
                    assert_ne!(coes[j][l], 0.0);
                }
            }
        }
    }
    smat
}
```
