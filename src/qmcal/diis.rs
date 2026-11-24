//! Pulay DIIS（Direct Inversion in the Iterative Subspace）
//!
//! 把最近若干步的 Fock 矩阵做线性组合来加速并稳定 SCF 收敛：组合系数由
//! “让误差向量（Fock 与密度矩阵的交换子）的模最小”确定，约束为系数之和为 1。
//!
//! 相比纯阻尼迭代，DIIS 通常能把迭代次数减半，并能消掉阻尼迭代常见的
//! 两态振荡。

use ndarray::Array2;

use crate::maths;

/// Pulay DIIS 外推器
pub struct Diis {
    /// 保留的历史步数上限
    pub max_len: usize,
    /// 历史 Fock 矩阵
    focks: Vec<Array2<f64>>,
    /// 历史误差向量
    errors: Vec<Array2<f64>>,
}

impl Diis {
    /// 创建 DIIS 外推器（`max_len` 为历史长度，常用 6~10）
    pub fn new(max_len: usize) -> Self {
        Self {
            max_len: max_len.max(2),
            focks: Vec::new(),
            errors: Vec::new(),
        }
    }

    /// 清空历史（SCF 发散后重新开始）
    pub fn reset(&mut self) {
        self.focks.clear();
        self.errors.clear();
    }

    /// SCF 误差向量：Fock 矩阵与密度矩阵的交换子 `F D S - S D F`
    ///
    /// 自洽时 Fock 与密度矩阵可对易，该量趋于零，适合作为收敛判据的“方向”。
    pub fn commutator(fock: &Array2<f64>, dmat: &Array2<f64>, smat: &Array2<f64>) -> Array2<f64> {
        let fds = fock.dot(dmat).dot(smat);
        let sdf = smat.dot(dmat).dot(fock);
        &fds - &sdf
    }

    /// 推入当前 Fock 矩阵与误差向量，返回外推后的 Fock 矩阵
    ///
    /// 历史不足（只有一步）或线性方程组数值条件太差时返回 `None`，
    /// 调用方应退回普通迭代（配合阻尼）。
    pub fn extrapolate(&mut self, fock: &Array2<f64>, err: &Array2<f64>) -> Option<Array2<f64>> {
        self.focks.push(fock.clone());
        self.errors.push(err.clone());
        while self.focks.len() > self.max_len {
            self.focks.remove(0);
            self.errors.remove(0);
        }
        loop {
            if self.focks.len() < 2 {
                return None;
            }
            if let Some(coef) = self.solve() {
                let mut fmat = Array2::<f64>::zeros(fock.raw_dim());
                for (i, fock_i) in self.focks.iter().enumerate() {
                    fmat.scaled_add(coef[i], fock_i);
                }
                return Some(fmat);
            }
            // 条件数太差：丢掉最旧的一步再试
            self.focks.remove(0);
            self.errors.remove(0);
        }
    }

    /// 求解 DIIS 线性方程组，得到组合系数
    ///
    /// 求解 `[[e_i·e_j, 1], [1^T, 0]] [c, λ] = [0, 1]`。
    fn solve(&self) -> Option<Vec<f64>> {
        let n = self.focks.len();
        let mut bmat = Array2::<f64>::zeros((n + 1, n + 1));
        for i in 0..n {
            for j in 0..n {
                bmat[(i, j)] = frob(&self.errors[i], &self.errors[j]);
            }
            bmat[(i, n)] = 1.0;
            bmat[(n, i)] = 1.0;
        }
        // 轻微正则化：误差向量线性相关时避免矩阵奇异
        let trace: f64 = (0..n).map(|i| bmat[(i, i)]).sum();
        let ridge = if trace > 0.0 {
            trace / n as f64 * 1.0e-10
        } else {
            1.0e-12
        };
        for i in 0..n {
            bmat[(i, i)] += ridge;
        }

        let mut rhs = Array2::<f64>::zeros((n + 1, 1));
        rhs[(n, 0)] = 1.0;
        let sol = maths::inv(&bmat).dot(&rhs);
        let coef: Vec<f64> = (0..n).map(|i| sol[(i, 0)]).collect();

        let finite = coef.iter().all(|v| v.is_finite());
        let bounded = coef.iter().map(|v| v.abs()).sum::<f64>() < 10.0;
        if finite && bounded { Some(coef) } else { None }
    }
}

/// 两个矩阵的 Frobenius 内积
fn frob(a: &Array2<f64>, b: &Array2<f64>) -> f64 {
    a.iter().zip(b.iter()).map(|(x, y)| x * y).sum()
}
