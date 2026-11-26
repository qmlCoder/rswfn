# 无约束优化算法

`rswfn::maths::optim` 提供一组求解无约束极小化问题的迭代算法：

```text
minimize F(x)
```

其中 `x` 是长度为 `n` 的自变量向量。目前包含：

- 最速下降法（Steepest Descent）
- 非线性共轭梯度法（Conjugate Gradient）
- 牛顿法（Newton Method）
- 阻尼牛顿法（Damped Newton Method）
- 拟牛顿法（BFGS）

所有算法共用同一套目标函数接口和收敛判据，通常只需要把算法名称替换一下即可切换。

## 目录结构

```text
src/optim/
├── mod.rs    # 公共 API：FuncEval、算法函数导出
├── utils.rs  # 数值导数、线搜索、向量运算等内部工具
├── sd.rs     # 最速下降法
├── cg.rs     # 非线性共轭梯度法
├── nt.rs     # 牛顿法
├── dn.rs     # 阻尼牛顿法
└── qn.rs     # 拟牛顿法（BFGS）
```

## 公共接口

每个算法模块都提供同名函数 `minimize`，签名一致：

```rust
pub fn minimize<F>(func: F, paras: &Vec<f64>, conv_thresh: f64) -> Result<Vec<f64>, String>
where
    F: Fn(&Vec<f64>) -> FuncEval,
```

- `func`：目标函数闭包，输入为自变量 `&Vec<f64>`，返回 `FuncEval`。
- `paras`：初始自变量。
- `conv_thresh`：收敛阈值。
- 返回值：收敛后的自变量；若失败或达到最大迭代次数仍未收敛，返回 `Err(String)`。

顶层还提供语义更明确的别名：

| 算法 | 顶层函数 | 模块函数 |
| --- | --- | --- |
| 最速下降法 | `optim::steepest_descent` | `optim::sd::minimize` |
| 非线性共轭梯度法 | `optim::conjugate_gradient` | `optim::cg::minimize` |
| 牛顿法 | `optim::newton` | `optim::nt::minimize` |
| 阻尼牛顿法 | `optim::damped_newton` | `optim::dn::minimize` |
| 拟牛顿法（BFGS） | `optim::quasi_newton` | `optim::qn::minimize` |

## 目标函数：`FuncEval`

`FuncEval` 保存一次函数求值的结果。函数值必须提供，梯度和 Hessian 可选：

```rust
use ndarray::Array2;
use rswfn::maths::optim::FuncEval;

pub struct FuncEval {
    pub value: f64,
    pub grad: Option<Vec<f64>>,
    pub hessian: Option<Array2<f64>>,
}
```

常用构造方式：

```rust
// 只提供函数值：算法自动用中心差分求梯度
FuncEval::from_value(value)

// 提供函数值和解析梯度
FuncEval::from_grad(value, grad)

// 提供函数值和解析 Hessian
FuncEval::from_hessian(value, hessian)

// 同时提供梯度和 Hessian
FuncEval::new(value, Some(grad), Some(hessian))
```

其中 `grad` 的长度必须等于自变量维数，`hessian` 的形状必须为 `(n, n)`。牛顿法和阻尼牛顿法会严格检查 Hessian；其他算法只在需要 Hessian 估计初始步长时使用它。若只提供 Hessian 而没有梯度，梯度仍会由函数值中心差分得到。

### 数值导数

当用户没有提供解析导数时，算法会自动做有限差分：

- 梯度：中心差分，步长约为 `EPSILON^(1/3) * max(1, |x_i|)`。
- Hessian：
  - 若已有解析梯度，则对梯度做中心差分，约需 `2n` 次函数求值；
  - 否则对函数值做二阶中心差分，约需 `O(n^2)` 次函数求值。

因此，若只能提供函数值，数值 Hessian 的代价会随维数平方增长；对于高维问题，优先提供解析梯度或 Hessian。

## 算法介绍

### 1. 最速下降法（SD）

文件：[`sd.rs`](sd.rs)

每步沿负梯度方向搜索：

```text
p_k = -g_k
x_{k+1} = x_k + alpha_k p_k
```

步长 `alpha_k` 通过 Armijo 回溯线搜索确定。若用户提供了 Hessian，算法会用它估计一个二次模型下的初始步长，但不会用 Hessian 改变搜索方向。

特点：

- 实现简单、单步计算量小、内存占用低。
- 对良态问题稳定，但对病态问题容易出现“锯齿”现象，收敛较慢。
- 通常为线性收敛。

适用场景：作为基线算法，或问题维度较高、只希望使用梯度信息时。

### 2. 非线性共轭梯度法（CG）

文件：[`cg.rs`](cg.rs)

使用 Polak–Ribière+ 公式构造共轭方向：

```text
beta_k = max(0, g_k · (g_k - g_{k-1}) / (g_{k-1} · g_{k-1}))
p_k = -g_k + beta_k p_{k-1}
```

线搜索采用强 Wolfe 条件：

```text
F(x + alpha p) <= F(x) + c1 * alpha * g·p
|g(x + alpha p)·p| <= c2 * |g·p|
```

若强 Wolfe 线搜索失败，则回退到 Armijo 回溯线搜索；若新方向不是下降方向，则自动重新开始（退化为负梯度方向）。

特点：

- 只需要梯度，不需要 Hessian。
- 内存占用为 `O(n)`，适合较大规模问题。
- 对二次型，在精确线搜索下理论上具有至多 `n` 步收敛性；非线性问题中通常表现为较快的线性收敛，在满足一定条件时可达超线性。

适用场景：中大规模问题、梯度容易获得但 Hessian 难以构造或存储的问题。

### 3. 牛顿法（NT）

文件：[`nt.rs`](nt.rs)

使用完整牛顿步：

```text
H_k p_k = -g_k
x_{k+1} = x_k + p_k
```

其中 `H_k` 是当前点的 Hessian。算法通过特征分解求解线性方程组，不显式求逆。

特点：

- 在极小点附近具有二次收敛速度。
- 不做线搜索，采用完整步长。
- 要求初始点足够接近极小点，且 Hessian 正定、可逆。
- 若 Hessian 奇异、牛顿方向不是下降方向，或步后函数值非有限，会返回 `Err`。

适用场景：已有可靠 Hessian，且初始点已经比较接近极小点的情况。

### 4. 阻尼牛顿法（DN）

文件：[`dn.rs`](dn.rs)

阻尼牛顿法在牛顿方向上进行线搜索。若 Hessian 非正定，则通过特征值平移构造正定矩阵：

```text
(H_k + shift * I) p_k = -g_k
x_{k+1} = x_k + alpha_k p_k
```

其中当 `lambda_min <= 0` 时，`shift = -lambda_min + 1e-8 * scale`；否则 `shift = 0`。`scale` 是 Hessian 特征值绝对值的最大值（至少为 1）。随后用 Armijo 回溯线搜索确定 `alpha_k`。

特点：

- 保留牛顿法的快速局部收敛性，同时比纯牛顿法更稳健。
- 当 Hessian 非正定时仍能构造下降方向。
- 如果修正后的方向仍不是下降方向，会回退到负梯度方向。
- 计算代价与牛顿法相近，需要对 Hessian 做特征分解。

适用场景：希望使用二阶信息，但初始点可能不在极小点附近，或 Hessian 可能非正定的情况。

### 5. 拟牛顿法（BFGS）

文件：[`qn.rs`](qn.rs)

拟牛顿法不直接计算 Hessian，而是用梯度差商逐步构造逆 Hessian 的近似矩阵 `B`。搜索方向为：

```text
p_k = -B_k g_k
```

每次迭代根据位移和梯度变化更新 `B`：

```text
s_k = x_{k+1} - x_k
y_k = g_{k+1} - g_k
rho_k = 1 / (y_k · s_k)

B_{k+1} = (I - rho_k s_k y_k^T) B_k (I - rho_k y_k s_k^T) + rho_k s_k s_k^T
```

线搜索采用强 Wolfe 条件，失败时回退到 Armijo 回溯线搜索。

实现中的保护措施：

- 初始 `B_0 = I`。
- 当 `y_k · s_k` 不满足曲率条件时，跳过本次更新。
- 若数值误差导致搜索方向不是下降方向，则把 `B` 重置为单位矩阵，并使用负梯度方向。
- 若更新后 `B` 出现非有限值，也会重置为单位矩阵。

特点：

- 核心迭代只需要梯度，不需要显式 Hessian。
- 通常具有超线性收敛速度，是通用性很好的选择。
- 需要存储 `n x n` 的逆 Hessian 近似，内存为 `O(n^2)`。
- Hessian 可选，仅用于估计线搜索初始步长。

适用场景：大多数中小规模无约束优化问题，尤其是希望比 SD/CG 更快、又不想提供 Hessian 的情况。

## 收敛判据与错误处理

所有算法使用相同的收敛判据：

```text
max_i |∂F/∂x_i| <= conv_thresh
```

即梯度的无穷范数小于等于阈值时认为收敛。

默认最大迭代次数为 `1000`。若达到最大迭代次数仍未收敛，返回带有当前梯度范数的 `Err`。常见错误包括：

- 自变量为空或包含非有限值；
- 收敛阈值不是有限正数；
- 梯度维数不匹配，或牛顿类方法的 Hessian 形状不匹配；
- 函数值或梯度包含非有限值，或牛顿类方法使用的 Hessian 包含非有限值；
- 线搜索失败；
- 牛顿方向不是下降方向（纯牛顿法）。

## 使用示例

### 只需要梯度的最速下降法

```rust
use rswfn::maths::optim::{FuncEval, steepest_descent};

// F(x, y) = x^2 + 2y^2
let func = |paras: &Vec<f64>| {
    let x = paras[0];
    let y = paras[1];
    FuncEval::from_grad(x * x + 2.0 * y * y, vec![2.0 * x, 4.0 * y])
};

let result = steepest_descent(func, &vec![1.0, 1.0], 1e-6).unwrap();
println!("{:?}", result);
```

### 使用解析 Hessian 的阻尼牛顿法

```rust
use ndarray::Array2;
use rswfn::maths::optim::{FuncEval, damped_newton};

let func = |paras: &Vec<f64>| {
    let x = paras[0];
    let y = paras[1];

    let value = x * x + 2.0 * y * y;
    let grad = vec![2.0 * x, 4.0 * y];
    let hessian = Array2::from_shape_vec(
        (2, 2),
        vec![2.0, 0.0, 0.0, 4.0],
    ).unwrap();

    FuncEval::new(value, Some(grad), Some(hessian))
};

let result = damped_newton(func, &vec![1.0, 1.0], 1e-8).unwrap();
println!("{:?}", result);
```

### 只提供函数值

```rust
use rswfn::maths::optim::{FuncEval, quasi_newton};

let func = |paras: &Vec<f64>| {
    let x = paras[0];
    let y = paras[1];
    FuncEval::from_value(x * x + 2.0 * y * y)
};

// 梯度会自动用中心差分近似
let result = quasi_newton(func, &vec![1.0, 1.0], 1e-6).unwrap();
println!("{:?}", result);
```

## 算法对比

| 算法 | 所需导数 | 主要方向 | 线搜索 | 收敛特点 | 适用场景 |
| --- | --- | --- | --- | --- | --- |
| 最速下降法 | 梯度 | `-g` | Armijo | 线性，简单稳健但可能很慢 | 基线、高维问题 |
| 共轭梯度法 | 梯度 | `-g + beta p` | 强 Wolfe / Armijo | 通常快于 SD，内存低于 QN | 大规模、内存受限 |
| 牛顿法 | 梯度 + Hessian | `-H^{-1} g` | 无 | 局部二次收敛 | 初值好、Hessian 可靠 |
| 阻尼牛顿法 | 梯度 + Hessian | 修正后的 `-H^{-1} g` | Armijo | 局部快速且更稳健 | Hessian 可能非正定 |
| 拟牛顿法（BFGS） | 梯度 | `-B g` | 强 Wolfe / Armijo | 超线性，通用 | 中小规模通用优化 |

## 选择建议

- 只想快速验证问题、或需要最简单的基线：使用最速下降法。
- 问题维度较高、内存受限、只愿意计算梯度：使用共轭梯度法。
- 已有解析 Hessian，且初始点接近极小点：使用牛顿法。
- 有 Hessian，但初始点可能较远或 Hessian 可能非正定：使用阻尼牛顿法。
- 只有梯度，希望获得比 SD/CG 更快的收敛：优先使用拟牛顿法（BFGS）。
- 如果只能提供函数值，且 `O(n^2)` 次函数求值不可接受，应避免牛顿类方法，改用 CG 或 BFGS。

## 测试

优化模块包含覆盖二次函数、Rosenbrock 函数、解析梯度、数值梯度、仅 Hessian 等场景的单元测试：

```bash
cargo test --lib optim
```

## 注意事项

- `conv_thresh` 是梯度的绝对尺度阈值。若函数本身尺度很大或很小，应相应调整阈值。
- 数值梯度存在截断误差和舍入误差，阈值不宜设置得比数值噪声还小。
- 纯牛顿法是局部方法，不保证从任意初始点收敛；需要更强的全局性时请使用阻尼牛顿法或 BFGS。
- 牛顿法和阻尼牛顿法需要对 Hessian 做特征分解，计算复杂度为 `O(n^3)`；BFGS 的逆 Hessian 近似需要 `O(n^2)` 内存。
- 所有算法都假设目标函数在迭代过程中返回有限值；若函数值、梯度或 Hessian 出现 `NaN`/`Inf`，会返回错误。
