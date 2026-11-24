// uks_x_lda
// rks_x_lda
// uks_x_b88
// rks_x_b88
// uks_x_pbe
// rks_x_pbe
// uks_x_pw91
// rks_x_pw91
// uks_c_vwn5
// rks_c_vwn5
// uks_c_p86
// rks_c_p86
// uks_c_lyp
// rks_c_lyp
// uks_c_pw91
// rks_c_pw91
// uks_c_pbe
// rks_c_pbe
// uks_xc_b97
// rks_xc_b97
// uks_xc_hcth407
// rks_xc_hcth407

## 输入参数

### 控制参数

- `ideriv`：导数阶数控制
  - `0`：只计算能量密度 `zk`
  - `1`：计算能量密度 + 一阶导数（势能）
  - `2`：计算能量密度 + 一阶导数 + 二阶导数

- `npt`：网格点的数量（数组长度）

### 密度变量（输入）

- `rhoa1(npt)`：α 自旋电子密度（每个网格点）
- `rhob1(npt)`：β 自旋电子密度（每个网格点）
- `sigmaaa1(npt)`：∇ρ*α · ∇ρ*α（α 自旋密度梯度）
- `sigmabb1(npt)`：∇ρ*β · ∇ρ*β（β 自旋密度梯度）
- `sigmaab1(npt)`：∇ρ*α · ∇ρ*β（自旋交叉密度梯度）

## 输出参数

### 能量密度（输出）

- `zk(npt)`：交换能密度（每个网格点的能量贡献）

### 一阶导数/势能（输出，ideriv ≥ 1）

- `vrhoa(npt)`：∂ε/∂ρ_α（α 自旋的交换势）
- `vrhob(npt)`：∂ε/∂ρ_β（β 自旋的交换势）
- `vsigmaaa(npt)`：∂ε/∂σ_αα
- `vsigmabb(npt)`：∂ε/∂σ_ββ
- `vsigmaab(npt)`：∂ε/∂σ_αβ

### 二阶导数（输出，ideriv ≥ 2）

- `v2rhoa2(npt)`：∂²ε/∂ρ_α²
- `v2rhob2(npt)`：∂²ε/∂ρ_β²
- `v2rhoab(npt)`：∂²ε/∂ρ*α∂ρ*β
- `v2rhoasigmaaa(npt)`：∂²ε/∂ρ*α∂σ*αα
- `v2rhoasigmaab(npt)`：∂²ε/∂ρ*α∂σ*αβ
- `v2rhoasigmabb(npt)`：∂²ε/∂ρ*α∂σ*ββ
- `v2rhobsigmabb(npt)`：∂²ε/∂ρ*β∂σ*ββ
- `v2rhobsigmaab(npt)`：∂²ε/∂ρ*β∂σ*αβ
- `v2rhobsigmaaa(npt)`：∂²ε/∂ρ*β∂σ*αα
- `v2sigmaaa2(npt)`：∂²ε/∂σ_αα²
- `v2sigmaaaab(npt)`：∂²ε/∂σ*αα∂σ*αβ
- `v2sigmaaabb(npt)`：∂²ε/∂σ*αα∂σ*ββ
- `v2sigmaab2(npt)`：∂²ε/∂σ_αβ²
- `v2sigmaabbb(npt)`：∂²ε/∂σ*αβ∂σ*ββ
- `v2sigmabb2(npt)`：∂²ε/∂σ_ββ²

## 物理含义

该函数实现的是 **Dirac 交换能泛函**（LDA 交换）：

```
ε_x^LDA(ρ_α, ρ_β) = -3/4 * (3/π)^(1/3) * (ρ_α^(4/3) + ρ_β^(4/3))
```

常数 `0.9305257363491` 实际对应 `3/4 * (3/π)^(1/3)`。

对于闭壳层系统（ρ*α = ρ*β = ρ/2），它简化为常用的 LDA 交换泛函形式：

```
ε_x = -3/4 * (3/π)^(1/3) * ρ^(4/3)
```

**注意**：虽然该函数接收了 sigma（梯度）相关的数组参数，但在纯 LDA 泛函中这些参数未被使用（导数为零）。这些参数的存在是为了与 GGA 泛函保持统一的接口规范。
