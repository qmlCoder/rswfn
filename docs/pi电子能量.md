# π 电子能量（POCV / MOCV）的计算原理

> 对应代码：`src/props/moleprop/energy.rs`
> （`Calculator::pi_pocv` / `pi_mocv` / `pi_pocv_parts` / `pi_mocv_parts` / `pi_decompose`）
> π 投影：`src/props/orbtprop/obtmat.rs`（`pi_pocv` / `pi_mocv` / `pocv` / `mocv`）、
> `src/props/orbtprop/mocv.rs`（`MocvDest::PiProperty`、`mocv_atos`）
> 能量泛函：`src/qmcal/ks.rs`（`rks_ele_eng`、`RKS::xc_eng`、`rsint::jk_mat`）
> 集成测试：`tests/energy.rs`（C6H6）

## 0. 结论速览

- **不能**把"分子轨道能量求和 $\sum_i n_i\varepsilon_i$"当作分子（或 π 电子）的能量——它重复计入了电子–电子相互作用。
- **不能**把 π 投影密度 $D_\pi$ 直接代入 KS 能量泛函当作"π 电子能量"——那是"把 6 个电子放进裸核场"的非物理体系，会给出 −70 Ha 之类的无意义大负数。
- 正确做法：把**完整密度**的能量泛函按 $D=D_\pi+D_\sigma$ 做**加法分解**
  $E_{\text{elec}}=E_\pi+E_\sigma+E_{\pi\sigma}$（三项之和严格等于完整电子能量），
  并把 **π 的可加贡献**定义为

$$
E_{\pi}^{\text{贡献}} = E_\pi + E_{\pi\sigma}.
$$

  它是"π 电子对电子能量的份额"（**电子能量量级**，大负数），而不是分子总能量；要与总能量比较还需再分摊核排斥能。

## 1. 为什么不能简单相加

**轨道能量和不是能量（双重计数）。** 对 HF：

$$
E_{\text{elec}} = \sum_i^{\text{occ}} n_i\varepsilon_i - V_{ee},
$$

其中 $V_{ee}$ 是总电子–电子排斥。$\sum \varepsilon_i$ 让每个电子都"感受"了所有其它电子，因而重复计入 $V_{ee}$。
KS-DFT 亦然（纯泛函）：

$$
E_{\text{elec}} = \sum_i^{\text{occ}} n_i\varepsilon_i
-\frac12\iint\frac{\rho(\mathbf r)\rho(\mathbf r')}{|\mathbf r-\mathbf r'|}\,d\mathbf r\,d\mathbf r'
-\int v_{xc}(\mathbf r)\rho(\mathbf r)\,d\mathbf r + E_{xc}.
$$

所以 $\sum n_i\varepsilon_i$ 只能作为**单电子指标**（如与 Hückel π 能量类比），**不是能量**。

**把 $D_\pi$ 代入完整泛函也不对。** KS 能量泛函

$$
E[D]=\operatorname{Tr}(DH)+\tfrac12\operatorname{Tr}\big(D\,J[D]\big)-\frac{a_x}{4}\operatorname{Tr}\big(D\,K[D]\big)+E_{xc}[\rho]
$$

只在 $D$ 是**整个分子**的密度、且核固定时才有意义。若代入的是 π 投影密度 $D_\pi$（约 6 个电子）：

- 核吸引项 $\operatorname{Tr}(D_\pi H)$ 用的是**全核、未屏蔽**（每个 C 核 $Z=6$）；
- 只有 π–π 库仑/交换，**缺少 σ 电子对 π 的屏蔽**；
- **没有加核排斥** $E_{\text{nuc}}$。

三者叠加使 $E[D_\pi]\approx -70$ Ha 这种大负数出现，但并无化学意义。

## 2. 完整电子能量

原子单位制，闭壳层（占据数 2），记号：

| 符号 | 含义 |
| --- | --- |
| $D$ | 密度矩阵（闭壳层 $D=2\,C_{\text{occ}}C_{\text{occ}}^\top$） |
| $H=T+V_{\text{ne}}$ | 核心哈密顿量（$T$ 动能，$V_{\text{ne}}$ 电子–核吸引） |
| $J[D]_{\mu\nu}=\sum_{\lambda\sigma}(\mu\nu\|\lambda\sigma)D_{\lambda\sigma}$ | 库仑矩阵 |
| $K[D]_{\mu\nu}=\sum_{\lambda\sigma}(\mu\lambda\|\nu\sigma)D_{\lambda\sigma}$ | 精确交换矩阵 |
| $a_x$ | 泛函的精确交换系数（B3LYP 为 0.2，PBE0 为 0.25，纯泛函为 0） |
| $E_{xc}[\rho]$ | DFT 半局域交换关联能（数值格点积分） |
| $E_{\text{nuc}}=\sum_{A<B}Z_AZ_B/R_{AB}$ | 核–核排斥能 |

$$
\boxed{\,E_{\text{elec}}=\operatorname{Tr}(DH)+\tfrac12\operatorname{Tr}(DJ[D])
-\frac{a_x}{4}\operatorname{Tr}(DK[D])+E_{xc}[\rho]\,}
$$
$$
E_{\text{total}}=E_{\text{elec}}+E_{\text{nuc}}
$$

$E_{\text{total}}$ 即高斯 `SCF Done` 报告的总能量，可作为自检基准。

## 3. π / σ / πσ 分解

令 π 投影把占据轨道系数投影到 π 子空间，得到 π 投影系数 $C_\pi$，进而

$$
D_\pi = 2\,C_\pi C_\pi^\top,\qquad D_\sigma := D - D_\pi .
$$

（$D_\pi=PDP$，$P$ 为 π 子空间投影；$D_\sigma$ 取**补集**以保证 $D=D_\pi+D_\sigma$ **精确成立**。）

由于 $J$、$K$ 对 $D$ **线性**（$J[D]=J_\pi+J_\sigma$，$K[D]=K_\pi+K_\sigma$），
库仑与精确交换部分可**精确**拆分；只有 $E_{xc}$ 非线性，采用 $E_{xc}[\rho_\pi]$、$E_{xc}[\rho_\sigma]$，并把差值归入交叉项。定义

$$
\begin{aligned}
E_\pi &= \operatorname{Tr}(D_\pi H)+\tfrac12\operatorname{Tr}(D_\pi J_\pi)-\frac{a_x}{4}\operatorname{Tr}(D_\pi K_\pi)+E_{xc}[\rho_\pi],\\
E_\sigma &= \operatorname{Tr}(D_\sigma H)+\tfrac12\operatorname{Tr}(D_\sigma J_\sigma)-\frac{a_x}{4}\operatorname{Tr}(D_\sigma K_\sigma)+E_{xc}[\rho_\sigma],\\
E_{\pi\sigma} &= \underbrace{\operatorname{Tr}(D_\pi J_\sigma)}_{\text{π–σ 库仑}}
-\underbrace{\frac{a_x}{2}\operatorname{Tr}(D_\pi K_\sigma)}_{\text{π–σ 交换}}
+\underbrace{\big(E_{xc}[\rho]-E_{xc}[\rho_\pi]-E_{xc}[\rho_\sigma]\big)}_{\text{XC 残差}}.
\end{aligned}
$$

**可加性证明（草图）**：把 $\tfrac12\operatorname{Tr}(DJ[D])$ 展开，

$$
\tfrac12\operatorname{Tr}(DJ[D])
=\tfrac12\operatorname{Tr}(D_\pi J_\pi)+\tfrac12\operatorname{Tr}(D_\sigma J_\sigma)+\operatorname{Tr}(D_\pi J_\sigma),
$$

精确交换项同理（交叉交换系数为 $a_x/2$），XC 残差按定义补齐，于是

$$
\boxed{\,E_\pi+E_\sigma+E_{\pi\sigma}=E_{\text{elec}}\,}
$$

严格成立。**π 的可加贡献**取

$$
E_{\pi}^{\text{贡献}}=E_\pi+E_{\pi\sigma}.
$$

## 4. π 投影方法：POCV 与 MOCV

两者都产出"π 部分的占据轨道系数"，差别在投影方式：

**POCV（P 轨道垂直方向投影法）**
- 对每个（重）原子求一个垂直于分子平面的 π 方向 $\hat{\mathbf n}$（`props::atomprop::direction`）；
- 只保留共轭原子上的 p 轨道分量，并把该原子的 $(p_x,p_y,p_z)$ 投影到 $\hat{\mathbf n}$ 上（`maths::vec_project`），其余分量置零；
- `pocv(&dirs, false, false)`：对占据轨道逐列投影。
- 返回系数为 **raw** 形式。

**MOCV（Molecular Orbital Coefficient Vector）**
- 为每个原子构造局部坐标系（LCS，`Stm`）；把分子轨道系数变换到各自 LCS；
- 由 `mocv_atos(mole, dest)` 选出 π 对称的轨道（如 `PX`、`D+1`、`D-2`、`F+1`/`F-2`/`F+3`…，`MocvDest::PiProperty`）；
- 只保留这些 π 型 AO 的系数。
- 返回系数为 **sph** 形式。

无论哪种，都在 `pi_decompose` 里通过 `Coefs::set_cmat` 写入并按 `basis().int_form()`（`car` 或 `sph`）转到积分形式，再构造成 $D_\pi$。

## 5. 实现流程

`Calculator` 提供 `pi_pocv` / `pi_mocv`（返回 π 贡献 $E_\pi+E_{\pi\sigma}$）与
`pi_pocv_parts` / `pi_mocv_parts`（返回完整分解 `PiEnergy`）。

1. **准备**（`prepare_log`）：仅接受 log/out 文件；用 `LogReader::read_gau_job()` 读出泛函与基组名，
   泛函经 `XcFunc::from_name` 识别；要求**闭壳层**且含轨道系数；基组缺失时按名字用 `rsbse` 加载。
2. **π 系数**：`obtmat::Calculator::pi_pocv()`（raw）或 `pi_mocv(&MocvDest::PiProperty)`（sph）。
3. **分解**（`pi_decompose`）：
   - 完整密度 $D=2C_{\text{occ}}C_{\text{occ}}^\top$；π 密度 $D_\pi=2C_\pi C_\pi^\top$；$D_\sigma=D-D_\pi$；
   - `RsInt::mat1e` 得 $H=T+V$，`RsInt::mat2e` 得双电子积分；
   - `rsint::jk_mat(D_π)`→$(J_\pi,K_\pi)$、`jk_mat(D_σ)`→$(J_\sigma,K_\sigma)$，相加得 $J_{\text{tot}},K_{\text{tot}}$；
   - `RKS::new(mole, func)` 在数值格点上求 $E_{xc}[\rho]$、$E_{xc}[\rho_\pi]$、$E_{xc}[\rho_\sigma]$，以及 $E_{\text{nuc}}$；
   - 按第 3 节公式计算 $E_\pi,E_\sigma,E_{\pi\sigma},E_{\text{elec}}$。
4. **输出**（`PiEnergy`）：

```rust
pub struct PiEnergy {
    pub e_pi: f64,         // E_π
    pub e_sigma: f64,      // E_σ
    pub e_cross: f64,      // E_πσ（= 库仑 + 交换 + XC 残差）
    pub e_cross_coul: f64, // Tr(D_π J_σ)
    pub e_cross_exch: f64, // -a_x/2·Tr(D_π K_σ)
    pub e_cross_xc: f64,   // E_xc[ρ]-E_xc[ρ_π]-E_xc[ρ_σ]
    pub e_elec: f64,       // E_elec
    pub e_nuc: f64,        // E_nuc
    pub e_total: f64,      // E_elec + E_nuc
}
// PiEnergy::pi() == e_pi + e_cross
```

## 6. 数值示例：C6H6（B3LYP/6-31G(d)）

单位 Hartree。log 中 `SCF Done E(RB3LYP) = -232.247385084`。

| 方法 | $E_\pi$ | $E_\sigma$ | $E_{\pi\sigma}$ | $E_{\text{elec}}$ | $E_{\text{nuc}}$ | $E_{\text{total}}$ | π 贡献 $E_\pi+E_{\pi\sigma}$ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| POCV | −70.3786 | −428.8168 | +63.8634 | −435.3321 | +203.0815 | **−232.2506** | **−6.5153** |
| MOCV | −71.5165 | −428.6089 | +64.7934 | −435.3321 | +203.0815 | **−232.2506** | **−6.7232** |

要点：

- $E_\pi+E_\sigma+E_{\pi\sigma}=E_{\text{elec}}$ **精确闭合**（测试断言误差 $<10^{-6}$ Ha）。
- $E_{\text{total}}$ 与 log 的 `SCF Done` 相差仅约 **3 mHa**（来自数值格点与积分/泛函微差），且与投影方法无关（POCV、MOCV 的 $E_{\text{elec}}/E_{\text{total}}$ 相同）——这是正确性的强验证。
- π 贡献约 **−6.5 Ha**（6 个 π 电子，约 1 Ha/π 电子量级）。

### 6.1 为什么 $E_{\pi\sigma}$ 为正

$E_{\pi\sigma}$ 的三个分项（C6H6）：

| 方法 | 库仑 $\operatorname{Tr}(D_\pi J_\sigma)$ | 交换 $-\tfrac{a_x}{2}\operatorname{Tr}(D_\pi K_\sigma)$ | XC 残差 | 合计 |
| --- | --- | --- | --- | --- |
| POCV | **+64.888** | −0.193 | −0.831 | +63.863 |
| MOCV | **+65.819** | −0.189 | −0.837 | +64.793 |

- **主导项是 π–σ 库仑排斥**：$\operatorname{Tr}(D_\pi J_\sigma)=\sum D_{\pi,\mu\nu}D_{\sigma,\lambda\sigma}(\mu\nu\|\lambda\sigma)$，二电子积分与两个密度矩阵皆正 → 两团负电荷的静电排斥，恒为正，且 π（6 电子）与 σ（36 电子）全程 $1/r_{12}$ 相互作用，量级最大。
- **交换**：交换积分正，但系数 $-\tfrac{a_x}{2}$（B3LYP $a_x=0.2$）带负号 → 小负。
- **XC 残差**：$x^{4/3}$ 超加性 $\Rightarrow E_{xc}[\rho]<E_{xc}[\rho_\pi]+E_{xc}[\rho_\sigma]$ → 小负（稳定化）。

因此净值为正。物理解读：$E_\pi$ 把 π 电子当成"裸核场里只有它自己"（过负），$E_{\pi\sigma}$ 中的 π–σ 库仑排斥把它抵回来，二者之和才是 π 的净贡献。**单独看 $E_\pi$ 或 $E_{\pi\sigma}$ 都会误导，只有和（及总能量）有物理意义。**

## 7. 约定与注意事项

- **适用范围**：仅 log/out 文件、闭壳层、且泛函可被 `XcFunc` 识别（HF 或 `m06l` 等暂不支持）。
- **XC 不可分离**：$E_{xc}[\rho]\neq E_{xc}[\rho_\pi]+E_{xc}[\rho_\sigma]$，本实现用 $E_{xc}[\rho_\pi]$、$E_{xc}[\rho_\sigma]$ 并将差值归入 $E_{\pi\sigma}$。$E_\pi$、$E_{\pi\sigma}$ **各自的数值依赖此约定**；稳健的是二者之和与 $E_{\text{elec}}$。
- **量纲**：π 贡献是"电子能量份额"（大负数），不是分子总能量；要与总能量同量级比较，需再按约定分摊核排斥能 $E_{\text{nuc}}$。
- **$D_\sigma$ 可能局部为负**（$D_\sigma=D-D_\pi$）：当 π 密度接近总密度时，个别格点 $\rho_\sigma$ 可能略小于 0；常规共轭体系（如苯）不会触发。若遇异常体系需留意。
- **性能**：每个方法需构建一次 100×74 数值格点并做 3 次 $E_{xc}$ 积分（$\rho,\rho_\pi,\rho_\sigma$）与 2 次 `jk_mat`，C6H6 集成测试约 55 s。
