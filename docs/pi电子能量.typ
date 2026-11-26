#set page(paper: "a4", margin: (x: 2cm, y: 2cm), numbering: "1")
#set text(font: ("Libertinus Serif", "Microsoft YaHei"), size: 10.5pt, lang: "zh")
#set par(justify: true, leading: 0.85em)
#show raw.where(block: true): set text(size: 8.5pt)
#show table: set text(size: 9pt)

= π 电子能量（POCV / MOCV）的计算原理

#block(fill: luma(242), inset: 10pt, radius: 3pt, width: 100%)[
  *对应代码*：`src/props/moleprop/energy.rs`（`Calculator::pi_pocv` / `pi_mocv` / `pi_pocv_parts` / `pi_mocv_parts` / `pi_decompose`）\
  *π 投影*：`src/props/orbtprop/obtmat.rs`（`pi_pocv` / `pi_mocv` / `pocv` / `mocv`）、`src/props/orbtprop/mocv.rs`（`MocvDest::PiProperty`、`mocv_atos`）\
  *能量泛函*：`src/qmcal/ks.rs`（`rks_ele_eng`、`RKS::xc_eng`、`rsint::jk_mat`）\
  *命令行*：`src/cli.rs`（`moleprop pi-energy`、`print_pi_energy`）\
  *集成测试*：`tests/energy.rs`（C6H6）
]

== 0. 结论速览

- *不能*把"分子轨道能量求和 $sum_i n_i epsilon_i$"当作分子（或 π 电子）的能量——它重复计入了电子–电子相互作用。
- *不能*把 π 投影密度 $D_pi$ 直接代入 KS 能量泛函当作"π 电子能量"——那是"把 6 个电子放进裸核场"的非物理体系，会给出 $-70$ Ha 之类的无意义大负数。
- 正确做法：把*完整密度*的能量泛函按 $D = D_pi + D_sigma$ 做*加法分解* $E_("elec") = E_pi + E_sigma + E_(pi sigma)$（三项之和严格等于完整电子能量），并把 *π 的可加贡献*定义为 $E_pi + E_(pi sigma)$。
- 它是"π 电子对电子能量的份额"（*电子能量量级*，大负数），而不是分子总能量；要与总能量比较还需再分摊核排斥能。

== 1. 为什么不能简单相加

*轨道能量和不是能量（双重计数）。* 对 HF：

$
  E_("elec") = sum_i^"occ" n_i epsilon_i - V_("ee"),
$

其中 $V_("ee")$ 是总电子–电子排斥。$sum epsilon_i$ 让每个电子都"感受"了所有其它电子，因而重复计入 $V_("ee")$。KS-DFT 亦然（纯泛函）：

$
  E_("elec") = sum_i^"occ" n_i epsilon_i
    - 1/2 integral integral frac(rho(bold(r)) rho(bold(r)'), abs(bold(r) - bold(r)')) d bold(r) d bold(r)'
    - integral v_("xc")(bold(r)) rho(bold(r)) d bold(r) + E_("xc").
$

所以 $sum n_i epsilon_i$ 只能作为*单电子指标*（如与 Hückel π 能量类比），*不是能量*。

*把 $D_pi$ 代入完整泛函也不对。* KS 能量泛函

$
  E[D] = op("Tr")(D H) + 1/2 op("Tr")(D J[D]) - a_x/4 op("Tr")(D K[D]) + E_("xc")[rho]
$

只在 $D$ 是*整个分子*的密度、且核固定时才有意义。若代入的是 π 投影密度 $D_pi$（约 6 个电子）：

- 核吸引项 $op("Tr")(D_pi H)$ 用的是*全核、未屏蔽*（每个 C 核 $Z = 6$）；
- 只有 π–π 库仑/交换，*缺少 σ 电子对 π 的屏蔽*；
- *没有加核排斥* $E_("nuc")$。

三者叠加使 $E[D_pi] approx -70$ Ha 这种大负数出现，但并无化学意义。

== 2. 完整电子能量

原子单位制，闭壳层（占据数 2），记号：

#table(
  columns: 2,
  align: (left, left),
  inset: 5pt,
  stroke: 0.5pt + luma(170),
  [*符号*], [*含义*],
  [$D$], [密度矩阵（闭壳层 $D = 2 C_("occ") C_("occ")^top$）],
  [$H = T + V_("ne")$], [核心哈密顿量（$T$ 动能，$V_("ne")$ 电子–核吸引）],
  [$J[D]_(mu nu) = sum_(lambda sigma) (mu nu | lambda sigma) D_(lambda sigma)$], [库仑矩阵],
  [$K[D]_(mu nu) = sum_(lambda sigma) (mu lambda | nu sigma) D_(lambda sigma)$], [精确交换矩阵],
  [$a_x$], [泛函精确交换系数（B3LYP 0.2，PBE0 0.25，纯泛函 0）],
  [$E_("xc")[rho]$], [DFT 半局域交换关联能（数值格点积分）],
  [$E_("nuc") = (Z_A Z_B) / R_(A B)$], [核–核排斥能（对 $A < B$ 求和）],
)

$
  E_("elec") = op("Tr")(D H) + 1/2 op("Tr")(D J[D]) - a_x/4 op("Tr")(D K[D]) + E_("xc")[rho]
$

$
  E_("total") = E_("elec") + E_("nuc")
$

$E_("total")$ 即高斯 `SCF Done` 报告的总能量，可作为自检基准。

== 3. π / σ / πσ 分解

令 π 投影把占据轨道系数投影到 π 子空间，得到 π 投影系数 $C_pi$，进而

$
  D_pi = 2 C_pi C_pi^top, quad D_sigma = D - D_pi.
$

（$D_pi = P D P$，$P$ 为 π 子空间投影；$D_sigma$ 取*补集*以保证 $D = D_pi + D_sigma$ *精确成立*。）

由于 $J$、$K$ 对 $D$ *线性*（$J[D] = J_pi + J_sigma$，$K[D] = K_pi + K_sigma$），库仑与精确交换部分可*精确*拆分；只有 $E_("xc")$ 非线性，采用 $E_("xc")[rho_pi]$、$E_("xc")[rho_sigma]$，并把差值归入交叉项。定义

$
  E_pi = op("Tr")(D_pi H) + 1/2 op("Tr")(D_pi J_pi) - a_x/4 op("Tr")(D_pi K_pi) + E_("xc")[rho_pi]
$

$
  E_sigma = op("Tr")(D_sigma H) + 1/2 op("Tr")(D_sigma J_sigma) - a_x/4 op("Tr")(D_sigma K_sigma) + E_("xc")[rho_sigma]
$

$
  E_(pi sigma) = op("Tr")(D_pi J_sigma) - a_x/2 op("Tr")(D_pi K_sigma) + (E_("xc")[rho] - E_("xc")[rho_pi] - E_("xc")[rho_sigma])
$

其中 $E_(pi sigma)$ 的三个分项：

- *π–σ 库仑* $op("Tr")(D_pi J_sigma)$：π 与 σ 电子云之间的静电排斥，通常为大正；
- *π–σ 交换* $-a_x/2 op("Tr")(D_pi K_sigma)$：通常为小负；
- *XC 残差* $E_("xc")[rho] - E_("xc")[rho_pi] - E_("xc")[rho_sigma]$：通常为小负。

*可加性证明（草图）*：把 $1/2 op("Tr")(D J[D])$ 展开，

$
  1/2 op("Tr")(D J[D]) = 1/2 op("Tr")(D_pi J_pi) + 1/2 op("Tr")(D_sigma J_sigma) + op("Tr")(D_pi J_sigma),
$

精确交换项同理（交叉交换系数为 $a_x / 2$），XC 残差按定义补齐，于是

$
  E_pi + E_sigma + E_(pi sigma) = E_("elec")
$

严格成立。*π 的可加贡献*取 $E_pi + E_(pi sigma)$。

== 4. π 投影方法：POCV 与 MOCV

两者都产出"π 部分的占据轨道系数"，差别在投影方式：

*POCV（P 轨道垂直方向投影法）*

- 对每个（重）原子求一个垂直于分子平面的 π 方向 $hat(bold(n))$（`props::atomprop::direction`）；
- 只保留共轭原子上的 p 轨道分量，并把该原子的 $(p_x, p_y, p_z)$ 投影到 $hat(bold(n))$ 上（`maths::vec_project`），其余分量置零；
- `pocv(&dirs, false, false)`：对占据轨道逐列投影；
- 返回系数为 *raw* 形式。

*MOCV（分子轨道系数向量）*

- 为每个原子构造局部坐标系（LCS，`Stm`）；把分子轨道系数变换到各自 LCS；
- 由 `mocv_atos(mole, dest)` 选出 π 对称的轨道（如 `PX`、`D+1`、`D-2`、`F+1` / `F-2` / `F+3` 等，`MocvDest::PiProperty`）；
- 只保留这些 π 型 AO 的系数；
- 返回系数为 *sph* 形式。

无论哪种，都在 `pi_decompose` 里通过 `Coefs::set_cmat` 写入并按 `basis().int_form()`（`car` 或 `sph`）转到积分形式，再构造成 $D_pi$。

== 5. 实现流程

`Calculator` 提供 `pi_pocv` / `pi_mocv`（返回 π 贡献 $E_pi + E_(pi sigma)$）与 `pi_pocv_parts` / `pi_mocv_parts`（返回完整分解 `PiEnergy`）。

+ *准备*（`prepare_log`）：仅接受 log/out 文件；用 `LogReader::read_gau_job()` 读出泛函与基组名，泛函经 `XcFunc::from_name` 识别；要求*闭壳层*且含轨道系数；基组缺失时按名字用 `rsbse` 加载。
+ *π 系数*：`obtmat::Calculator::pi_pocv()`（raw）或 `pi_mocv(&MocvDest::PiProperty)`（sph）。
+ *分解*（`pi_decompose`）：
  - 完整密度 $D = 2 C_("occ") C_("occ")^top$；π 密度 $D_pi = 2 C_pi C_pi^top$；$D_sigma = D - D_pi$；
  - `RsInt::mat1e` 得 $H = T + V$，`RsInt::mat2e` 得双电子积分；
  - `rsint::jk_mat(D_pi)` 得 $(J_pi, K_pi)$、`jk_mat(D_sigma)` 得 $(J_sigma, K_sigma)$，相加得 $J_("tot"), K_("tot")$；
  - `RKS::new(mole, func)` 在数值格点上求 $E_("xc")[rho]$、$E_("xc")[rho_pi]$、$E_("xc")[rho_sigma]$，以及 $E_("nuc")$；
  - 按第 3 节公式计算 $E_pi, E_sigma, E_(pi sigma), E_("elec")$。
+ *输出*（`PiEnergy`）：

```rust
pub struct PiEnergy {
    pub e_pi: f64,         // E_pi
    pub e_sigma: f64,      // E_sigma
    pub e_cross: f64,      // E_pi_sigma（= 库仑 + 交换 + XC 残差）
    pub e_cross_coul: f64, // Tr(D_pi J_sigma)
    pub e_cross_exch: f64, // -a_x/2·Tr(D_pi K_sigma)
    pub e_cross_xc: f64,   // E_xc[rho]-E_xc[rho_pi]-E_xc[rho_sigma]
    pub e_elec: f64,       // E_elec
    pub e_nuc: f64,        // E_nuc
    pub e_total: f64,      // E_elec + E_nuc
}
// PiEnergy::pi() == e_pi + e_cross
```

== 6. 数值示例：C6H6（B3LYP/6-31G(d)）

单位 Hartree。log 中 `E(RB3LYP) = -232.247385084`。

#table(
  columns: 8,
  align: center,
  inset: 4pt,
  stroke: 0.5pt + luma(170),
  [*方法*], [$E_pi$], [$E_sigma$], [$E_(pi sigma)$], [$E_("elec")$], [$E_("nuc")$], [$E_("total")$], [*π 贡献*],
  [POCV], [$-70.3786$], [$-428.8168$], [$+63.8634$], [$-435.3321$], [$+203.0815$], [*$-232.2506$*], [*$-6.5153$*],
  [MOCV], [$-71.5165$], [$-428.6089$], [$+64.7934$], [$-435.3321$], [$+203.0815$], [*$-232.2506$*], [*$-6.7232$*],
)

要点：

- $E_pi + E_sigma + E_(pi sigma) = E_("elec")$ *精确闭合*（测试断言误差 $< 10^(-6)$ Ha）。
- $E_("total")$ 与 log 的 `SCF Done` 相差仅约 *3 mHa*（来自数值格点与积分/泛函微差），且与投影方法无关（POCV、MOCV 的 $E_("elec")$ / $E_("total")$ 相同）——这是正确性的强验证。
- π 贡献约 *$-6.5$ Ha*（6 个 π 电子，约 1 Ha/π 电子量级）。

=== 6.1 为什么 $E_(pi sigma)$ 为正

$E_(pi sigma)$ 的三个分项（C6H6）：

#table(
  columns: 5,
  align: center,
  inset: 4pt,
  stroke: 0.5pt + luma(170),
  [*方法*], [*库仑* $op("Tr")(D_pi J_sigma)$], [*交换* $-a_x/2 op("Tr")(D_pi K_sigma)$], [*XC 残差*], [*合计*],
  [POCV], [*$+64.888$*], [$-0.193$], [$-0.831$], [$+63.863$],
  [MOCV], [*$+65.819$*], [$-0.189$], [$-0.837$], [$+64.793$],
)

- *主导项是 π–σ 库仑排斥*：$op("Tr")(D_pi J_sigma) = sum D_(pi, mu nu) D_(sigma, lambda sigma) (mu nu | lambda sigma)$，二电子积分与两个密度矩阵皆正，两团负电荷的静电排斥恒为正；且 π（6 电子）与 σ（36 电子）全程 $1/r_12$ 相互作用，量级最大。
- *交换*：交换积分正，但系数 $-a_x/2$（B3LYP $a_x = 0.2$）带负号，为小负。
- *XC 残差*：$x^(4/3)$ 超加性，$E_("xc")[rho] < E_("xc")[rho_pi] + E_("xc")[rho_sigma]$，为小负（稳定化）。

因此净值为正。物理解读：$E_pi$ 把 π 电子当成"裸核场里只有它自己"（过负），$E_(pi sigma)$ 中的 π–σ 库仑排斥把它抵回来，二者之和才是 π 的净贡献。*单独看 $E_pi$ 或 $E_(pi sigma)$ 都会误导，只有和（及总能量）有物理意义。*

== 7. 约定与注意事项

- *适用范围*：仅 log/out 文件、闭壳层、且泛函可被 `XcFunc` 识别（HF 或 `m06l` 等暂不支持）。
- *XC 不可分离*：$E_("xc")[rho] != E_("xc")[rho_pi] + E_("xc")[rho_sigma]$，本实现用 $E_("xc")[rho_pi]$、$E_("xc")[rho_sigma]$ 并将差值归入 $E_(pi sigma)$。$E_pi$、$E_(pi sigma)$ *各自的数值依赖此约定*；稳健的是二者之和与 $E_("elec")$。
- *量纲*：π 贡献是"电子能量份额"（大负数），不是分子总能量；要与总能量同量级比较，需再按约定分摊核排斥能 $E_("nuc")$。
- *$D_sigma$ 可能局部为负*（$D_sigma = D - D_pi$）：当 π 密度接近总密度时，个别格点 $rho_sigma$ 可能略小于 0；常规共轭体系（如苯）不会触发。若遇异常体系需留意。
- *性能*：每个方法需构建一次 100×74 数值格点并做 3 次 $E_("xc")$ 积分（$rho, rho_pi, rho_sigma$）与 2 次 `jk_mat`，C6H6 集成测试约 55 s。

== 8. 命令行与库调用

=== 8.1 命令行

```bash
# π 电子能量分解（默认 POCV）
rswfn moleprop pi-energy moles/C6H6.out

# 改用 MOCV
rswfn moleprop pi-energy moles/C6H6.out --method mocv
```

输出包含 π/σ/πσ 三项、交叉项的库仑/精确交换/XC 残差分解、完整电子能量、
核排斥能、总能量，以及 π 的可加贡献 $E_pi + E_(pi sigma)$。

=== 8.2 库调用

```rust
use rswfn::base::Mole;
use rswfn::props::moleprop::energy::Calculator;

let mole = Mole::from_file("moles/C6H6.out");
let caler = Calculator::new(&mole);

// 只要 π 的可加贡献 E_π + E_πσ
let e_pi = caler.pi_pocv().unwrap();

// 完整分解
let parts = caler.pi_pocv_parts().unwrap();
println!("E_π = {}", parts.e_pi);
println!("E_πσ = {}", parts.e_cross);
println!("E_elec = {}", parts.e_elec);
```

`--method mocv` 对应 `pi_mocv` / `pi_mocv_parts`。
