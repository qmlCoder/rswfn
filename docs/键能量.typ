#set page(paper: "a4", margin: (x: 2cm, y: 2cm), numbering: "1")
#set text(font: ("Libertinus Serif", "Microsoft YaHei"), size: 10.5pt, lang: "zh")
#set par(justify: true, leading: 0.85em)
#show raw.where(block: true): set text(size: 8.5pt)
#show table: set text(size: 9pt)

= 键能量（bond / rest / cross）的计算原理

#block(fill: luma(242), inset: 10pt, radius: 3pt, width: 100%)[
  *对应代码*：`src/props/bondprop/energy.rs`（接口：`Calculator::bond_energy` / `bond_parts`、`BondEnergy`）、`src/props/bondprop/energy_v1.rs`（实现：`bond_parts` / `decompose`）\
  *公共准备*：`src/props/moleprop/energy.rs`（`prepare_log`，读取 log 中的泛函与基组）\
  *能量泛函*：`src/qmcal/ks.rs`（`RKS::xc_eng`、`rsint::jk_mat`）\
  *命令行*：`src/cli.rs`（`bondprop energy`、`print_bond_energy`）\
  *集成测试*：`tests/energy.rs`（`test_bond_energy_c6h6`）
]

== 0. 结论速览

- *键密度* $D_("bond")$：只保留两个原子的原子轨道（AO）系数得到的密度块，即 $D$ 在 $A union B$ 上的块 $P_(A B) D P_(A B)$。
- *不能*把 $D_("bond")$ 单独代入完整 KS 泛函当作"键能"——它缺少其余电子的屏蔽，核吸引却仍用*全部*原子核，必然给出无意义的大负数。
- 正确做法：把*完整密度*按 $D = D_("bond") + D_("rest")$ 做*加法分解*
  $E_("elec") = E_("bond") + E_("rest") + E_("cross")$（三项之和严格等于完整电子能量）。
- *键的可加贡献*定义为 $E_("bond") + E_("cross")$，其余部分为 $E_("rest")$，二者之和恰为 $E_("elec")$。
- 它是"该两原子密度块对电子能量的份额"（*电子能量量级*，大负数），*不是键解离能*。

== 1. 为什么不能把键密度直接代入完整泛函

KS 能量泛函

$
  E[D] = op("Tr")(D H) + 1/2 op("Tr")(D J[D]) - a_x/4 op("Tr")(D K[D]) + E_("xc")[rho]
$

只在 $D$ 是*整个分子*的密度、且核固定时才有意义。若只保留两个原子的系数得到 $D_("bond")$ 并直接代入：

- 核吸引项 $op("Tr")(D_("bond") H)$ 使用的是*全部分子核、未屏蔽*（例如每个 C 核 $Z = 6$）；
- 只有键内电子的库仑/交换，*缺少其余电子对键的屏蔽*；
- 若再叠加*整个分子*的核排斥 $E_("nuc")$，各项量级完全不匹配。

其结果是一个"把小密度放进全部裸核场"的非物理能量，数值可达数十 Hartree 的大负数，没有化学意义。这与 π 能量分解中"不能把 $D_pi$ 单独代入泛函"是同一个坑。

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

== 3. bond / rest / cross 分解

设选定的原子对为 $(A, B)$，令 $P_(A B)$ 为"只保留 $A$、$B$ 的原子轨道、其余原子轨道清零"的对角投影矩阵。保留两个原子的占据轨道系数 $C_("bond") = P_(A B) C_("occ")$，得到

$
  D_("bond") = 2 C_("bond") C_("bond")^top = P_(A B) D P_(A B),
$

即 $D$ 在 $A union B$ 上的*密度块*。其余部分取补密度

$
  D_("rest") := D - D_("bond"),
$

于是 $D = D_("bond") + D_("rest")$ *精确成立*。（注意 $D_("rest")$ 并非"其余原子的密度块"，它还包含键与其余部分之间的交叉块。）

由于 $J$、$K$ 对 $D$ *线性*，库仑与精确交换部分可*精确*拆分；只有 $E_("xc")$ 非线性，采用 $E_("xc")[rho_("bond")]$、$E_("xc")[rho_("rest")]$，并把差值归入交叉项。定义

$
  E_("bond") = op("Tr")(D_("bond") H) + 1/2 op("Tr")(D_("bond") J_("bond")) - a_x/4 op("Tr")(D_("bond") K_("bond")) + E_("xc")[rho_("bond")]
$

$
  E_("rest") = op("Tr")(D_("rest") H) + 1/2 op("Tr")(D_("rest") J_("rest")) - a_x/4 op("Tr")(D_("rest") K_("rest")) + E_("xc")[rho_("rest")]
$

$
  E_("cross") = op("Tr")(D_("bond") J_("rest")) - a_x/2 op("Tr")(D_("bond") K_("rest")) + (E_("xc")[rho] - E_("xc")[rho_("bond")] - E_("xc")[rho_("rest")])
$

其中 $E_("cross")$ 的三个分项：

- *键–其余库仑* $op("Tr")(D_("bond") J_("rest"))$：键电子与其余电子的静电排斥，通常为大正；
- *键–其余交换* $-a_x/2 op("Tr")(D_("bond") K_("rest"))$：通常为小负；
- *XC 残差* $E_("xc")[rho] - E_("xc")[rho_("bond")] - E_("xc")[rho_("rest")]$：通常为小负。

*可加性证明（草图）*：把 $1/2 op("Tr")(D J[D])$ 按 $D = D_("bond") + D_("rest")$ 展开，

$
  1/2 op("Tr")(D J[D]) = 1/2 op("Tr")(D_("bond") J_("bond")) + 1/2 op("Tr")(D_("rest") J_("rest")) + op("Tr")(D_("bond") J_("rest")),
$

精确交换项同理（交叉交换系数为 $a_x/2$），XC 残差按定义补齐，于是

$
  E_("bond") + E_("rest") + E_("cross") = E_("elec")
$

严格成立。*键的可加贡献*取

$
  E_("bond")^"贡献" = E_("bond") + E_("cross"), quad E_("rest")^"贡献" = E_("rest"),
$

二者之和恰为 $E_("elec")$。

== 4. 键密度的构造

1. 取积分形式 `form = mole.basis().int_form()`（`"car"` 或 `"sph"`），占据轨道索引 `occ = coefs().obt_idxs(true)`；
2. 由 `coefs().get_cmat(form, mole)` 取该形式的系数矩阵，选出占据列 $C_("occ")$，得完整密度 $D = 2 C_("occ") C_("occ")^top$；
3. 由 `basis().atos().atms(form)` 得到每行 AO 所属原子；把*不属于 $A$、$B$* 的行清零，得到 $C_("bond")$；
4. $D_("bond") = 2 C_("bond") C_("bond")^top$，$D_("rest") = D - D_("bond")$。

两点说明：

- $D_("bond")$ 就是 $D$ 的 $(A union B)$ 主块（含 $A$–$B$ 交叉块），它是半正定的（$2 (P C)(P C)^top$）；而 $D_("rest")$ 一般*不是*半正定矩阵。
- 这里*不做新的 SCF*，也不要求保留下来的系数归一或正交；被保留的是"密度块"，不是一个独立归一化的波函数。

== 5. 实现流程

`Calculator` 提供 `bond_energy(atm1, atm2)`（返回键贡献 $E_("bond") + E_("cross")$）与 `bond_parts(atm1, atm2)`（返回完整分解 `BondEnergy`），原子序号从 0 开始。

+ *准备*（复用 `moleprop::energy::Calculator::prepare_log`）：仅接受 log/out 文件；用 `LogReader::read_gau_job()` 读出泛函与基组名，泛函经 `XcFunc::from_name` 识别；要求*闭壳层*且含轨道系数；基组缺失时按名字用 `rsbse` 加载。
+ *校验*：原子序号必须落在 $[0, N_"atm")$ 内，且两原子不同。
+ *分解*（`energy_v1::decompose`）：
  - 按第 4 节构造 $D$、$D_("bond")$、$D_("rest")$；
  - `RsInt::mat1e` 得 $H = T + V$，`RsInt::mat2e` 得双电子积分；
  - `rsint::jk_mat(D_("bond"))` 得 $(J_("bond"), K_("bond"))$、`jk_mat(D_("rest"))` 得 $(J_("rest"), K_("rest"))$，相加得 $J_("tot"), K_("tot")$；
  - `RKS::new(mole, func)` 在数值格点上求 $E_("xc")[rho]$、$E_("xc")[rho_("bond")]$、$E_("xc")[rho_("rest")]$，以及 $E_("nuc")$；
  - 按第 3 节公式计算 $E_("bond"), E_("rest"), E_("cross"), E_("elec")$。
+ *输出*（`BondEnergy`）：

```rust
pub struct BondEnergy {
    pub e_bond: f64,         // E_bond
    pub e_rest: f64,         // E_rest
    pub e_cross: f64,        // E_cross（= 库仑 + 交换 + XC 残差）
    pub e_cross_coul: f64,   // Tr(D_bond J_rest)
    pub e_cross_exch: f64,   // -a_x/2·Tr(D_bond K_rest)
    pub e_cross_xc: f64,     // E_xc[rho]-E_xc[rho_bond]-E_xc[rho_rest]
    pub e_elec: f64,         // E_elec
    pub e_nuc: f64,          // E_nuc（整个分子）
    pub e_total: f64,        // E_elec + E_nuc
}
// BondEnergy::bond() == e_bond + e_cross
```

== 6. 数值示例：C6H6（B3LYP/6-31G(d)）

单位 Hartree。log 中 `E(RB3LYP) = -232.247385084`。取苯环上相邻的 C–C 键（原子 0–1）与一根 C–H 键（原子 0–6；序号从 0 开始）：

#table(
  columns: 8,
  align: center,
  inset: 4pt,
  stroke: 0.5pt + luma(170),
  [*键*], [$E_("bond")$], [$E_("rest")$], [$E_("cross")$], [$E_("elec")$], [$E_("nuc")$], [$E_("total")$], [*键贡献*],
  [C0–C1], [$-168.171825$], [$-351.591949$], [$+84.431702$], [$-435.332072$], [$+203.081472$], [*$-232.250600$*], [*$-83.740123$*],
  [C0–H6], [$-93.973199$], [$-396.748853$], [$+55.389980$], [$-435.332072$], [$+203.081472$], [*$-232.250600$*], [*$-38.583219$*],
)

交叉项 $E_("cross")$ 的分项：

#table(
  columns: 5,
  align: center,
  inset: 4pt,
  stroke: 0.5pt + luma(170),
  [*键*], [*库仑* $op("Tr")(D_("bond") J_("rest"))$], [*交换* $-a_x/2 op("Tr")(D_("bond") K_("rest"))$], [*XC 残差*], [*合计*],
  [C0–C1], [*$+85.047721$*], [$-0.232037$], [$-0.383982$], [$+84.431702$],
  [C0–H6], [*$+55.764223$*], [$-0.138286$], [$-0.235956$], [$+55.389980$],
)

要点：

- 两种键都满足 $E_("bond") + E_("rest") + E_("cross") = E_("elec")$ *精确闭合*（测试断言误差 $< 10^(-6)$ Ha）。
- 两次分解给出完全相同的 $E_("elec") = -435.3321$、$E_("nuc") = +203.0815$、$E_("total") = -232.2506$，且与 π 能量分解（POCV / MOCV）的结果*完全一致*——它们分解的是同一个完整密度，这是实现正确性的强验证。
- $E_("total")$ 与 log 的 `SCF Done` 相差仅约 *3 mHa*（来自数值格点与积分/泛函微差）。
- 交叉项均以*键–其余库仑排斥*为主（C–C $+85.05$，C–H $+55.76$），交换与 XC 残差为小负，机制与 $E_(pi sigma) > 0$ 相同。
- 键贡献分别约 $-83.74$ Ha（C–C）与 $-38.58$ Ha（C–H），且各自与 $E_("rest")$ 之和都等于 $E_("elec")$。这是"该两原子密度块对电子能量的份额"，*不是键解离能*：它包含核心电子与全部核的未屏蔽吸引，随后又被 $E_("cross")$ 中的库仑排斥大量抵消。
- *C–H 的键贡献绝对值明显小于 C–C，不能据此认为 C–H 键更弱*：两者的密度块电子数与空间范围不同（H 只有一个价电子），比较的是"电子能量份额"而非键强度。真正的化学键能需要另外做能量差（如分子与碎片的 $E_("total")$ 之差）。

== 7. 约定与注意事项

- *适用范围*：仅 log/out 文件、闭壳层、且泛函可被 `XcFunc` 识别（HF 或 `m06l` 等暂不支持）。与 π 能量分解共用同一套 `prepare_log` 限制。
- *XC 不可分离*：$E_("xc")[rho] != E_("xc")[rho_("bond")] + E_("xc")[rho_("rest")]$，本实现把差值归入 $E_("cross")$。$E_("bond")$、$E_("cross")$ *各自的数值依赖此约定*；稳健的是三者之和与 $E_("elec")$。
- *量纲*：键贡献是"电子能量份额"（大负数），不是键解离能，也不是分子总能量；要与总能量同量级比较，需再按约定分摊核排斥能 $E_("nuc")$。
- *$D_("rest")$ 可能局部为负*：$D_("rest") = D - D_("bond")$ 一般不是半正定矩阵，在键密度占优的格点 $rho_("rest")$ 可能略小于 0；常规体系（如苯的 C–C 键）不会触发异常。若遇特殊体系需留意。
- *这是"一对原子 vs 其余"的二分解*：不同原子对的分解*彼此不可加*（每个结果都是同一 $E_("elec")$ 的一种二分法），不能把多个键的 `bond()` 相加当作总键能。
- *与 π 分解的关系*：两者模板相同（$D = D_X + D_"其余"$ 的精确加法分解），区别只在"分支"的选取——π 取轨道子空间投影，键取 $A union B$ 的 AO 块投影。
- *性能*：需构建一次 100×74 数值格点并做 3 次 $E_("xc")$ 积分（$rho, rho_("bond"), rho_("rest")$）与 2 次 `jk_mat`；C6H6 集成测试（debug，C–C 与 C–H 各一次）约 50 s。

== 8. 命令行与库调用

=== 8.1 命令行

```bash
# C1–C2 键的能量分解（原子序号从 1 开始）
rswfn bondprop energy 1 2 moles/C6H6.out
```

输出包含键/其余/交叉三项、交叉项的库仑/精确交换/XC 残差分解、完整电子能量、
核排斥能、总能量，以及键的可加贡献 $E_("bond") + E_("cross")$。

=== 8.2 库调用

```rust
use rswfn::base::Mole;
use rswfn::props::bondprop::energy::Calculator;

let mole = Mole::from_file("moles/C6H6.out");
let caler = Calculator::new(&mole);

// 只要键的可加贡献 E_bond + E_cross
let e_bond = caler.bond_energy(0, 1).unwrap();

// 完整分解
let parts = caler.bond_parts(0, 1).unwrap();
println!("E_bond = {}", parts.e_bond);
println!("E_cross = {}", parts.e_cross);
println!("E_elec = {}", parts.e_elec);
```

注意库调用使用 0 基原子序号，命令行使用 1 基。
