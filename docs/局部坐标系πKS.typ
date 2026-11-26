#set page(paper: "a4", margin: (x: 2cm, y: 2cm), numbering: "1")
#set text(font: ("Libertinus Serif", "Microsoft YaHei"), size: 10.5pt, lang: "zh")
#set par(justify: true, leading: 0.85em)
#show raw.where(block: true): set text(size: 8.5pt)
#show table: set text(size: 9pt)

= 局部坐标系下的 π-only Kohn-Sham（LKS）

#block(fill: luma(242), inset: 10pt, radius: 3pt, width: 100%)[
  *对应代码*：`src/qmcal/lks.rs`（`LKS` / `LKSResult`）\
  *服务与命令行*：`src/service/qmcal.rs`（`lks`）、`src/cli.rs`（`QmcalCmds::Lks`）\
  *局部坐标系*：`src/props/atomprop/direction.rs`（`LCS` / `RLCS`）\
  *变换矩阵*：`src/props/orbtprop/mocv.rs`（`get_tmat_car`）\
  *局部积分*：`src/rsint/mod.rs`（`mat1e_lcs` / `mat2e_lcs`）\
  *导出*：`src/writer/fch.rs`（`FchWriter`）、`src/writer/mfn.rs`（`MfnWriter`）\
  *集成测试*：`tests/qmcalc.rs`（`test_lks_c2h4`）
]

== 0. 结论速览

- 在*原子局部坐标系*中计算单/双电子积分，只保留 π 对称的原子轨道，得到远小于全基组的 *π 约化基*；
- 在该约化基上做闭壳层 KS 自洽场，*Fock 矩阵与系数矩阵都只有 $n_pi times n_pi$*；
- DIIS 外推、本征求解、格点交换关联积分同样在 $n_pi$ 维进行；
- 这是 *π 子系统模型*，能量不含 σ 电子与核内层，*不是* 分子总能量（与 `docs/pi电子能量.md` 是同一个坑）；
- 双电子积分的四级变换仍复用 `mat2e_lcs`（全维变换后再截断），因此*积分耗时并未减少*；变小的是 SCF 内循环。

== 1. 局部坐标系（LCS）

对每个原子构造正交右手系（`Stm`，`mat` 的每一列是一个基向量）：

- $bold(e)_x$：分子平面外法向 $bold(n)$；
- $bold(e)_z$：沿某条键的方向；
- $bold(e)_y = bold(e)_x times bold(e)_z$。

于是局部 $x$ 轴指向平面外，*π 轨道就是局部 x 指数为奇数的分量*。无法确定法向的原子（如线性分子端点）退回随机坐标系（`RLCS`）。

== 2. π 约化基

壳层 $s$ 的局部笛卡尔基 $chi^l$ 与全局基 $chi^g$ 的关系为

$ chi^l = chi^g T^T, quad v^l = T v^g, quad S^l = T S^g T^T, $

其中 $T = "get_tmat_car"(d, L)$（$d$ 为壳层笛卡尔维数，$L$ 为 LCS 旋转矩阵）。局部积分由 `mat1e_lcs` / `mat2e_lcs` 给出后，只保留*局部 x 指数为奇数*的原子轨道索引集 $cal(P)$：

$ S_pi = S^l [cal(P), cal(P)], quad H_pi = H^l [cal(P), cal(P)], quad (p q | r s)_pi = (p q | r s)^l [cal(P)^4]. $

对 `p` 壳层是 $p_x$；对笛卡尔 `d` 壳层是 $d_(x y), d_(x z)$；对 `f` 是 $f_(x^3), f_(x y^2), f_(x z^2), f_(x y z)$，与 `mocv_syms(PiProperty)` 选出的 π 子空间一致。分裂价基组每个原子可能有多条 `p` 壳层，全部计入。

== 3. 约化自洽场

Roothaan–Hall 方程在 π 约化基上求解：

$ F_pi C_pi = S_pi C_pi epsilon_pi, quad D_pi = 2 C_(pi,"occ") C_(pi,"occ")^T, $

$ F_pi = H_pi + J[D_pi] - a_x/2 K[D_pi] + V_("xc")[rho_pi]. $

其中 $J, K$ 由 π 双电子积分与 $D_pi$ 构造（`rsint::jk_mat`）。交换关联势 $V_("xc")$ 在数值格点上积分：格点处 π 原子轨道的值/梯度由全局笛卡尔 AO 值经同一 $T$ 旋转后取 π 列。

迭代初猜取 $H_pi$ 的核心哈密顿量，用 Pulay DIIS 加速；收敛判据为相邻两次能量变化小于 `tol`（默认 $10^(-8)$ Ha）。

== 4. π 电子数

- `LKS::npi_elec = Some(n)` 时直接使用 $n$（必须为偶数）；
- 为 `None` 时按“拥有 `p` 壳层的原子数”估计：每个提供 `p` 壳层的原子贡献一个 π 电子（乙烯 2，苯 6）；
- 对吡咯 N（孤对贡献 2）、带电体系等该估计可能偏差，应显式指定。

占据轨道数 $n_("occ") = n_pi / 2$，按能量从低到高填充。

== 5. 能量与物理意义

电子能量沿用完整 KS 泛函表达式，但代入的是 π 约化密度：

$ E_pi = op("Tr")(D_pi H_pi) + 1/2 op("Tr")(D_pi J[D_pi]) - a_x/4 op("Tr")(D_pi K[D_pi]) + E_("xc")[rho_pi]. $

#block(fill: luma(242), inset: 10pt, radius: 3pt, width: 100%)[
  *警告*：$E_pi$ 是“把若干 π 电子放进裸核场”的模型量，缺少 σ 电子的屏蔽，数值可达数十 Hartree 的大负数，*不是* 分子总能量。核排斥能 $E_("nuc")$ 是完整分子的值，不能与之简单相加当作总能量。详见 `docs/pi电子能量.md`。
]

它的用途是以很小的矩阵代价获得 *π 分子轨道与轨道能*（以及 π 密度），在定性分析和调试 π 子空间时有用。

== 5.1 打印与导出

- *局部 π 系数打印*（`--print-mo`）：输出行 = π 原子轨道、列 = π 分子轨道的系数表，
  行标签只含 π 成分，如 `C1:PX`、`C3:XY`（`PX` 表示局部 x，即面外 π 分量）。
- *展开回正常（全局）系数*（`LKSResult::expand_coefs` / `to_mole`）：
  $C_("glob") = P dot C_pi$，其中 $P$ 把每个局部 π 原子轨道展开到全局笛卡尔基。
- *导出*（`-o`）：`.fch/.fchk` 或 `.mwfn`。文件描述的是 π 子体系：
  基函数数 $N_("basis")$ 为完整基组，轨道数 $N_("indbasis") = n_pi < N_("basis")$，
  电子数写作 π 电子数。因此它与核电荷总数不自洽，仅供可视化 π 轨道用。

== 6. 数值示例

#table(
  columns: (auto, auto, auto, auto, auto),
  align: center,
  inset: 5pt,
  stroke: 0.5pt + luma(170),
  [*体系 / 基组*], [*全基*], [*π 约化基*], [*π 电子*], [*SCF 迭代*],
  [乙烯 B3LYP/6-31G(d)], [38], [8], [2], [4],
  [苯 B3LYP/STO-3G], [36], [6], [6], [5],
)

苯/STO-3G 的 π 轨道能级（Hartree）：$-11.444$（a2u）、$-11.403$（e1g，二重简并）、$-11.319$（e2u，二重简并）、$-11.302$（b2g）。其能级次序与简并与苯的 π 分子轨道理论一致——这是局部坐标变换取到真正“平面外 π 分量”的有力验证。

== 7. 约定与注意事项

- *适用*：平面/共轭体系的闭壳层（$alpha = beta$）。非平面或 sp3 体系里“局部 x”不再对应 π，结果无意义。
- *积分耗时*：`mat2e_lcs` 仍是全维四级变换（甚至比 `mat2e` 更贵），π 截断发生在其后；若要把积分本身也做小，需要另写“只变换到 π 指标”的受限四级变换。
- *不导出 fch*：π 约化基不对应完整基组，导出（`.fch`/`.mwfn`）时 $N_("indbasis") < N_("basis")$、电子数为 π 电子数，与核电荷不自洽，仅供可视化；亦可用 `LKSResult::expand_coefs` 把系数展开回全局笛卡尔基（$C_("glob") = P dot C_pi$）自行处理。
- *电子数*：见第 4 节，杂原子体系请显式指定。

== 8. 命令行

```bash
rswfn qmcal lks -i mole.fch -b 6-31G(d) -f b3lyp
rswfn qmcal lks -i mole.fch --npi-elec 6 --print-mo -o pi.mwfn
```

输出 π 约化基大小、π 电子数、π 电子能量、交换相关能、核排斥能，以及 π 分子轨道能量与各阶段耗时；
`--print-mo` 打印局部 π 系数表，`-o` 导出 `.fch/.fchk` 或 `.mwfn`（系数已展开回全局基）。
