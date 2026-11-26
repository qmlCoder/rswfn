#set page(paper: "a4", margin: (x: 2cm, y: 2cm), numbering: "1")
#set text(font: ("Libertinus Serif", "Microsoft YaHei"), size: 10.5pt, lang: "zh")
#set par(justify: true, leading: 0.85em)
#show raw.where(block: true): set text(size: 8.5pt)
#show table: set text(size: 9pt)

= RESP 电荷的计算原理

#block(fill: luma(242), inset: 10pt, radius: 3pt, width: 100%)[
  *对应代码*：`src/props/atomprop/_resp.rs`（`Calculator` / `RespOptions` / `RespResult` / `mk_points` / `esp_values` / `normal_equations` / `fit_stage` / `identify_ch_groups` / `unit_sphere_points`）\
  *静电势*：`src/props/gridprop/density.rs`（`nuc_esp` / `ele_esp` / `tot_esp`）、`src/rsint/esp.rs`（`EspCalculator`）\
  *命令行*：`src/cli.rs`（`atomprop charge resp`、`print_resp_charges`）\
  *参考实现*：Multiwfn 3.8 `population.f90` 中的 `RESP`、`RESPiter`、`setMKpt`、`unitspherept`、`setESPfitvdwr`\
  *集成测试*：`tests/atomprop.rs`（`test10` 乙烯、`test11` 苯）
]

== 0. 结论速览

- RESP 是*静电势拟合电荷*：在分子周围生成一批拟合点，用原子电荷 $q_i$ 的库仑势
  $sum_i q_i slash r_i$ 去逼近这些点上的真实静电势 $V_k$。
- 纯最小二乘会让埋藏原子获得不物理的大电荷，因此引入*双曲限制*（hyperbolic
  restraint）把电荷拉向 0。限制强度用参数 $a$ 控制，$b$ 控制限制函数在
  $q approx 0$ 附近的平滑程度。
- 标准 RESP 采用*两阶段*策略：第一阶段用弱限制拟合所有原子；第二阶段固定
  sp³ 碳/亚甲基碳以外的原子，只对其余 CH 部分用强限制再拟合，并约束同一碳上的
  氢电荷相等。
- 总电荷约束 $sum_i q_i = Q_"mol"$ 始终保留。
- 本实现与 Multiwfn 3.8 的默认参数一致（MK 网格、$a_1=0.0005$、$a_2=0.001$、
  $b=0.1$），对 C6H6 和苯酚的电荷与 Multiwfn 相差约 $10^(-4)$ e 量级
  （最大 $1.1 times 10^(-4)$ e）。

== 1. 从静电势到线性最小二乘

=== 1.1 拟合方程

设有 $N_"pt"$ 个拟合点 $bold(r)_k$，其上的静电势为 $V_k$；分子有 $N_"atom"$ 个
原子，中心为 $bold(R)_i$，待求电荷为 $q_i$。记

$
  r_(k i) = abs(bold(r)_k - bold(R)_i),
  quad A_(k i) = 1 / r_(k i).
$

若电荷模型完全准确，应有

$
  V_k = sum_i A_(k i) q_i.
$

由于 $N_"pt"$ 远大于 $N_"atom"$，一般无精确解，只能最小化残差平方和：

$
  chi^2 (bold(q))
  = sum_k (V_k - sum_i A_(k i) q_i)^2.
$

=== 1.2 法方程

对 $q_j$ 求偏导并令其为零：

$
  (partial chi^2) / (partial q_j)
  = -2 sum_k A_(k j) (V_k - sum_i A_(k i) q_i)
  = 0.
$

整理得到线性方程组

$
  M bold(q) = bold(b),
$

其中

$
  M_(i j)
  = sum_k A_(k i) A_(k j)
  = sum_k 1 / (r_(k i) r_(k j)),
  quad
  b_i
  = sum_k A_(k i) V_k
  = sum_k V_k / r_(k i).
$

$M$ 是 $N_"atom" times N_"atom"$ 的对称矩阵。程序里先用一次遍历构造 $M$ 和
$bold(b)$，两个 RESP 阶段共用同一套法方程。

== 2. 双曲限制与迭代求解

=== 2.1 目标函数

直接最小化 $chi^2$ 时，$M$ 往往接近奇异，拟合电荷对拟合点的微小变化很敏感，
容易给出 $plus.minus 2$ e 量级的埋藏原子电荷。RESP 在目标函数中加入双曲惩罚：

$
  F(bold(q))
  = 1/2 sum_k (V_k - sum_i A_(k i) q_i)^2
    + a sum_(i in.not H) (sqrt(q_i^2 + b^2) - b).
$

- 惩罚项只施加在*非氢原子*上；氢原子通过第二阶段的等价约束控制。
- 当 $abs(q_i)$ 远大于 $b$ 时，$sqrt(q_i^2+b^2)-b approx abs(q_i)-b$，惩罚近似线性；
  当 $abs(q_i) < b$ 时，惩罚近似 $q_i^2 slash (2b)$，是光滑的二次型，避免
  在 $q_i=0$ 处出现不可导的 $abs(q_i)$。
- Multiwfn 的默认 $b=0.1$ e，$a_1=0.0005$（第一阶段），$a_2=0.001$（第二阶段）。

#block(fill: luma(245), inset: 8pt, radius: 2pt, width: 100%)[
  *注意因子 2*：若把惩罚写成 $(a slash 2) sum (sqrt(q^2+b^2)-b)$，则代码里的 `a`
  需要乘 2。本实现与 Multiwfn 保持一致，采用上面的 $a sum (...)$ 形式。
]

=== 2.2 非线性项的处理

对 $F$ 求偏导：

$
  (partial F) / (partial q_j)
  = -sum_k A_(k j) (V_k - sum_i A_(k i) q_i)
    + a q_j / sqrt(q_j^2 + b^2)
  = 0.
$

即

$
  sum_i M_(j i) q_i
  + a / sqrt(q_j^2 + b^2) q_j
  = b_j.
$

这是一个非线性方程组。程序采用不动点迭代：用上一轮的 $q_i^((k))$ 计算对角项，
把

$
  Lambda_i^((k)) = a / sqrt((q_i^((k)))^2 + b^2)
$

加到 $M$ 的对角元上，然后解线性方程组

$
  (M + op("diag")(Lambda^((k)))) bold(q)^((k+1)) = bold(b).
$

迭代直到 $max_i abs(q_i^((k+1)) - q_i^((k))) < epsilon$，默认
$epsilon = 10^(-6)$ e、最大 300 次。若 $a=0$，则只解一次线性方程。

=== 2.3 为什么这样能工作

双曲限制的作用可以理解为“软约束”：电荷越大，对角惩罚越小，但始终存在；
它不会像硬约束那样把电荷钉死，而是在“拟合 ESP”和“保持电荷合理”之间取平衡。
第一阶段的弱限制给出整体合理的电荷分布；第二阶段的强限制只作用于 CH 部分，
用来稳定脂肪碳/氢的电荷，同时固定其它原子以保持与第一阶段一致。

== 3. 等式约束

=== 3.1 总电荷约束

分子总电荷必须守恒：

$
  sum_i q_i = Q_"mol",
  quad
  Q_"mol" = sum_A Z_A^"eff" - N_"elec".
$

$Z_A^"eff"$ 是有效核电荷（赝势体系也自洽），$N_"elec"$ 是总电子数。
程序用 `mole.charge()` 得到 $Q_"mol"$。

=== 3.2 额外电荷约束

第二阶段需要把非 CH 原子固定在第一阶段电荷上：

$
  q_i = q_i^"stage1",
  quad i in.not "active".
$

更一般地，可以写成若干线性约束

$
  sum_(i in S_c) q_i = Q_c,
  quad c = 1, dots, N_"cons".
$

=== 3.3 增广线性系统

把等式约束用拉格朗日乘子并入法方程，得到对称不定系统

$
  mat(
    delim: "(",
    M + Lambda, bold(1), C^top;
    bold(1)^top, 0, 0;
    C, 0, 0
  )
  mat(delim: "(", bold(q); lambda; bold(mu))
  =
  mat(delim: "(", bold(b); Q_"mol"; bold(Q)_"cons")
$

其中 $bold(1)$ 是全 1 列向量，$C$ 的每一行对应一个电荷约束（该行在约束涉及的
原子位置为 1），$bold(Q)_"cons"$ 是各约束的目标值。解出的 $bold(q)$ 就是该阶段
的电荷，乘子 $lambda$、$bold(mu)$ 不参与后续计算。

程序实际构造的是“收缩后”的矩阵，见下一节。

== 4. 等价约束与矩阵收缩

如果一组原子 $g$ 必须具有相同电荷，则令 $q_i = tilde(q)_g$（$i in g$），
把法方程的行和列按组求和：

$
  tilde(M)_(g h) = sum_(i in g) sum_(j in h) M_(i j),
  quad
  tilde(b)_g = sum_(i in g) b_i.
$

约束行/列也要相应收缩：

- 总电荷约束中，组 $g$ 的系数是组大小 $abs(g)$；
- 电荷约束 $c$ 中，组 $g$ 的系数是“组 $g$ 中属于 $S_c$ 的原子个数”。

这与 Multiwfn 的做法等价：Multiwfn 先构造未收缩的增广矩阵，再把属于同一等价组
的行相加、列相加；本实现直接按组累加，避免显式构造未收缩矩阵。

== 5. MK 拟合点

=== 5.1 原子半径

MK 网格使用一套专门用于 ESP 拟合的 vdW 半径（单位 Å）：

#table(
  columns: 9,
  align: center,
  inset: 4pt,
  stroke: 0.5pt + luma(170),
  [*H*], [*C*], [*N*], [*O*], [*F*], [*P*], [*S*], [*Cl*], [*Z > 17*],
  [1.20], [1.50], [1.50], [1.40], [1.35], [1.80], [1.75], [1.70], [UFF $slash$ 1.2],
)

表中只列出常见元素；其余 H–Cl 元素仍使用 Multiwfn 的 MK 专用表，
$Z > 17$ 的元素使用 UFF 半径除以 1.2。
程序内部坐标是 Bohr，因此这些 Å 值会乘 `BOHR`（Å → Bohr）后使用。

=== 5.2 层与点数

每个原子生成 4 层球面点，半径倍数为

$
  s in {1.4, 1.6, 1.8, 2.0}.
$

第 $i$ 个原子、第 $s$ 层的球半径为

$
  R_(i s) = s dot R_i^"vdW".
$

球面上的期望点数为

$
  n_(i s)
  = "nint"(4 pi R_(i s)^2 rho_"MK"),
  quad
  rho_"MK" = 6 "点/Å"^2.
$

在 Bohr 单位下 $rho_"MK" = 6 slash "BOHR"^2 approx 1.68017 "Bohr"^(-2)$。
`nint` 表示四舍五入。

=== 5.3 球面撒点：`unitspherept`

Multiwfn 不用 Lebedev 网格，而是用简单的经纬度采样。给定期望点数 $n$：

```rust
n_equ  = int(sqrt(n * pi));       // 每个纬圈的最大点数
n_vert = n_equ / 2;               // 纬圈数（整数除法）

for i_v in 0..=n_vert {
    theta   = pi * i_v / n_vert;
    r_xy    = sin(theta);
    z       = cos(theta);
    n_xy    = max(1, int(n_equ * r_xy));

    for i_h in 1..=n_xy {
        phi = 2 * pi * i_h / n_xy;
        push(cos(phi) * r_xy, sin(phi) * r_xy, z);
        if points.len() == n { return points; }
    }
}
```

注意实际返回的点数*可能少于* $n$：循环在超过 $n$ 时提前结束，且
$n_"xy" = max(1, "int"(n_"equ" r_"xy"))$ 会带来取整损失。程序逐行复刻了这套
逻辑，因此与 Multiwfn 的球面点集在多数层上逐点一致。

=== 5.4 剔除内部点

对某个原子 $i$ 生成的候选点 $bold(p)$，若它落在任何*其它*原子 $j$ 的最内层
vdW 球内，则丢弃：

$
  abs(bold(p) - bold(R)_j) < 1.4 R_j^"vdW"
  quad ==> quad "discard".
$

这一步保证拟合点位于分子“外表面”附近，而不是埋在其他原子内部。程序只考虑
真实原子（本实现没有额外拟合中心）。

=== 5.5 伪代码

```rust
for i in 0..natom {
    for s in [1.4, 1.6, 1.8, 2.0] {
        R = s * vdw[i];
        n = round(4 * pi * (R / BOHR)^2 * 6.0);
        for dir in unitspherept(n) {
            p = center[i] + R * dir;
            if all(dist(p, center[j]) >= 1.4 * vdw[j] for j != i) {
                points.push(p);
            }
        }
    }
}
```

对 C6H6，本实现生成 9066 个点，Multiwfn 生成 9086 个点。差异来自不同平台的
`sin/cos` 末位差异与单位换算常数（Multiwfn 用 `b2a = 0.5291772`，本实现用更
精确的 `BOHR`），属于 0.1% 量级的采样差异。

== 6. 静电势的计算

拟合点上的 $V_k$ 用解析高斯积分计算，而不是实空间格点积分：

$
  V(bold(r))
  = sum_A Z_A^"eff" / abs(bold(r) - bold(R)_A)
    - integral (rho(bold(r)')) / abs(bold(r) - bold(r)') dif bold(r)'.
$

- 核项：遍历原子累加 $Z_A^"eff" slash abs(bold(r)-bold(R)_A)$；落在核上时返回
  `1e3` 以避免无穷大（与 Multiwfn 一致）。
- 电子项：由密度矩阵的高斯型基函数解析积分得到，包含电子负号，结果为负值。
- 两项相加即 `tot_esp`，单位 Hartree。坐标全部是 Bohr，$1 slash r$ 直接给出
  Hartree，不需要额外常数。

详见 `docs/静电势.md` 与 `src/rsint/esp.rs`。

== 7. 标准两阶段 RESP

=== 7.1 第一阶段：弱限制

- 拟合中心：所有原子。
- 等价约束：无（每个原子独立）。
- 电荷约束：只有总电荷。
- $a_1 = 0.0005$，$b = 0.1$。
- 惩罚只加在非氢原子上。
- 迭代到 $10^(-6)$ e 或 300 次。

得到 `stage1` 电荷。

=== 7.2 第二阶段：强限制 + CH 约束

首先识别“活跃”碳原子：

- 碳原子的配位数等于 4（sp³ 碳，包括 CH₃、sp³ CH₂、sp³ CH 等）；
- 或者碳原子恰好连接 2 个 H 和 1 个 C（=CH₂ 或 -CH₂-）。

活跃集合 = 这些碳原子 + 直接连在它们上面的氢原子。

第二阶段：

- 非活跃原子固定为第一阶段电荷：$q_i = q_i^"stage1"$；
- 同一个活跃碳上的氢原子共享电荷（等价约束）；
- 活跃碳本身仍作为独立拟合中心（带有 $a_2$ 的双曲限制）；
- $a_2 = 0.001$，$b = 0.1$；
- 总电荷约束保留。

若没有活跃碳（例如苯、苯酚），第二阶段直接跳过，最终电荷等于第一阶段结果。

=== 7.3 参数汇总

#table(
  columns: 4,
  align: (left, center, center, left),
  inset: 4pt,
  stroke: 0.5pt + luma(170),
  [*参数*], [*默认值*], [*代码字段*], [*作用*],
  [层半径倍数], [$1.4, 1.6, 1.8, 2.0$], [`layer_scales`], [MK 四层球面],
  [点密度], [6 点/Å²], [`pt_density`], [每 Å² 球面点数],
  [第一阶段限制], [$0.0005$], [`hyper_a1`], [弱双曲限制],
  [第二阶段限制], [$0.001$], [`hyper_a2`], [强双曲限制],
  [tightness], [$0.1$], [`hyper_b`], [双曲函数参数 $b$],
  [最大迭代], [$300$], [`max_iter`], [电荷收敛迭代上限],
  [收敛阈值], [$10^(-6)$], [`conv`], [最大电荷变化],
  [阶段开关], [true], [`two_stage`], [是否执行第二阶段],
)

=== 7.4 流程

```rust
// 1. 生成 MK 点并计算 ESP
let points = caler.mk_points(&opts);
let esp    = caler.esp_values(&points, EspType::Total);

// 2. 构造法方程（两个阶段共用）
let (M, b) = normal_equations(&centers, &points, &esp);

// 3. 第一阶段：所有原子独立
let stage1 = fit_stage(M, b, total_charge, singletons, [], a1, b);

// 4. 第二阶段：固定非 CH 原子，约束 CH 氢等价
let (ch_list, groups) = identify_ch_groups(&mole);
if !ch_list.is_empty() {
    let constraints = atoms_not_in(ch_list)
        .map(|i| (vec![i], stage1[i]));
    let charges = fit_stage(M, b, total_charge, groups, constraints, a2, b);
} else {
    let charges = stage1;
}

// 5. RMSE / RRMSE
let (rmse, rrmse) = fit_errors(&points, &esp, &centers, &charges);
```

== 8. 数值验证

用 WSL 中的 Multiwfn 3.8 `Multiwfn_noGUI` 对同一批 `.fch` 文件运行
“7 → 18 → 1” 的 RESP 流程，与本实现比较。

#table(
  columns: 4,
  align: (left, center, center, center),
  inset: 5pt,
  stroke: 0.5pt + luma(170),
  [*分子*], [*MK 点数（本实现 / Multiwfn）*], [*最大电荷偏差*], [*RMSE（本实现 / Multiwfn）*],
  [C6H6], [9066 / 9086], [$5.5 times 10^(-5)$ e], [$1.2513 times 10^(-3)$ / $0.001251$],
  [苯酚], [9574 / 9580], [$1.1 times 10^(-4)$ e], [$1.6839 times 10^(-3)$ / $0.001684$],
)

C6H6 的代表性原子（单位 e）：

#table(
  columns: 4,
  align: (left, center, center, center),
  inset: 4pt,
  stroke: 0.5pt + luma(170),
  [*原子*], [*本实现*], [*Multiwfn*], [*差值*],
  [C1], [$-0.107960$], [$-0.107916$], [$-4.4 times 10^(-5)$],
  [C2], [$-0.111595$], [$-0.111639$], [$+4.5 times 10^(-5)$],
  [C5], [$-0.113418$], [$-0.113363$], [$-5.5 times 10^(-5)$],
  [H7], [$+0.110097$], [$+0.110070$], [$+2.7 times 10^(-5)$],
  [H8], [$+0.110911$], [$+0.110926$], [$-1.5 times 10^(-5)$],
  [H11], [$+0.111965$], [$+0.111922$], [$+4.3 times 10^(-5)$],
)

乙烯的两阶段测试（`test10`）：

- 阶段 1：C $= -0.297456$，H $= +0.148728$；
- 阶段 2：C $= -0.297035$，H $= +0.148517$；
- 总电荷为 0，同一 CH₂ 上的两个 H 严格相等，说明等价约束生效。

#block(fill: luma(245), inset: 8pt, radius: 2pt, width: 100%)[
  *说明*：Multiwfn 无法读取当前 `FchWriter` 写出的乙烷 `.fch`（基组块格式问题），
  因此乙烷类 sp³ 体系的第二阶段暂未做数值对照；第二阶段的约束/固定逻辑是按
  Multiwfn 源码逐行实现的。
]

== 9. 实现映射与用法

=== 9.1 代码结构

#table(
  columns: 2,
  align: (left, left),
  inset: 4pt,
  stroke: 0.5pt + luma(170),
  table.header([*文件 / 函数*], [*作用*]),
  [`src/props/atomprop/_resp.rs`], [RESP 主实现],
  [`Calculator::resp` / `resp_with`], [入口：默认 / 自定义参数],
  [`Calculator::mk_points`], [MK 拟合点],
  [`Calculator::esp_values`], [并行计算拟合点 ESP],
  [`normal_equations`], [构造 $M$ 与 $bold(b)$],
  [`fit_stage`], [带双曲限制与等式约束的迭代],
  [`identify_ch_groups`], [识别第二阶段活跃 CH，生成等价组],
  [`unit_sphere_points`], [复刻 `unitspherept`],
  [`resp_vdw_radius` / `uff_radius`], [MK / UFF 半径],
  [`src/props/gridprop/density.rs`], [解析 ESP（`tot_esp`）],
  [`src/cli.rs`], [命令行入口与结果打印],
  [`tests/atomprop.rs`], [C2H4 / C6H6 集成测试],
)

=== 9.2 命令行

```bash
# 默认标准两阶段
rswfn atomprop charge resp moles/C6H6.fch

# 自定义 MK 点密度
rswfn atomprop charge resp moles/phenol.fch --pt-density 6.0

# 只做第一阶段（弱限制）
rswfn atomprop charge resp moles/C2H4.out --one-stage
```

输出包含：阶段 1 电荷、最终 RESP 电荷、总电荷、ESP 拟合 RMSE/RRMSE、
两阶段迭代次数和收敛状态。

=== 9.3 库调用

```rust
use rswfn::props::atomprop::resp::{Calculator, EspType, RespOptions};
use rswfn::base::Mole;

let mole = Mole::from_file("moles/C2H4.out");
let caler = Calculator::new(&mole);

// 默认参数
let res = caler.resp().unwrap();
println!("charges = {:?}", res.charges);
println!("RMSE = {:.3e}", res.rmse);

// 自定义参数
let opts = RespOptions {
    pt_density: 6.0,
    hyper_a1: 0.0005,
    hyper_a2: 0.001,
    hyper_b: 0.1,
    esp_type: EspType::Total,
    ..RespOptions::default()
};
let res = caler.resp_with(&opts).unwrap();
```

== 10. 单位、边界与限制

- *单位*：内部坐标 Bohr，静电势 Hartree，电荷 e。MK 半径表用 Å 存储，使用时
  乘 `BOHR` 转到 Bohr。
- *总电荷*：用 `mole.charge()`（有效核电荷和减总电子数），因此带电分子和
  ECP 体系也能正确约束。
- *仅 MK 网格*：尚未实现 CHELPG、外部拟合点、额外拟合中心和多构象加权。
- *仅单构象*：当前不读取构象列表。
- *约束*：`fit_stage` 内部支持任意等价组和线性电荷约束，但 CLI 只暴露了
  标准两阶段流程；自定义约束尚未做成命令行选项。
- *导出*：CLI 只打印结果，尚未写 Amber 风格 `.chg` 文件。
- *数值条件*：法方程 $M$ 条件数可能较大；当前用 LU 求逆。对超大体系，未来可
  改用 QR/SVD 或正则化最小二乘。
- *成键判定*：第二阶段用 `atoms.nebs()`，其判据是共价半径和乘 1.1；Multiwfn
  默认是 1.15。对键长异常的体系，CH 识别可能有差异。

== 11. 参考文献

+ Bayly, C. I.; Cieplak, P.; Cornell, W. D.; Kollman, P. A.
  A well-behaved electrostatic potential based method using charge restraints for
  deriving atomic charges: the RESP model.
  _J. Phys. Chem._ *1993*, 97, 10269–10280.
+ Cieplak, P.; Cornell, W. D.; Bayly, C.; Kollman, P. A.
  Application of the multimolecule and multiconformational RESP methodology to
  biopolymers.
  _J. Comput. Chem._ *1995*, 16, 1357–1377.
+ Singh, U. C.; Kollman, P. A.
  An approach to computing electrostatic charges for molecules.
  _J. Comput. Chem._ *1984*, 5, 129–145.
+ Rappé, A. K.; Casewit, C. J.; Colwell, K. S.; Goddard, W. A.; Skiff, W. M.
  UFF, a full periodic table force field for molecular mechanics and molecular
  dynamics simulations.
  _J. Am. Chem. Soc._ *1992*, 114, 10024–10039.
+ Bondi, A. van der Waals volumes and radii.
  _J. Phys. Chem._ *1964*, 68, 441–451.
+ Lu, T.; Chen, F. Multiwfn: a multifunctional wavefunction analyzer.
  _J. Comput. Chem._ *2012*, 33, 580–592.
+ Lu, T. A comprehensive electron wavefunction analysis toolbox for chemists,
  Multiwfn.
  _J. Chem. Phys._ *2024*, 161, 082503.
+ Tian Lu, _Multiwfn manual_（RESP 章节），version 3.8.
