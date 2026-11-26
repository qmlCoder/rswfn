#set page(width: 21cm, height: 29.7cm, margin: (x: 2.2cm, y: 2.0cm))
#set text(font: ("New Computer Modern", "Microsoft YaHei"), size: 10.5pt, lang: "zh")
#set par(justify: true, leading: 0.72em)
#set heading(numbering: "1.1")

#align(center)[
  #text(size: 20pt, weight: "bold")[终端里的分子结构预览]
  #v(4pt)
  #text(size: 10pt, fill: rgb("#666666"))[rswfn 的 #raw("rswfn view") 子命令与 TUI #raw("F4") 预览窗口]
]

#v(6pt)

终端里没有像素，只有"字符格"。但一个字符格可以携带一个字符、一个前景色和一个
背景色，而半块字符（U+2580）恰好只画格子的上半部分——于是一个格子就能显示
#strong[两种独立的颜色]，等效于两个纵向像素。把这个技巧、正交投影和一点光照凑在
一起，纯文本终端也能显示有立体感的球棍模型。

#figure(
  image("preview/tui.png", width: 100%),
  caption: [TUI 里的预览窗口（`F4` 打开）。画面来自真实的终端缓冲，112×38 字符格],
)

#heading("用法")

命令行预览（非交互，可直接管道、截图或写进脚本）：

```bash
rswfn view C6H6.fch                        # 默认按当前终端大小铺满
rswfn view phenol.fch --size 100x44         # 指定画布尺寸（列×行）
rswfn view C6H6.fch --style space           # ball / space / wire
rswfn view C6H6.fch --rotate 30,-15 --zoom 1.2
rswfn view C6H6.fch --color always > shot.ans
rswfn view C6H6.fch --json                  # 原子坐标（Å）与成键，便于脚本处理
```

#table(
  columns: (auto, 1fr),
  inset: 6pt,
  stroke: 0.5pt + rgb("#cccccc"),
  align: (left, left),
  table.header([参数], [说明]),
  [`--size WxH`], [画布尺寸；省略时取终端大小（留一行给光标）],
  [`--style`], [`ball` 球棍、`space` 空间填充、`wire` 线框],
  [`--rotate yaw,pitch`], [绕竖直轴 / 水平轴的初始旋转角，单位度],
  [`--zoom`], [在"自动铺满"基础上的缩放倍数],
  [`--no-labels` / `--no-depth`], [关闭元素符号 / 关闭远处原子的景深雾化],
  [`--color`], [`auto` 只在终端上色、`always`、`never`（退化为灰度字符）],
)

TUI 预览窗口的按键：

#table(
  columns: (auto, 1fr),
  inset: 5pt,
  stroke: 0.5pt + rgb("#cccccc"),
  align: (left, left),
  table.header([按键], [作用]),
  [方向键 / `hjkl`], [旋转视角],
  [`+` / `-`], [缩放],
  [`s`], [球棍 → 空间填充 → 线框 循环],
  [`t` / `d`], [元素符号开关 / 景深开关],
  [`1` / `2` / `3`], [正视 / 侧视 / 俯视],
  [`r`], [复位视角],
  [`Esc` / `F4`], [关闭预览],
)

#heading("渲染管线")

#heading(level: 2)[主轴对齐]

分子文件里的坐标朝向是任意的，直接投影经常看到"侧着"的苯环。预览在打开时对坐标
协方差矩阵做一次特征分解（3×3 对称矩阵，用 Jacobi 旋转，不引入额外的线代依赖），
把摊得最开的两个主轴摆到屏幕的 x / y，最窄的方向留给深度。平面分子于是正对镜头，
非平面分子也总能拿到一个信息量最大的视角；再叠加一个很小的初始倾斜（约 20°/12°）
保留立体感。特征向量的正负号用坐标的三次矩定下来，保证同一分子每次打开朝向一致。

#heading(level: 2)[投影与自适应缩放]

坐标先减去几何中心，再按主轴旋转、按 `yaw` / `pitch` 旋转，最后做正交投影
（不做透视：终端分辨率本来就低，透视只会让近处的球畸形）。缩放系数取

$ k = 0.44 times min(W, 2H) \/ E_max, $

其中 $W$、$H$ 是字符格的列数与行数，$E_max$ 是投影后每个原子"坐标 + 球半径"的
最大外廓。行数要乘 2，因为每个字符格纵向有两个像素；而 $E_max$ 必须把球半径算
进去，否则旋转到某个角度时边上的原子会被画布边缘切掉（这一点有单元测试守着）。

#heading(level: 2)[像素化、着色与遮挡]

每个像素记录一个 RGB 值，初始为深色背景：

- 球体：按像素到球心的距离求球面法线，做 Lambert 漫反射，再加一点高光；
- 化学键：按圆柱侧面法线着色，沿键方向在两端原子的 CPK 颜色之间插值，
  并把圆柱裁剪到球面为止，避免"糖葫芦"；
- 景深：把深度归一化后向背景色混合，远处原子自然变暗；
- 遮挡：画家算法，按深度从远到近绘制，键的深度取两端平均值。

颜色来自 `src/datas/elems.rs` 里的 CPK 十六进制值，与常见的分子可视化程序一致
（C 灰、H 白、O 红、N 蓝、S 黄……）。

#heading(level: 2)[自适应输出]

同一份 `Frame` 可以导出三种形态：ratatui 的 `Line`（TUI）、ANSI 真彩色（CLI）、
按亮度映射的 ASCII 灰度（无色终端或重定向到文件）。三者共用同一份像素数据，
所以截图和实际运行看到的画面完全一致。

#figure(
  image("preview/styles.png", width: 100%),
  caption: [球棍、空间填充、线框三种模型，同一分子、同一视角],
)

#heading("分辨率与适用边界")

一个 $W times H$ 的终端，用半块字符相当于 $W times 2H$ 个彩色像素；用盲文
（Braille，2×4 点）能到 $2W times 4H$ 个单色点。按常见的 100×40 字符终端估算：

#table(
  columns: (auto, auto, 1fr),
  inset: 6pt,
  stroke: 0.5pt + rgb("#cccccc"),
  align: (left, left, left),
  table.header([体系规模], [可读性], [说明]),
  [10–30 原子], [很好], [环、取代基、顺反与构象差异都能一眼分辨],
  [30–100 原子], [尚可], [需要旋转、缩放看局部，标签会开始拥挤],
  [100 原子以上], [有限], [建议改用盲文骨架、局部放大，或直接用图形程序],
)

渲染成本可以忽略：100×40 的画布是 4000 个像素、几十个球和键，单帧在毫秒量级，
只在预览窗口打开时重绘。

#heading("和其他工具的分工")

终端预览的定位是"扫一眼"：SSH 进服务器、批量检查一叠输出文件、写脚本时确认
几何结构是否合理。精细的立体感、等值面、动画轨迹仍然交给 wmview、VMD、PyMOL
这类图形程序；两者数据来源相同，互不冲突。

#figure(
  image("preview/complex.png", width: 78%),
  caption: [Ni(NH3)2(OH)2 的球棍模型：CPK 配色 + 景深，100×46 字符格],
)

#heading("代码位置")

- `src/viewer.rs`：投影、光栅化、着色与三种输出格式，附带单元测试；
- `src/cli.rs`：`rswfn view` 的参数与出口；
- `src/tui.rs`：`F4` 弹层、按键处理与绘制；
- `src/base/atoms.rs`：`Atoms::bonds()` 提供构造时自动判定的成键列表。
