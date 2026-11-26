# 变更记录

## 2026-11-26 发布 v0.1.3（修复 Linux/musl 构建）

- v0.1.2 已上传 crates.io，但 GitHub Release 的构建在 `x86_64-unknown-linux-musl`
  失败：`rsbse 0.12.1` 硬依赖 `native-tls`，Linux 上落到 `openssl-sys`，需要系统 libssl
  开发包（`musl link check` 也在同一个地方挂）。
- 修复：`Cargo.toml` 增加 Linux 专属依赖
  `openssl-sys = { version = "0.9", features = ["vendored"] }`，让 OpenSSL 从源码编译，
  `cargo install` 和 musl 发布构建都不再依赖目标机上的 OpenSSL。`Cargo.lock` 新增
  `openssl-src 300.6.1+3.6.3`。
- 版本号 0.1.2 → 0.1.3：crates.io 上的 0.1.2 不可覆盖，且缺少上述修复。
- `cargo install rswfn` 取最高版本，所以 0.1.3 自动取代 0.1.2，无需 yank。
- README 里的版本示例同步更新为 0.1.3。

## 2026-11-26 发布 v0.1.2

- 版本号 `0.1.1` → `0.1.2`（crates.io 上 0.1.1 已于 2026-11-21 发布，不能重复发布）。
- crates.io 包的 `exclude` 新增 `/*.pdf`，避免仓库根目录的临时 PDF 被打进 crate。
- README 里 crates.io 与安装脚本示例的版本号（0.1.0 / 0.1.1）统一更新为 0.1.2。
- 发布前校验：`cargo package --list` 确认不含 `docs/`、`moles/`、PDF；
  `cargo test --offline --bin rswfn` 74 项全部通过。
- GitHub Release 由 `.github/workflows/release.yml` 在推送 `v0.1.2` tag 时自动构建
  四个平台产物（Windows x64、Linux x64 gnu/musl、Linux aarch64）并上传。

## 2026-11-26 预览旋转缩放稳定性排查（release 二进制滞后）

- 用户反馈“横向旋转时分子仍忽大忽小”。排查发现代码修复本身正确，问题在于 `which rswfn`
  指向的 `target/release/rswfn.exe`（构建于 11:55）早于缩放修复，而含修复的 debug 二进制
  构建于 12:03，两者表现不同：
  - debug（含修复）：yaw 0/30/60/89° 时渲染高度恒为 24，宽度 57/53/44/32 自然收窄。
  - 旧 release：高度 24/28/28/28，即旧的投影二维包围盒自适应，高度随 yaw 变化。
- 处理：`cargo clean -p rswfn --release` 后重建 release（耗时约 15m35s，仅剩既有的
  lib/bin 同名 `rswfn.pdb` 冲突警告）。新 release（12:23）实测高度恒为 24，确认修复生效。
- 结论：无需再改代码；重启 TUI 使用新二进制即可。`q`/`Ctrl+C` 退出，`Esc` 不退出。
  若运行的是别处已安装的二进制，用 `cargo install --path .` 或直接运行 `target/release/rswfn.exe`。

## 2026-11-25 终端分子结构预览（`rswfn view` / TUI `F4`）

- 新增 `src/viewer.rs`：把分子坐标投影成"每个字符格两个像素"的半块字符画。
  - 投影：先用坐标协方差矩阵的主轴（自带的 Jacobi 旋转，不引入线代依赖）把分子摆正，
    再按四元数 Trackball 旋转（绕屏幕轴累加、无俯仰限制）并做正交投影；自适应缩放把球半径也算进外廓，旋转时不越界。
  - 着色：球面法线做 Lambert 漫反射 + 少量高光，化学键按圆柱着色并裁剪到球面，
    远处原子向背景雾化（景深）；元素颜色取 `datas::elems` 里的 CPK 十六进制值。
  - 遮挡：画家算法按深度排序，键的深度取两端平均值。
  - 输出：`Frame::lines()` 给 ratatui，`Frame::ansi()` 给 CLI，`Frame::plain()`
    给无色终端或重定向输出，三者共用同一份像素数据。
- 新增 CLI 子命令 `rswfn view <file>`，参数：`--size 宽x高`、`--style ball|space|wire`、
  `--rotate yaw,pitch`、`--zoom`、`--no-labels`、`--no-depth`、`--color auto|always|never`；
  `--json` 输出原子坐标（Å）与成键关系。
- TUI：`F4`（或菜单/结果页的 `p`）打开全屏结构预览；方向键旋转、`+/-` 缩放、
  `s` 切换球棍/空间填充/线框、`t` 元素标签、`d` 景深、`1/2/3` 正视/侧视/俯视、`r` 复位、
  `Esc` 关闭预览（不会退出程序）。弹层只覆盖内容区，不再压住顶部信息条和底部状态栏。
- 旋转改为四元数 Trackball：`View` 不再保存 yaw/pitch，而是累积四元数；方向键/`hjkl`
  绕当前屏幕轴旋转，越过 90° 后可继续翻转（参考 three.js TrackballControls），不会再卡在极点。
  新增 `trackball_rotation_passes_poles`、`quat_rotation_preserves_length` 测试。
- 自适应缩放改用旋转不变的包围球（质心到原子 + 绘制半径），不再用投影后的二维包围盒。
  之前横向旋转时投影宽度变化会带着整体比例变，看起来忽大忽小；现在 `fit_scale` 只依赖包围球、
  画布尺寸和 zoom，旋转时大小完全稳定。新增 `auto_scale_is_rotation_invariant` 测试。
- 读取分子信息的后台线程顺带缓存 `viewer::Molecule` 快照，打开预览不再二次读文件；
  `Header::from_path` 收敛了这段读取逻辑。
- `src/base/atoms.rs`：新增 `Atoms::bonds()`，返回构造时按距离判据生成的成键列表。
- 测试：14 个 viewer 单元测试（投影不越界、元素标签、真彩色转义、缩放/旋转/模型差异、
  从 `moles/C6H6.fch` 读取）+ TUI 预览窗口渲染测试；设置 `RSWFN_SHOT=<路径>` 时
  `tui::tests::dumps_viewer_screenshot` 会把真实终端缓冲导出成 ANSI，供文档截图用。

## 2026-11-25 CLI/TUI 结果结构化与 JSON 输出

- 新增 `src/report.rs`：`Report` / `Table` / `Column` / `Field` 结果模型。文本渲染实现
  `Display`（按 CJK 显示宽度对齐），并支持 `to_json()` / `to_json_compact()`。
- 新增 `src/service/`：`atomprop`、`bondprop`、`moleprop`、`orbtprop`、`gridprop`、
  `qmcal`、`file` 任务层，统一返回 `Result<Report, String>`；CLI 与 TUI 共用同一份实现。
- `src/cli.rs`：原有 `print_*` 与计算逻辑移入 service；CLI 只负责参数映射与输出。
  新增全局 `--json` 开关（`rswfn <命令> ... --json` 或 `rswfn --json <命令> ...`）。
- `src/tui.rs`：改用 service，修复 `input()` 不打印提示导致的“假死”；加载一次分子后可
  反复分析，计算出错只提示不退出。
- 计算过程中的 SCF 迭代、拓扑格点、读取进度等输出从 stdout 改到 stderr，保证 `--json`
  的 stdout 只包含 JSON（`qmcal/hf.rs`、`qmcal/ks.rs`、`props/gridprop/topolog.rs`、
  `props/orbtprop/obtmat.rs`、`qmcal/dhmo.rs`、`reader/{mde,out}.rs`、`base/atoms.rs`）。
- 可预期错误改为返回错误信息并以退出码 1 结束（文件不存在/类型不支持/原子或轨道序号越界/
  缺少基组或轨道系数等）。
- 新增 `Report` 的 JSON 往返与文本渲染单元测试。

说明：

- 文本输出格式有小幅调整（表格宽度自适应、总电子数等移入字段区），数值本身不变。
- `Report` 目前是二进制端的模块；若后续需要给库使用者或 Python 绑定复用，可平移到
  `src/lib.rs` 导出。
- `[profile.release] panic` 当时仍为 `abort`（后续改为 `unwind`，见下一节）；计算器内部的
  `panic!` 尚未全部改为 `Result`。

## 2026-11-25 全屏 TUI（ratatui）

- `src/tui.rs` 重写为全屏界面：顶部分子信息条，左侧任务菜单，右侧欢迎页/参数表单/计算进度/
  结果视图，底部状态栏与按键提示。
- 表单由 `field_specs()` 元数据生成，覆盖 CLI 的全部 30 个分析/文件任务；必填文件、数字、
  选项在提交前校验，错误以弹窗显示。
- 计算在后台线程执行，`catch_unwind` 兜底；`Esc` 可把长任务转入后台继续跑，完成后按 `F3`
  查看结果，`F2` 查看运行日志（SCF 迭代等）。
- 新增文件浏览器（`o` / 表单内 `Ctrl+O`）、帮助（`F1`）、日志（`F2`）、错误弹窗；结果可导出
  文本 / JSON（`e` / `j`）。
- 进度输出统一改用 `log`：CLI 安装 stderr 日志器（stdout 保持纯结果/JSON），TUI 收集到内存
  日志面板，避免后台线程写屏幕破坏界面。
- `Cargo.toml`：新增 `ratatui = "0.29"`；release profile 从 `panic = "abort"` 改为 `unwind`，
  否则 worker 线程的 panic 会直接终止整个进程。
- `src/main.rs`：无参数启动时调用 `tui::run()`。
- 新增测试：文本输入编辑、表单默认值与文件校验、任务表唯一性，以及用 `TestBackend` 对欢迎页
  和全部 30 个表单做无终端渲染。

说明：

- 进程内计算无法真正中止，`Esc` 只是把任务转入后台；若要彻底取消，需要给 SCF/格点循环加取消令牌。
- TUI 与 CLI 共用 `src/service/`，新增分析任务时只需在 `TASKS` / `field_specs()` / `build_task()`
  各加一条，即可同时出现在命令线和界面中。

## 2026-11-25 release 构建崩溃排查

- 现象：`cargo build --release` 编译 bin 时 rustc 以 `0xc0000005`（STATUS_ACCESS_VIOLATION）
  退出。这是**编译器/LLVM 崩溃**，不是代码错误；同一条命令重试可以成功，属于非确定性故障。
- 原因判断：lib 与 bin 同名 `rswfn`，release 同时开了 `debug = true` 和 `strip = true`，两者都会写
  `target/release/deps/rswfn.pdb`；在 `lto = true` + `codegen-units = 1` 下，这个 PDB 同名冲突
  容易表现为 rustc 访问违例（构建日志里一直有 “output filename collision ... rswfn.pdb” 警告）。
  而且 `strip = true` 本来就会剥掉调试信息，`debug = true` 没有意义。
- 处理：release profile 改为 `debug = false`，消除 PDB 同名冲突；重新构建成功，PDB 警告消失，
  产出的 `rswfn.exe`（8.8 MB）运行正常。
- 备选：若仍偶发崩溃，把 `lto = true` 改为 `lto = "thin"`、`codegen-units = 16`，构建时间会显著
  下降（fat LTO + 单 CGU 约 17–20 分钟），性能差异通常很小。
- 需要性能剖析时：临时改为 `debug = true` + `strip = false`，但要注意 lib/bin 同名 PDB 的问题。

## 2026-11-25 TUI：文件列表、CLI 层级树、弹出层边框

- 左侧新增文件列表面板（与 `o` 弹出的浏览器共用 `FilePanel`）：显示目录与支持的格式，
  `Enter` 打开文件/进入目录，`Backspace` 上级，当前文件高亮。
- 分析菜单改为真正的树：层级严格对应 `src/cli.rs` 的子命令结构
  （atomprop → charge → mulliken/lowdin/…；bondprop → order → …；qmcal → hf/dft 等）。
  `←/→` 或 `h/l` 折叠/展开，`Enter/空格` 展开或打开任务，默认折叠。
- 每个任务增加 `cli` 字段，表单底部和状态栏会显示对应的命令行，例如
  `rswfn atomprop charge mulliken <file>`。
- `Tab` 在文件列表与分析菜单之间切换，当前面板边框高亮。
- 弹出层（文件浏览器）改为带边框的窗口，内容放在边框内部；帮助/日志/错误原本已有边框。
- `rswfn <文件>` / `rswfn <目录>` / `rswfn --tui [路径]` 直接进入 TUI；目录作为文件列表根目录。
  其他参数仍走 CLI；`--help`、`--json` 等以 `-` 开头的参数不受影响。
- 新增测试：树覆盖全部任务、折叠/展开、文件列表、初始目录/文件、浏览器边框、文件面板渲染。
- `view`（`src/viewer.rs`，F4/p 结构预览）也作为树中的顶级节点「结构预览（view）」；Enter 打开
  预览，状态栏显示 `rswfn view <file>`。该模块与 `E:\code\rswfn-viewer` 工作树逐字节一致，
  14 个 viewer 测试全部通过。
- 菜单中 `Esc` 不再退出程序（只提示退出方式），退出用 `q` 或 `Ctrl+C`；弹窗/表单/结果的
  `Esc` 仍只关闭或返回。新增 `esc_in_menu_does_not_quit` 测试。
- 工作树 `E:\code\rswfn-viewer` 里的 `demo/`（PNG 截图 + `render.py`）没有并入主仓库；
  README/notes 未引用这些图片。
