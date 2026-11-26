# 变更记录

## 2026-11-26 预览兼容不支持真彩色的终端（Xshell）

- 现象：在 Xshell 里结构预览整屏变成随机彩色横条，看不出分子。
- 原因：`Frame::ansi()` / `Frame::lines()` 无条件输出 24 位真彩色 SGR
  （`38;2;r;g;b;48;2;r;g;b`），并做了游程压缩；终端不支持真彩色时这些参数会被当成
  普通 SGR / 256 色索引重新解释，于是每条色带/色块都变成随机的调色板颜色。
  TUI 主题用的是命名色（Cyan/Green/…），所以只有预览区域崩坏。
- 新增 `viewer::Depth`（True / Ansi256 / Mono）与 `Depth::detect()`：读
  `COLORTERM`、`TERM`、`TERM_PROGRAM`、`WT_SESSION` 判断能力，TERM 为空时按平台兜底。
- 新增 `Rgb::ansi256()`（6×6×6 色立方 + 24 级灰阶取最近色）与 `Rgb::color_in(depth)`。
- `Frame` 新增 `ansi_in(depth)` / `lines_in(depth)`；`ansi()` / `lines()` 保留为真彩色版本。
  Mono 模式复用 `plain()` 的 ASCII 灰度斜坡（抽出为 `mono_char`）。
- CLI `run()` 根据检测结果选择：真彩色 → `ansi_in(True)`，256 色 → `ansi_in(Ansi256)`，
  无色 → `plain()`。
- TUI 预览：新增 `c` 键在 真彩色 → 256 色 → ASCII 之间循环，弹层底色与页脚同步跟随；
  打开预览时默认用 `Depth::detect()`。
- 验证：`TERM=xterm-256color` 时输出 12 行 `38;5;`、0 处 `38;2;`；`TERM=dumb` 时 0 个
  转义序列；新增 `ansi256_quantizes_to_palette` 测试，75 项全部通过。

## 2026-11-26 OSS 一行安装（latest.txt + 镜像版安装脚本）

- `install.sh` / `install.ps1`：新增「镜像 `latest.txt` 解析最新版本」。
  - `install.sh`：`resolve_latest_from_mirror()`，remote 模式下先取 `${BASE_URL}/latest.txt`
    （一行版本号），取不到再退回 GitHub 的 releases/latest。
  - `install.ps1`：对应 `Get-MirrorLatestVersion`（`$remoteBase/latest.txt`）。
  - 只影响「未显式指定版本」的情况，`--version` / `-Version` 行为不变，CI 的安装器
    自测（本地目录 + 显式版本）不受影响。
- `.github/scripts/upload_oss.py`：
  - 新增 `--installers`：把仓库里的 `install.sh` / `install.ps1` 改写成「默认下载基址
    指向本镜像」后上传到 `{prefix}/install.sh`、`{prefix}/install.ps1`（同时改写脚本内
    的示例 raw.githubusercontent URL；sh 强制 LF、ps1 强制 CRLF）。
  - 新增 `--latest`：上传 `{prefix}/latest.txt`（一行版本号）。
  - 版本目录用 immutable 长缓存，根部文件用 `max-age=60`，避免 `latest.txt` 被缓存住。
  - 位置参数改为 `<version> [dist-dir]`：dist 可省，只发 installers / latest.txt。
- `scripts/release.sh` 与 `.github/workflows/release.yml` 调用时加 `--installers --latest`；
  `mirror-oss.yml` 刻意只补版本目录（回填旧版本不应把 `latest.txt` 退回）。
- 已实测（两个平台都跑通）：
  - Linux（WSL）：`curl -fsSL https://wmview.oss-cn-hangzhou.aliyuncs.com/rswfn/install.sh | sh`
    → 解析到 0.1.3、下载 musl 包、SHA256 校验、自检通过。
  - Windows：同源 `install.ps1` 同样解析到 0.1.3、下载 zip、校验、自检、安装。
- 记录一个与本改动无关的环境坑：从 Zed/Git Bash 调 `powershell`（5.1）时，
  `PSModulePath` 继承了 PS7 的 `C:\Program Files\PowerShell\7\Modules`，PS 5.1 会先
  命中版本不兼容的模块，导致 `Get-FileHash` / `Expand-Archive` 报「无法识别为 cmdlet」。
  用标准 PSModulePath 运行即正常；正常开一个 PowerShell 窗口不受影响。

## 2026-11-26 局部坐标系 π-only Kohn-Sham（`qmcal lks`）

- 新增 `qmcal::lks`（`src/qmcal/lks.rs`）：`LKS` / `LKSResult`。
  - 为每个原子构造局部坐标系（`atomprop::direction::Calculator::LCS`，无坐标系时退回 `RLCS`）。
  - 用 `RsInt::mat1e_lcs` / `mat2e_lcs` 在局部坐标系（笛卡尔）下计算单/双电子积分，
    再按「局部 x 指数为奇数」挑出 π 成分（`p_x`、`d_xy`、`d_xz`、`f_xxx`、`f_xyz`、…），
    截取得到 `nπ × nπ` 的 S/H/双电子积分；Fock、DIIS、对角化、格点 XC 都在 `nπ` 维上进行。
  - 格点上的 π 原子轨道取值/梯度由全局笛卡尔 AO 值经各壳层 `get_tmat_car` 旋转后取 π 列，
    与 `mat1e_lcs` 的变换约定（`v^l = T·v^g`）一致。
  - π 电子数：`LKS::npi_elec` 可显式指定；缺省按「拥有 p 壳层的原子数」估计（乙烯→2、苯→6），
    对吡咯等含孤对杂原子体系可能偏少，应显式指定。
  - `LKSResult::expand_coefs` 将 π 系数展开回全局笛卡尔基（`C_glob = P·C_π`）。
  - `LKSResult::local_labels` / `mo_table_local`：打印 **局部坐标系下** 的 π 分子轨道系数表，
    行标签只含 π 成分（如 `C1:PX`、`C3:XY`），列为 π 分子轨道并附轨道能。
  - `LKSResult::to_mole`：构造一个「全局笛卡尔基 + π 占据数」的 `Mole`（系数为展开后的正常
    全局轨道系数，`nπ` 个 MO，电子数=π 电子数），供 `FchWriter` / `MfnWriter` 导出。
- CLI 新增 `rswfn qmcal lks -i <file> [-b <basis>] [-f <func>] [--npi-elec N] [-o out.fch|out.mwfn] [--print-mo]`
  （`src/cli.rs` 的 `QmcalCmds::Lks`、`src/service/qmcal.rs::lks`），输出 π 约化基大小、
  π 电子数、π 电子能量、π 分子轨道能量与耗时；`--print-mo` 打印局部 π 系数表，
  `-o` 导出 `.fch/.fchk` 或 `.mwfn`（系数已展开回全局基，便于 Multiwfn 等可视化程序查看）。
- 新增 `writer::MfnWriter`（`src/writer/mfn.rs`）：生成与 `MfnReader` 兼容的 Multiwfn 波函数文件；
  `Nbasis` 为基函数数、`Nindbasis` 为轨道数（π-only 时前者大于后者）。
- 注意导出文件描述的是 **π 子体系**：`Nindbasis = nπ < Nbasis`，电子数写作 π 电子数，
  故与核电荷总数不自洽（仅供可视化 π 轨道）；能量写入的是 π 模型电子能量。
- 注意：这是 **π 子系统模型**，能量不含 σ 电子与核内层，*不是* 分子总能量
  （原理与坑见 `docs/pi电子能量.md`）。积分变换仍复用 `mat2e_lcs`（全维四级变换后再截断），
  因此双电子积分耗时未减少；真正变小的是 Fock/系数/DIIS/本征求解与格点 XC。
- 新增集成测试 `tests/qmcalc.rs::test_lks_c2h4`（乙烯 B3LYP/6-31G(d)：nπ=8、π 电子数 2、
  SCF 收敛、`expand_coefs` 形状正确、局部系数表含 `PX` 标签、fch / mwfn 导出回读
  轨道数与 π 电子数一致）。实测苯 STO-3G 的 π 轨道能级次序与简并
  （a2u < e1g 双重简并 < e2u 双重简并 < b2g）与理论一致。

## 2026-11-26 核-核排斥能下沉到 `Atoms`

- `base::atoms::Atoms` 新增 `nuc_eng()`：按 `Σ_{i<j} Z_i Z_j / r_ij` 计算核-核排斥能，
  距离用 `dist(i, j)` 现算（与 `dmat()` 一致，避免 `add()` 后 `dists` 变陈旧）。
- `qmcal::ks::RKS::nuc_eng` 与 `qmcal::hf::RHF::nuc_eng` 改为直接调用 `self.mole.atoms.nuc_eng()`，
  删除两份重复实现；`UKS` 原已委托给 `RKS`，无需改动。
- 仍沿用原子序数 `nuc` 而非有效核电荷 `crg`，数值与改动前完全一致；
  CLI 的 `qmcal hf/dft` 已输出核排斥能，无需新增命令。

## 2026-11-26 发行改为本地交叉编译（Windows + WSL）

- 新增 `scripts/release.sh`：一条命令完成「编译 → 打包 → 分别上传 OSS 与 GitHub」。
  - Windows 目标（`x86_64-pc-windows-msvc`）本机编译，
    `RUSTFLAGS=-C target-feature=+crt-static`（与 CI 一致，不依赖 vcruntime140.dll）；
  - 3 个 Linux 目标通过 `wsl -d <distro>` 交给 `scripts/build-wsl.sh` 编译；
  - 打包结构与 CI 一致：zip 为扁平结构（`rswfn.exe` + README/LICENSE），tar.gz 带顶层目录；
  - `SHA256SUMS` 改用 Python 生成：Git Bash 的 `sha256sum` 在 Windows 上输出
    `hash *name`（星号 + 单空格），不符合 `<64hex><两空格><name>` 的契约；
  - 默认创建并推送 `v<version>` tag；`--skip-build` / `--no-oss` / `--no-github` /
    `--no-tag` / `--targets` / `--wsl` / `--wsl-target-dir` 可选。
- 新增 `scripts/build-wsl.sh`：在 WSL 里逐目标设置
  `CC_*` / `CARGO_TARGET_*_LINKER`（`musl-gcc`、`aarch64-linux-gnu-gcc`，vendored OpenSSL 也用它），
  编译后把 tar.gz 打到 `dist/`；默认 `CARGO_TARGET_DIR=~/.cache/rswfn-target`
  （Linux 原生文件系统，比 `/mnt/e` 快很多）。
- `.github/workflows/release.yml`：去掉 `push: tags` 触发，改为仅 `workflow_dispatch`，
  避免与本地发布同时写同一份 Release 资产；组件保留作为本地工具链不可用时的兜底。
- 环境侧准备（已在本机执行）：WSL 里 `rustup target add aarch64-unknown-linux-gnu` 与
  `apt install gcc-aarch64-linux-gnu`。
- 验证范围：aarch64 交叉编译（`cargo build --locked --target aarch64-unknown-linux-gnu`）
  成功，产物为 `ELF 64-bit LSB pie executable, ARM aarch64`；Windows 打包与
  `SHA256SUMS` 格式通过契约校验 `^[0-9a-f]{64}  [^/]+$`；WSL 委派链路（路径映射、
  参数与 env 传递）已验证。未跑完整的 4 目标发布（首次全量构建每目标约 10–20 分钟）。
- 两个坑：
  - Git Bash 把 `E:\code\rswfn` 传给原生程序时反斜杠被 MSYS 吃掉（`wslpath` 返回
    `E:coderswfn`），所以改用 `REPO_WSL=/mnt$PWD`；
  - 脚本必须是 LF：CRLF 会让 WSL 的 dash 在 `set -eu` 上报 `Illegal option -`，
    release.sh 调用前会先 `sed -i 's/\r//g' scripts/*.sh` 兜底（`.gitattributes`
    已有 `*.sh text eol=lf`）。

## 2026-11-26 把 OSS 凭据收到 .env

- 新增 `.env`（仓库根，已加入 `.gitignore`）记录 `OSS_ACCESS_KEY_ID` /
  `OSS_ACCESS_KEY_SECRET` / `OSS_ENDPOINT` / `OSS_BUCKET` / `OSS_PREFIX`；
  另加可提交的模板 `.env.example`。
- `.github/scripts/upload_oss.py`：新增 `load_dotenv()`，取值顺序为
  **进程环境变量 > 仓库根 `.env`**，因此 CI 里靠 GitHub Secrets 注入的环境变量
  优先，本地直接跑脚本则自动读 `.env`（支持 `export ` 前缀、引号、BOM/CRLF）。
- `Cargo.toml` 的 `exclude` 增加 `/.env` 与 `/.env.*`，防止凭据误入 crates.io 包。
- 文档：`docs/install.md` 第 11 节补充 `.env` 用法——
  `gh secret set -f .env --repo qmlCoder/rswfn` 一次性同步到 Actions secrets，
  以及本地用 `gh release download` + `upload_oss.py` 回填。
- 实测：在未导出环境变量的子进程里跑脚本，成功回填 `oss://wmview/rswfn/v0.1.3/`。

## 2026-11-26 发布资产同步到阿里云 OSS

- 新增 `.github/scripts/upload_oss.py`：用官方 `oss2` SDK 把 `dist/` 里的发布资产上传到
  `oss://wmview/rswfn/v<version>/`（4 个归档 + `SHA256SUMS`）。
  - 目录名用 `v<version>`，与 `install.sh` / `install.ps1` 的镜像约定
    `<base-url>/v<version>/<asset>` 一致，所以 OSS 地址可以直接当 `--base-url` 用。
  - 按扩展名设置 `Content-Type`，并带 `Cache-Control: public, max-age=31536000, immutable`
    （版本目录内容不会再变）。
- `.github/workflows/release.yml`：`release` job 补上 `actions/checkout`（拿脚本），
  并在创建 GitHub Release 后多一步 `Upload to Aliyun OSS`。
  - 因为 GitHub 不允许在 step 的 `if` 里直接引用 secrets，凭据放在 job 级 `env`
    （`OSS_ACCESS_KEY_ID` / `OSS_ACCESS_KEY_SECRET` / `OSS_ENDPOINT`，可选
    `OSS_BUCKET` / `OSS_PREFIX`），再用 `env` 判断：没配置就静默跳过，不影响
    GitHub Release 与分叉仓库。
- 新增 `.github/workflows/mirror-oss.yml`：`workflow_dispatch` 手动把某个已有 Release
  镜像到 OSS（填 tag，留空则取 latest），用于回填 v0.1.3 这类之前没镜像的版本；
  重跑幂等。手动触发时缺凭据直接报错，不静默跳过。
- 文档：`docs/install.md` 第 5 节新增「官方阿里云 OSS 镜像」（目录结构、公共读与
  “必须显式给版本”的限制），第 11 节发布清单加入 OSS 校验与 secrets 说明；README 补一行。

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
