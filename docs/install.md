# rswfn 安装指南

`rswfn` 有两种安装途径：

- **预编译二进制**：Windows 使用 `install.ps1`，Linux 使用 `install.sh`，一条命令从
  GitHub Releases 下载并安装，不需要 Rust 工具链，也不需要管理员权限。
- **从源码编译**：用 `cargo install rswfn` 从 crates.io 安装（需要 Rust 工具链）。

> **版本说明**：预编译二进制来自 GitHub Releases。如果仓库还没有任何 Release
> （例如刚引入安装脚本、尚未发布首个版本），两条脚本命令会因为找不到 Release
> 而失败，请先用 `cargo install rswfn`；首个 Release（`v0.1.x`）发布后即可使用
> 两条一行命令。`rswfn` 0.1.0 已发布到 crates.io，可以随时 `cargo install`。

## 1. 安装方式对比

| 方式 | 需要 Rust 工具链 | 支持平台 | 说明 |
| --- | --- | --- | --- |
| `install.sh` / `install.ps1` | 否 | Windows x86_64；Linux x86_64（默认全静态 musl，可选 gnu）、aarch64（gnu，glibc ≥ 2.39） | 下载预编译二进制，校验 SHA256 后安装 |
| `cargo install rswfn` | 是 | Cargo 支持的所有平台（含 macOS、Windows ARM64） | 从 crates.io 拉源码编译，耗时较长；当前版本 0.1.0 |
| 源码构建 | 是 | 同上 | 克隆仓库后 `cargo build --release --locked` |

Release 中提供的资产：

| 资产 | 平台 | 备注 |
| --- | --- | --- |
| `rswfn-<version>-x86_64-pc-windows-msvc.zip` | Windows x86_64 | 使用静态 CRT，单个可执行文件，不依赖 `vcruntime140.dll` |
| `rswfn-<version>-x86_64-unknown-linux-musl.tar.gz` | Linux x86_64 | 全静态链接，任意发行版/容器可运行；安装脚本默认选它 |
| `rswfn-<version>-x86_64-unknown-linux-gnu.tar.gz` | Linux x86_64 | 动态链接 glibc，要求 glibc ≥ 2.39；用 `--target` 指定 |
| `rswfn-<version>-aarch64-unknown-linux-gnu.tar.gz` | Linux aarch64 | 动态链接 glibc，要求 glibc ≥ 2.39；暂无 aarch64 musl 产物 |
| `SHA256SUMS` | — | 每个归档一行 `<哈希><两个空格><文件名>` |

没有预编译产物的平台（macOS、Windows ARM64、Linux aarch64 且 glibc < 2.39 等）请使用
`cargo install rswfn` 或源码构建。安装脚本在 macOS 上会直接报错并给出该提示。

## 2. Linux：一键安装

```sh
curl -fsSL https://raw.githubusercontent.com/qmlCoder/rswfn/main/install.sh | sh
```

脚本会：

- 识别 CPU 架构，并据此选择 target（不需要探测系统 libc）：
  - `x86_64` / `amd64` → **`x86_64-unknown-linux-musl`**（全静态，任何发行版或容器都能
    运行，包括 Alpine、CentOS 7、Debian 11 等，不依赖 glibc 版本）；
  - `aarch64` / `arm64` → `aarch64-unknown-linux-gnu`（暂无 musl 产物，要求 glibc ≥ 2.39）；
  - 其他架构直接报错，并提示用 `--target` 指定目标平台。
- 从 GitHub Releases 下载对应归档，并与 `SHA256SUMS` 校验（见第 6 节）；
- 安装到 `~/.local/bin`；`$HOME` 未设置或不可写时回退 `/usr/local/bin`（仅在该目录可写时），
  两者都不行则报错并提示 `--dir`；**不会调用 sudo**；
- **不会修改任何 shell 配置文件**（`~/.bashrc`、`~/.zshrc`、`~/.profile` 等都不动）；
  如果安装目录不在 `PATH` 中，会打印需要追加的 `export PATH=...` 行（见第 9 节）。

安装完成后执行 `rswfn --version` 验证。

### 可选：安装动态链接的 gnu 版本

追求一点点性能、且系统 glibc ≥ 2.39 时，可以改为安装 gnu 产物：

```sh
curl -fsSL https://raw.githubusercontent.com/qmlCoder/rswfn/main/install.sh | sh -s -- --target x86_64-unknown-linux-gnu
```

如果指定的 target 没有对应资产（例如 aarch64 + musl），脚本会报 404 并提示改用其他
`--target`。想回到默认的静态版，去掉 `--target` 重新执行一键命令即可覆盖安装
（如果之前用 `RSWFN_TARGET` 指定过，记得同时清除该环境变量）。

### 先审阅再执行

```sh
curl -fsSL https://raw.githubusercontent.com/qmlCoder/rswfn/main/install.sh -o install.sh
# 查看 install.sh 的内容
sh install.sh
```

### 传递参数

以 `curl ... | sh` 方式运行时，参数要写在 `sh -s --` 后面：

```sh
curl -fsSL https://raw.githubusercontent.com/qmlCoder/rswfn/main/install.sh | sh -s -- --version 0.1.1
curl -fsSL https://raw.githubusercontent.com/qmlCoder/rswfn/main/install.sh | sh -s -- --dir "$HOME/bin"
curl -fsSL https://raw.githubusercontent.com/qmlCoder/rswfn/main/install.sh | sh -s -- --uninstall
```

也可以先下载脚本，再带参数反复运行（两者等价）：

```sh
curl -fsSL https://raw.githubusercontent.com/qmlCoder/rswfn/main/install.sh -o install.sh
sh install.sh --version 0.1.1
sh install.sh --dir "$HOME/bin"
sh install.sh --uninstall
```

## 3. Windows：一键安装

在 PowerShell 中运行：

```powershell
irm https://raw.githubusercontent.com/qmlCoder/rswfn/main/install.ps1 | iex
```

脚本会：

- 使用 `x86_64-pc-windows-msvc` 产物（目前只提供 x86_64；ARM64 会明确报错）。产物使用
  **静态 CRT**，是单个 `rswfn.exe`，不需要额外安装 `vcruntime140.dll` 等运行库；
- 下载归档并与 `SHA256SUMS` 校验（见第 6 节）；
- 安装到 `%LOCALAPPDATA%\Programs\rswfn`（`LOCALAPPDATA` 缺失时回退
  `%USERPROFILE%\AppData\Local\Programs\rswfn`，不存在则创建）；
- 默认把安装目录追加到**用户级** PATH（已存在则不重复添加），不写系统 PATH、不需要
  管理员权限；加 `-NoPath` 可跳过这一步；
- 安装完成后请**重开终端**，再执行 `rswfn --version` 验证。

`irm ... | iex` 不能直接传参数。需要指定版本、目录或卸载时，先下载脚本再用
`powershell -File` 运行：

```powershell
irm https://raw.githubusercontent.com/qmlCoder/rswfn/main/install.ps1 -OutFile install.ps1
powershell -ExecutionPolicy Bypass -File .\install.ps1 -Version 0.1.1
powershell -ExecutionPolicy Bypass -File .\install.ps1 -NoPath
powershell -ExecutionPolicy Bypass -File .\install.ps1 -Uninstall
```

> 如果直接运行 `.\install.ps1 ...` 被 PowerShell 执行策略拦截，上面
> `-ExecutionPolicy Bypass` 的写法只对当前这次进程生效，不会修改系统设置。
> 脚本兼容 Windows PowerShell 5.1 和 PowerShell 7。

## 4. 参数与环境变量

两个脚本语义一致：**命令行 > 环境变量 > 默认值**。

| 用途 | Linux（`install.sh`） | Windows（`install.ps1`） | 环境变量 |
| --- | --- | --- | --- |
| 显示帮助 | `-h`、`--help` | `-h`、`-Help` | — |
| 指定版本 | `--version <x.y.z>`（可带 `v` 前缀） | `-Version <x.y.z>` | `RSWFN_VERSION` |
| 安装目录 | `--dir <path>` | `-Dir <path>` | `RSWFN_INSTALL_DIR` |
| 下载基址 | `--base-url <url\|path>` | `-BaseUrl <url\|path>` | `RSWFN_BASE_URL` |
| 跳过 SHA256 校验 | `--skip-verify` | `-SkipVerify` | `RSWFN_SKIP_VERIFY=1` |
| 指定 target | `--target <triple>` | `-Target <triple>` | `RSWFN_TARGET` |
| 卸载 | `--uninstall` | `-Uninstall` | — |
| 不写 PATH | 不适用（sh 从不修改 rc 文件） | `-NoPath` | — |

- 版本号可以省略。省略时脚本查询 GitHub API 的最新 release；查询失败且未显式指定
  版本时会报错，并提示用 `--version x.y.z` / `-Version x.y.z`。
- `--target` / `-Target` 用于强制指定目标平台，如 `x86_64-unknown-linux-gnu`、
  `x86_64-unknown-linux-musl`。
- 退出码：成功为 `0`，失败为非 `0`（`install.sh` 的用法错误为 `2`；`install.ps1`
  以抛错方式失败，`irm | iex` 场景下不会关闭当前终端）。

环境变量示例：

```sh
curl -fsSL https://raw.githubusercontent.com/qmlCoder/rswfn/main/install.sh | RSWFN_VERSION=0.1.1 sh
curl -fsSL https://raw.githubusercontent.com/qmlCoder/rswfn/main/install.sh | RSWFN_INSTALL_DIR="$HOME/bin" sh
curl -fsSL https://raw.githubusercontent.com/qmlCoder/rswfn/main/install.sh | RSWFN_SKIP_VERIFY=1 sh
```

注意环境变量要写在 `sh` 一侧（`curl ... | RSWFN_VERSION=0.1.1 sh`）；写在 `curl` 前面
只会影响 curl 进程，不会传给安装脚本。

Windows：

```powershell
$env:RSWFN_VERSION = "0.1.1"
$env:RSWFN_INSTALL_DIR = "D:\tools\rswfn"
irm https://raw.githubusercontent.com/qmlCoder/rswfn/main/install.ps1 -OutFile install.ps1
powershell -ExecutionPolicy Bypass -File .\install.ps1
```

## 5. 离线安装与镜像

下载地址规则是 `${BASE}/v${VERSION}/${ASSET}`，`BASE` 默认为
`https://github.com/qmlCoder/rswfn/releases/download`。用 `--base-url` / `-BaseUrl`
（或环境变量 `RSWFN_BASE_URL`）把它替换成镜像 URL 或本地目录即可。

`BASE` 支持三种形式：网络 URL（`https://…` / `http://…`）、`file:///…` 本地路径，
以及普通目录路径（本地目录，用于离线安装、镜像和 CI 自测）。

本地或内网镜像的目录结构要和 GitHub Release 保持一致：

```text
rswfn-mirror/
└── v0.1.1/
    ├── SHA256SUMS
    ├── rswfn-0.1.1-aarch64-unknown-linux-gnu.tar.gz
    ├── rswfn-0.1.1-x86_64-pc-windows-msvc.zip
    ├── rswfn-0.1.1-x86_64-unknown-linux-gnu.tar.gz
    └── rswfn-0.1.1-x86_64-unknown-linux-musl.tar.gz
```

Linux（本地目录、`file://` 或 HTTP(S) 镜像都支持）：

```sh
sh install.sh --version 0.1.1 --base-url /mnt/usb/rswfn-mirror
sh install.sh --version 0.1.1 --base-url file:///mnt/usb/rswfn-mirror
sh install.sh --version 0.1.1 --base-url https://mirror.example.com/rswfn/releases/download
```

Windows：

```powershell
powershell -ExecutionPolicy Bypass -File .\install.ps1 -Version 0.1.1 -BaseUrl D:\rswfn-mirror
powershell -ExecutionPolicy Bypass -File .\install.ps1 -Version 0.1.1 -BaseUrl file:///D:/rswfn-mirror
```

使用本地目录（含 `file://` 形式）作为下载基址时，脚本无法从中查询最新版本，
**必须显式指定版本**（`--version` / `-Version` 或 `RSWFN_VERSION`）；完全离线的环境同理。

## 6. 校验机制（SHA256SUMS）

每个 Release 都附带 `SHA256SUMS`，其中每行形如
`<64 位小写十六进制哈希><两个空格><资产文件名>`。脚本下载归档后会与之比对：

- 不匹配：立即失败、不安装，并打印期望值与实际值；
- Release 中缺少 `SHA256SUMS`（老 Release）：打印警告后继续安装；
- `--skip-verify` / `-SkipVerify` / `RSWFN_SKIP_VERIFY=1` 可跳过校验，**不建议**使用。

手动校验（Linux）：

```sh
version=0.1.1
asset=rswfn-${version}-x86_64-unknown-linux-musl.tar.gz
base=https://github.com/qmlCoder/rswfn/releases/download/v${version}
curl -fsSLO "${base}/${asset}"
curl -fsSLO "${base}/SHA256SUMS"
sha256sum -c SHA256SUMS --ignore-missing
```

`--ignore-missing` 会忽略 `SHA256SUMS` 中未下载的其他平台文件；输出出现
`rswfn-0.1.1-x86_64-unknown-linux-musl.tar.gz: OK` 即校验通过。

Windows：

```powershell
$version = "0.1.1"
$asset = "rswfn-$version-x86_64-pc-windows-msvc.zip"
$base = "https://github.com/qmlCoder/rswfn/releases/download/v$version"
Invoke-WebRequest "$base/$asset" -OutFile $asset -UseBasicParsing
Invoke-WebRequest "$base/SHA256SUMS" -OutFile SHA256SUMS -UseBasicParsing
(Get-FileHash $asset -Algorithm SHA256).Hash.ToLower()
Select-String $asset SHA256SUMS
```

两行中的哈希一致即校验通过。也可以使用
`certutil -hashfile .\rswfn-0.1.1-x86_64-pc-windows-msvc.zip SHA256` 得到同样结果。

如果校验失败，通常是下载中断或经过了会篡改内容的代理/镜像：重新下载即可；仍失败时
不要使用 `--skip-verify`，改用官方直链或 `cargo install rswfn`。

## 7. 升级

重复执行同一条安装命令即为原地升级，安装器会原子替换旧二进制，不残留临时文件：

```sh
curl -fsSL https://raw.githubusercontent.com/qmlCoder/rswfn/main/install.sh | sh
```

升级到指定版本：

```sh
curl -fsSL https://raw.githubusercontent.com/qmlCoder/rswfn/main/install.sh | sh -s -- --version 0.1.1
```

Windows 同理（重新运行 `irm ... | iex`，或 `powershell -ExecutionPolicy Bypass -File .\install.ps1 -Version 0.1.1`）。
升级后用 `rswfn --version` 确认版本。

## 8. 卸载

Linux：

```sh
sh install.sh --uninstall
```

安装时用了 `--dir` 或 `RSWFN_INSTALL_DIR` 的话，卸载时也要给出同一个目录：

```sh
sh install.sh --dir "$HOME/bin" --uninstall
```

脚本只删除 `<安装目录>/rswfn`，不会动其他文件。

Windows：

```powershell
powershell -ExecutionPolicy Bypass -File .\install.ps1 -Uninstall
```

`-Uninstall` 会删除 `rswfn.exe`，目录为空时一并删除目录；如果安装目录就是默认目录，
还会从用户 PATH 中移除该目录条目。用 `-Dir` 自定义目录安装的，卸载后请自行检查 PATH。

手动卸载：Linux 删除 `~/.local/bin/rswfn`；Windows 删除
`%LOCALAPPDATA%\Programs\rswfn` 目录。

## 9. PATH 处理

Linux：脚本**不会修改** `~/.bashrc`、`~/.zshrc`、`~/.profile` 等文件，只会在安装目录
不在 `PATH` 中时打印提示。把下面一行加到对应 shell 的启动文件即可（以 bash 为例）：

```sh
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc
```

只想在当前终端临时生效：`export PATH="$HOME/.local/bin:$PATH"`。用 `command -v rswfn`
检查是否生效。Linux 上安装目录也可能是 `/usr/local/bin`（`$HOME` 不可用时的回退），
该目录通常已在 `PATH` 中。

Windows：默认写入**用户级** PATH，需要重开终端才会生效；不想写 PATH 时加 `-NoPath`，
之后用完整路径或自行添加。用 `where.exe rswfn` 检查是否生效。

## 10. 常见问题

### `rswfn: command not found`

1. 确认已安装：Linux 查看 `~/.local/bin/rswfn`，Windows 查看
   `%LOCALAPPDATA%\Programs\rswfn\rswfn.exe`。
2. 检查 PATH：Linux `command -v rswfn`，Windows `where.exe rswfn`。
3. Linux 把 `~/.local/bin` 加入 PATH（见第 9 节）；Windows 重开终端，或确认安装时是否
   用了 `-NoPath`。

### Linux 上出现 `GLIBC_2.xx not found`

默认安装的 x86_64 musl 产物是**全静态**的，不依赖系统 glibc，正常情况下不会出现这个
错误。出现该错误说明装的是 gnu 产物，可能原因：

- 你曾用 `--target x86_64-unknown-linux-gnu` 或 `RSWFN_TARGET` 指定了 gnu 版本，而系统
  glibc 低于 2.39。去掉 `--target`（并清除 `RSWFN_TARGET`）重新执行一键命令，即可覆盖
  安装为默认的静态 musl 版本；
- 你在 aarch64 上，而该架构目前只有 gnu 产物（要求 glibc ≥ 2.39）。系统较旧时请用
  `cargo install rswfn` 或源码构建。

### macOS / Windows ARM64 / 其他平台

目前只提供 Windows x86_64 和 Linux（x86_64 musl/gnu、aarch64 gnu）的预编译产物。macOS
与其他架构请用 `cargo install rswfn`（需要 Rust 工具链）；Windows ARM64 上安装脚本会
明确报错；Linux aarch64 且 glibc < 2.39 的系统需要从源码构建。

### 公司代理 / GitHub 访问不稳定

- Linux 的 `curl` / `wget` 遵循 `HTTPS_PROXY` / `HTTP_PROXY` 环境变量，在终端里设置好
  代理后再执行安装命令。
- Windows 上脚本通过 `Invoke-WebRequest` 下载，走系统代理（系统设置「代理」中的配置）。
- GitHub 下载不稳定时：用 `--base-url` / `-BaseUrl` 指向内网镜像或本地离线包（见第 5
  节），或改用 `cargo install rswfn`。

### 没有管理员权限 / 不想装到系统目录

默认安装目录都在用户目录下（Linux `~/.local/bin`、Windows
`%LOCALAPPDATA%\Programs\rswfn`），不需要管理员权限。要装到别处时，用 `--dir` / `-Dir`
指定一个当前用户可写的目录即可。

### 安装脚本提示无法确定版本

脚本会查询 GitHub API 的最新 release。网络受限或仓库还没有任何 release 时查询会失败，
此时需要显式指定版本：`--version 0.1.1` / `-Version 0.1.1`。在首个 release 发布之前，
请使用 `cargo install rswfn`。

### 想先检查脚本再运行

```sh
curl -fsSL https://raw.githubusercontent.com/qmlCoder/rswfn/main/install.sh -o install.sh
# 阅读 install.sh 后
sh install.sh
```

Windows 同理：先 `irm ... -OutFile install.ps1`，再用 `powershell -ExecutionPolicy Bypass -File .\install.ps1` 运行。

## 11. 维护者发布清单

安装脚本依赖 GitHub Releases 上的资产与 `SHA256SUMS`，发布流程如下：

1. **改版本号**：修改 `Cargo.toml` 里的 `version`（首个 release 建议 `0.1.1`，因为
   `0.1.0` 已发布到 crates.io），并同步更新 `Cargo.lock`（跑一次 `cargo check` 或
   `cargo build --release` 即可；release workflow 用 `--locked` 构建）。
2. **提交并打 tag**：tag 必须与 `Cargo.toml` 的版本一致，形如 `v0.1.1`，否则 workflow
   会报错退出。

   ```sh
   git add Cargo.toml Cargo.lock
   git commit -m "release v0.1.1"
   git push origin main
   git tag v0.1.1
   git push origin v0.1.1
   ```

3. **等待 Actions 完成**：`.github/workflows/release.yml` 会构建 4 个产物并创建
   Release，附上 `SHA256SUMS`：

   ```text
   rswfn-0.1.1-x86_64-pc-windows-msvc.zip
   rswfn-0.1.1-x86_64-unknown-linux-gnu.tar.gz
   rswfn-0.1.1-x86_64-unknown-linux-musl.tar.gz
   rswfn-0.1.1-aarch64-unknown-linux-gnu.tar.gz
   SHA256SUMS
   ```

   可用 `gh release view v0.1.1 --repo qmlCoder/rswfn` 确认资产齐全。
4. **在干净机器上验证两条一行命令**：准备一台 Linux 和一台 Windows（都不预装 Rust），
   分别运行第 2、3 节的一行命令，再确认 `rswfn --version` 输出预期版本。
5. **（可选）同步 crates.io**：GitHub Release 与 crates.io 发布相互独立；如果希望
   `cargo install rswfn` 也能装到新版本，需要另外执行 `cargo publish`。

注意：在 Actions 页面手动触发（`workflow_dispatch`）只会构建产物，不会创建 Release，
适合在正式打 tag 前验证构建是否通过。
