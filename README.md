# rswfn

[![crates.io](https://img.shields.io/crates/v/rswfn.svg)](https://crates.io/crates/rswfn)
[![docs.rs](https://docs.rs/rswfn/badge.svg)](https://docs.rs/rswfn)
[![license](https://img.shields.io/crates/l/rswfn.svg)](#license)

**rswfn** (Rust Wavefunction) is a wavefunction analysis library and
command-line tool for quantum chemistry, written in Rust. It reads
electronic-structure data produced by common quantum chemistry programs and
performs a wide range of wavefunction analyses — population analysis, bonding
analysis, aromaticity, orbital analysis, and real-space (grid) function
analysis.

It can be used both as a **library** and as a **CLI / TUI** application.

## Features

- **Readers** for Gaussian (`.gjf`, `.log`, `.out`, `.fch`, `.fchk`), ORCA
  (`.inp`, `.out`), Molden (`.molden`), cube (`.cub`), `.mol`, `.sdf`, `.xyz`,
  and a custom `.dmi` / `.dmo` format that carries per-atom local coordinate
  frames.
- **Atomic properties**: Mulliken / Löwdin / Hirshfeld charges, spin
  population, π-electron population (the POCV and MOCV orbital-transformation
  methods), Fukui functions and dual descriptors, and free valence.
- **Bond properties**: Mayer and Wiberg bond orders, per-orbital bond-order
  contributions, σ / π / δ / φ decomposition, and π bond orders.
- **Molecular & orbital properties**: π-based aromaticity (PISD), multipole
  moments, orbital composition, and per-orbital π-electron counts.
- **Real-space / grid analysis**: electron density and its derivatives,
  sign(λ₂)ρ, RDG, IRI, ELF / LOL, kinetic-energy densities, electrostatic
  potential, information entropy, Fukui-function grids, and topology search
  for critical points.
- **Quantum-chemical calculations**: closed-shell Hartree-Fock (RHF) and
  Kohn-Sham DFT from a geometry plus basis set. The DFT functional is
  selectable (B3LYP, PBE0, PBE, BLYP, BP86, SVWN5), and the resulting
  wavefunction is exportable as `.fch`.
- **Integrals** over Gaussian basis sets: overlap, kinetic, nuclear-attraction
  and two-electron repulsion.
- **Writers** for cube (`.cub`), `.xyz`, `.gjf`, `.fch`, `.mol` / `.dmo` and
  `.obj`.
- **Terminal structure preview**: a truecolor ball-and-stick / space-filling /
  wireframe picture of the molecule, available interactively in the TUI (`F4`)
  and non-interactively as `rswfn view <file>`.

## Installation

### One-line install (prebuilt binary, no Rust toolchain needed)

Windows (PowerShell):

```powershell
irm https://raw.githubusercontent.com/qmlCoder/rswfn/main/install.ps1 | iex
```

Linux (x86_64 and aarch64):

```sh
curl -fsSL https://raw.githubusercontent.com/qmlCoder/rswfn/main/install.sh | sh
```

The scripts download a prebuilt binary from GitHub Releases, verify it against
`SHA256SUMS`, and install it:

| Platform | Archive | Default install directory |
| --- | --- | --- |
| Windows x86_64 (MSVC) | `rswfn-<version>-x86_64-pc-windows-msvc.zip` | `%LOCALAPPDATA%\Programs\rswfn` (appended to the user `PATH`) |
| Linux x86_64 | `rswfn-<version>-x86_64-unknown-linux-musl.tar.gz` (fully static) | `~/.local/bin` |
| Linux aarch64 | `rswfn-<version>-aarch64-unknown-linux-gnu.tar.gz` | `~/.local/bin` |

On Linux x86_64 the installer picks the fully static **musl** build, so it runs
on any distribution or container — Alpine, CentOS 7, Debian 11, … — regardless
of the installed glibc version. The aarch64 build is glibc-based and requires
glibc ≥ 2.39. Advanced users who prefer the dynamically linked build can install
it with `--target x86_64-unknown-linux-gnu` (slightly faster; it targets
glibc ≥ 2.39, so older systems may print a version warning there — the static
build is the safe default). The Windows binary uses the static CRT: it is a
single self-contained `.exe` with no `vcruntime140.dll` dependency.

> **Note:** the one-line installers read the release assets from GitHub
> Releases. If no release is available yet — or if you prefer building from
> source — use `cargo install rswfn` (version 0.1.3 is published on crates.io).
> Only Windows x86_64 and Linux x86_64/aarch64 have prebuilt binaries; Windows
> ARM64 and macOS have none.

Options (pin a version, choose the install directory, force a Linux target,
uninstall, or install from a mirror):

```sh
# Linux: arguments after `sh -s --` are passed to the installer
base=https://raw.githubusercontent.com/qmlCoder/rswfn/main/install.sh
curl -fsSL "$base" | sh -s -- --version 0.1.3
curl -fsSL "$base" | sh -s -- --dir "$HOME/bin"
curl -fsSL "$base" | sh -s -- --target x86_64-unknown-linux-gnu
curl -fsSL "$base" | sh -s -- --uninstall
curl -fsSL "$base" | sh -s -- --version 0.1.3 --base-url /mnt/usb/rswfn-mirror
```

```powershell
# Windows: download the script first, then run it with options
irm https://raw.githubusercontent.com/qmlCoder/rswfn/main/install.ps1 -OutFile install.ps1
powershell -ExecutionPolicy Bypass -File .\install.ps1 -Version 0.1.3
powershell -ExecutionPolicy Bypass -File .\install.ps1 -Dir "D:\tools\rswfn" -NoPath
powershell -ExecutionPolicy Bypass -File .\install.ps1 -Uninstall
powershell -ExecutionPolicy Bypass -File .\install.ps1 -Version 0.1.3 -BaseUrl D:\rswfn-mirror
```

The environment variables `RSWFN_VERSION`, `RSWFN_INSTALL_DIR`,
`RSWFN_BASE_URL`, `RSWFN_TARGET` and `RSWFN_SKIP_VERIFY` do the same job; the
command line takes precedence. See
[docs/install.md](https://github.com/qmlCoder/rswfn/blob/main/docs/install.md)
for the complete reference.

On Linux the installer never edits your shell configuration: it installs to
`~/.local/bin` and prints the `export PATH=...` line to add if that directory is
missing from your `PATH`. On Windows it appends the install directory to your
user `PATH` by default (pass `-NoPath` to skip it); reopen the terminal
afterwards.

If you prefer to inspect the script before running it:

```sh
curl -fsSL https://raw.githubusercontent.com/qmlCoder/rswfn/main/install.sh -o install.sh
# review install.sh, then run it
sh install.sh
```

### With Cargo (Rust toolchain required)

The same command works on every platform Cargo supports (Windows, Linux,
macOS):

```sh
cargo install rswfn
```

This builds `rswfn` from source; crates.io currently publishes version 0.1.3.

### As a library

```sh
cargo add rswfn
```

## Command-line usage

Running `rswfn` with no arguments launches the interactive TUI. With arguments
it runs as a CLI:

```sh
# Mulliken / Löwdin / Hirshfeld charges and spin population
rswfn atomprop charge mulliken C6H6.fch
rswfn atomprop charge lowdin C6H6.fch
rswfn atomprop charge hirshfeld C6H6.fch
rswfn atomprop charge spin C6H6.fch

# π-electron population (POCV / MOCV)
rswfn atomprop pi-pop C6H6.fch --method pocv

# Fukui functions and dual descriptor (needs the N, N+1 and N-1 wavefunctions)
rswfn atomprop fukui phenol.out phenol_N+1.out phenol_N-1.out

# Mayer / Wiberg bond orders
rswfn bondprop order mayer C6H6.fch

# σ/π/δ/φ decomposition of the bond order between atoms 1 and 2
rswfn bondprop decom 1 2 C6H6.fch

# Bond energy decomposition (bond / rest / cross) between atoms 1 and 2
# (closed-shell hybrid functional; reads the functional and basis from the log)
rswfn bondprop energy 1 2 C6H6.out

# Aromaticity (PISD)
rswfn moleprop aromacity C6H6.fch

# π electron energy decomposition (POCV, or --method mocv)
rswfn moleprop pi-energy C6H6.out

# Composition of molecular orbital 21
rswfn orbtprop composition 21 C6H6.fch

# Export a real-space function on a grid to a cube file
rswfn gridprop density C6H6.fch rho.cub --func density --step 0.2
# Electrostatic potential in a.u. (analytic Gaussian integrals, Multiwfn convention)
rswfn gridprop potential C6H6.fch esp.cub

# Closed-shell Hartree-Fock from a file that carries geometry + basis set
rswfn qmcal hf -i C6H6.fch -o C6H6_hf.fch

# Kohn-Sham DFT on the same kind of input (B3LYP by default)
rswfn qmcal dft -i C6H6.fch -o C6H6_b3lyp.fch

# Pick another exchange-correlation functional
rswfn qmcal dft -i C6H6.fch --func pbe0 -o C6H6_pbe0.fch

# Extract coordinates, convert formats, extract supporting information
rswfn coord C6H6.fch coords.xyz
rswfn trans C6H6.fch C6H6.gjf
rswfn si scan1.log scan2.log -o si.txt

# Preview the molecular structure in the terminal
# (colour comes from CPK colours; `--style space` and `--style wire` also work)
rswfn view C6H6.fch --size 100x40 --style ball --rotate 20,12
rswfn view phenol.fch --style space --color always
```

Run `rswfn --help`, or `rswfn <subcommand> --help`, for the full set of options.

## Library usage

```rust
use rswfn::props::atomprop::charge::Calculator;
use rswfn::base::Mole;

fn main() {
    // Read a wavefunction from a file; the format is detected from the extension.
    let mole = Mole::from_file("C6H6.fch");

    // Mulliken charges
    let charges = Calculator::new(&mole).mulliken();
    println!("{charges:?}");
}
```

## Feature flags

- `use-rsint` — use the pure-Rust `rsint` crate (a transliteration of libcint)
  for integrals. Enabled by default and builds on stable Rust.
- `use-rsbase` — resolve basis sets through the [`rsbse`](https://crates.io/crates/rsbse)
  Basis Set Exchange database instead of the embedded `src/rsbse/txts` files
  (`cargo build --features use-rsbase`). The data is not bundled with the crate:
  install it with `bse fetch-data` (from `cargo install rsbse`), or point
  `BSE_DATA_DIR` at an existing data directory.
- `use-rsxcs` — use the external pure-Rust `rsxcs` 0.2 crate (a Rust port of a
  subset of libxc 7.1.2) for exchange–correlation functionals. Enabled by
  default; disable it (for example with `--no-default-features`) to use the
  built-in DFTxclib translation in `src/rsxcs/DFTxclib.rs` instead.

All three backends are selected independently. With `use-rsint`, `use-rsbase`
and `use-rsxcs` disabled, `rswfn` builds entirely from its in-repo
implementations:

```sh
cargo build --no-default-features
```

## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT license ([LICENSE-MIT](LICENSE-MIT))

at your option.

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in this crate by you, as defined in the Apache-2.0 license, shall
be dual licensed as above, without any additional terms or conditions.

## 中文简介

`rswfn` 是一个用 Rust 编写的量子化学**波函数分析库与命令行工具**，可以读取
Gaussian、ORCA、Molden 等程序的输出文件，进行各类波函数分析，包括：原子电荷
布居（Mulliken / Löwdin / Hirshfeld）、π 电子布居（POCV / MOCV）、福井函数与
双描述符、自由价、键级（Mayer / Wiberg）及其分解、芳香性、轨道成分分析，以及
实空间函数（密度、静电势、ELF、RDG 等）的格点计算与导出，以及基于坐标与基
组的闭壳层 Hartree-Fock（RHF）与 Kohn-Sham DFT 计算（B3LYP、PBE0、PBE、
BLYP、BP86、SVWN5 可选）。另外可以在终端里直接预览分子结构：`rswfn view 文件`
或在 TUI 里按 `F4`，支持球棍 / 空间填充 / 线框三种模型、方向键旋转、真彩色 CPK
配色，不需要打开任何图形程序。既可以作为库使用，也可以直接作为命令行 / TUI 程序
运行。

### 安装

- **Windows**（PowerShell）：`irm https://raw.githubusercontent.com/qmlCoder/rswfn/main/install.ps1 | iex`
- **Linux**：`curl -fsSL https://raw.githubusercontent.com/qmlCoder/rswfn/main/install.sh | sh`
- **已安装 Rust 工具链**（Windows / Linux / macOS 通用）：`cargo install rswfn`

安装脚本从 GitHub Releases 下载预编译二进制、校验 SHA256 后安装，不需要 Rust
工具链和管理员权限。Linux x86_64 默认安装**全静态 musl** 产物，在任意发行版或容器
（Alpine、CentOS 7、Debian 11 等）中都能直接运行，不依赖 glibc 版本；aarch64 目前只
提供 gnu 产物，要求 glibc ≥ 2.39；高级用户可用 `--target x86_64-unknown-linux-gnu`
安装动态链接的 gnu 版本（性能略好，但同样要求 glibc ≥ 2.39）。Windows 产物使用静态
CRT，是单个 `.exe`，不依赖 `vcruntime140.dll`。

注意：最新版本 0.1.3 已发布到 crates.io，预编译二进制来自 GitHub Releases（首个
release 为 `v0.1.1`）。Windows 目前只提供 x86_64 产物（ARM64 暂不支持），默认安装到
`%LOCALAPPDATA%\Programs\rswfn` 并写入用户 PATH，需重开终端生效（`-NoPath` 可跳过）。
Linux 默认安装到 `~/.local/bin`，脚本不会修改 shell 配置文件，需要自行把该目录加入
PATH。指定版本、自定义目录、离线安装、卸载等用法详见
[docs/install.md](docs/install.md)。
