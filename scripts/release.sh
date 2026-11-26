#!/bin/sh
# rswfn 本地发布脚本：交叉编译 → 打包 → 分别上传阿里云 OSS 与 GitHub Release。
#
# 在 Windows 的 Git Bash 里运行：Windows 目标本机编译，3 个 Linux 目标交给 WSL；
# 在 Linux 上运行：Linux 目标本机编译（Windows 目标会被跳过）。
#
# 一次性准备：
#   Windows  rustup target add x86_64-pc-windows-msvc
#   WSL      rustup target add aarch64-unknown-linux-gnu
#            sudo apt install musl-tools gcc-aarch64-linux-gnu
#   Python   python -m pip install oss2          # 只有上传 OSS 时才需要
#   凭据     cp .env.example .env 并填入 OSS AccessKey（见 .env.example）
#   gh       gh auth status                      # 需要仓库写权限
#
# 用法：
#   sh scripts/release.sh                        # 编译 + 打包 + 双上传
#   sh scripts/release.sh --skip-build           # 跳过编译，只重新打包 + 上传
#   sh scripts/release.sh --no-upload            # 只编译 + 打包
#   sh scripts/release.sh --no-oss               # 不传 OSS
#   sh scripts/release.sh --targets "x86_64-pc-windows-msvc"
#   sh scripts/release.sh --version 0.1.4
#   sh scripts/release.sh --wsl Ubuntu-22.04 --wsl-target-dir /home/me/rswfn-target
#
# 产出的资产与 .github/workflows/release.yml 完全一致：
#   rswfn-<ver>-x86_64-pc-windows-msvc.zip          （扁平结构）
#   rswfn-<ver>-x86_64-unknown-linux-gnu.tar.gz     （含顶层目录）
#   rswfn-<ver>-x86_64-unknown-linux-musl.tar.gz
#   rswfn-<ver>-aarch64-unknown-linux-gnu.tar.gz
#   SHA256SUMS

set -eu

cd "$(dirname "$0")/.."

log() { printf '==> %s\n' "$*"; }
warn() { printf '==> WARNING: %s\n' "$*" >&2; }
die() {
    printf '==> ERROR: %s\n' "$*" >&2
    exit 1
}

usage() {
    sed -n '2,30p' "$0" | sed 's/^# \{0,1\}//'
}

VER=
TARGETS='x86_64-pc-windows-msvc x86_64-unknown-linux-gnu x86_64-unknown-linux-musl aarch64-unknown-linux-gnu'
DO_BUILD=1
DO_UPLOAD=1
DO_OSS=1
DO_GH=1
DO_TAG=1
DIST=dist
WSL_DISTRO=Ubuntu-22.04
# WSL 里编译到 Linux 原生文件系统：在 /mnt/e（DrvFs）上编译会慢很多
WSL_TARGET_DIR='~/.cache/rswfn-target'

while [ $# -gt 0 ]; do
    case "$1" in
        --version)
            [ $# -ge 2 ] || die '--version 需要参数'
            VER=$2
            shift 2
            ;;
        --targets)
            [ $# -ge 2 ] || die '--targets 需要参数'
            TARGETS=$2
            shift 2
            ;;
        --dist)
            [ $# -ge 2 ] || die '--dist 需要参数'
            DIST=$2
            shift 2
            ;;
        --wsl)
            [ $# -ge 2 ] || die '--wsl 需要参数'
            WSL_DISTRO=$2
            shift 2
            ;;
        --wsl-target-dir)
            [ $# -ge 2 ] || die '--wsl-target-dir 需要参数'
            WSL_TARGET_DIR=$2
            shift 2
            ;;
        --skip-build)
            DO_BUILD=0
            shift
            ;;
        --no-upload)
            DO_UPLOAD=0
            shift
            ;;
        --no-oss)
            DO_OSS=0
            shift
            ;;
        --no-github)
            DO_GH=0
            shift
            ;;
        --no-tag)
            DO_TAG=0
            shift
            ;;
        -h | --help)
            usage
            exit 0
            ;;
        *) die "未知参数: $1（用 --help 看用法）" ;;
    esac
done

# --- 版本与 tag -------------------------------------------------------------

command -v python >/dev/null 2>&1 || die '找不到 python（打包 zip 需要）'

CARGO_VER=$(sed -n 's/^version *= *"\(.*\)"/\1/p' Cargo.toml | head -n 1)
[ -n "$CARGO_VER" ] || die '无法从 Cargo.toml 读出版本号'
[ -n "${VER:-}" ] || VER=$CARGO_VER
[ "$VER" = "$CARGO_VER" ] || die "版本不一致：--version $VER，Cargo.toml $CARGO_VER"
TAG=v$VER
log "版本 $VER（tag $TAG）"

if [ -n "$(git status --porcelain --untracked-files=no 2>/dev/null)" ]; then
    warn '工作区有未提交的改动，产物将对应工作区内容而不是某个提交'
fi

if [ "$DO_TAG" = 1 ]; then
    if git rev-parse -q --verify "refs/tags/$TAG" >/dev/null 2>&1; then
        log "tag $TAG 已存在"
    else
        log "创建 tag $TAG"
        git tag "$TAG"
    fi
    log "推送 tag $TAG 到 origin"
    git push origin "$TAG"
fi

# --- 目标分组 ---------------------------------------------------------------

HOST_OS=linux
case "$(uname -s)" in
    MINGW* | MSYS* | CYGWIN*) HOST_OS=windows ;;
esac

WIN_TARGETS=
LINUX_TARGETS=
for t in $TARGETS; do
    case "$t" in
        *windows*) WIN_TARGETS="$WIN_TARGETS $t" ;;
        *) LINUX_TARGETS="$LINUX_TARGETS $t" ;;
    esac
done

# --- 编译 -------------------------------------------------------------------

if [ "$DO_BUILD" = 1 ]; then
    for t in $WIN_TARGETS; do
        [ "$HOST_OS" = windows ] || die "在非 Windows 主机上无法编译 $t"
        log "编译 $t（本机，静态 CRT）"
        # 与 release.yml 一致：静态链接 CRT，产物不依赖 vcruntime140.dll
        RUSTFLAGS='-C target-feature=+crt-static' cargo build --release --locked --target "$t"
    done

    if [ -n "$LINUX_TARGETS" ]; then
        if [ "$HOST_OS" = windows ]; then
            log "Linux 目标交给 WSL（$WSL_DISTRO）：$LINUX_TARGETS"
            # 不用 wslpath：Git Bash 把 E:\code\rswfn 当参数传给原生程序时反斜杠会被
            # MSYS 吃掉，得到 E:coderswfn。Git Bash 的 /e/code/rswfn 对应 WSL 的
            # /mnt/e/code/rswfn，直接前接 /mnt 即可。
            case "$PWD" in
                /?/*) REPO_WSL=/mnt$PWD ;;
                *) die "无法把 $PWD 映射成 WSL 路径（期望形如 /e/code/rswfn）" ;;
            esac
            # 整条命令作为一个参数传给 bash -lc，避免 MSYS 再改写里面的路径。
            # 先去掉可能的 CRLF：脚本在 Windows 上被编辑器改成 CRLF 后，WSL 的 dash 会在
            # `set -eu` 上报 "Illegal option -"。
            # $LINUX_TARGETS 自带前导空格，直接拼在版本号后面
            wsl -d "$WSL_DISTRO" -- bash -lc "cd '$REPO_WSL' || exit 1; sed -i 's/\r//g' scripts/*.sh 2>/dev/null || true; CARGO_TARGET_DIR=$WSL_TARGET_DIR RSWFN_DIST='$DIST' sh scripts/build-wsl.sh $VER$LINUX_TARGETS"
        else
            RSWFN_DIST="$DIST" sh scripts/build-wsl.sh "$VER" $LINUX_TARGETS
        fi
    fi
else
    log '跳过编译，复用已有产物'
fi

# --- 打包（Windows zip） ----------------------------------------------------
#
# Linux 的 tar.gz 由 build-wsl.sh 在 WSL 里打包（tar 元数据更干净）。

mkdir -p "$DIST"

for t in $WIN_TARGETS; do
    [ "$HOST_OS" = windows ] || continue
    bin="target/$t/release/rswfn.exe"
    [ -f "$bin" ] || die "找不到 $bin：先去掉 --skip-build 编译，或检查 target 名"

    stage="$DIST/rswfn-$VER-$t"
    asset="rswfn-$VER-$t.zip"
    rm -rf "$stage"
    mkdir -p "$stage"
    cp "$bin" "$stage/"
    cp README.md LICENSE-MIT LICENSE-APACHE "$stage/"

    log "打包 $asset"
    rm -f "$DIST/$asset"
    # CI 用 7z 从 stage 目录内部打包，zip 是扁平结构；这里用 zipfile 保持同样结构，
    # 同时避开 Compress-Archive 在非交互控制台下 Write-Progress 报错的问题。
    python -c 'import os, sys, zipfile
stage, out = sys.argv[1], sys.argv[2]
with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as zf:
    for name in sorted(os.listdir(stage)):
        zf.write(os.path.join(stage, name), name)' "$stage" "$DIST/$asset" ||
        die "打包 $asset 失败"
    rm -rf "$stage"
done

# --- 校验和 -----------------------------------------------------------------

log '生成 SHA256SUMS'
# 不用 sha256sum：Windows 上的 coreutils 会输出 `hash *name`（星号 + 单空格），
# 而契约要求 `<64hex><两个空格><name>`，且要保证哈希的是真实字节。
python -c 'import hashlib, pathlib, sys
dist = pathlib.Path(sys.argv[1])
rows = []
for path in sorted(dist.iterdir()):
    if path.is_file() and path.suffix in (".zip", ".gz"):
        rows.append((path.name, hashlib.sha256(path.read_bytes()).hexdigest()))
if not rows:
    sys.exit(f"{dist} 里没有任何归档，没什么可发布的")
with (dist / "SHA256SUMS").open("w", encoding="utf-8", newline="\n") as fh:
    for name, digest in rows:
        fh.write(f"{digest}  {name}\n")' "$DIST" || die '生成 SHA256SUMS 失败'
cat "$DIST/SHA256SUMS"

# --- 上传 -------------------------------------------------------------------

if [ "$DO_UPLOAD" = 1 ]; then
    if [ "$DO_OSS" = 1 ]; then
        log '上传到阿里云 OSS'
        python -c 'import oss2' 2>/dev/null || die '缺少 oss2：python -m pip install oss2'
        python scripts/upload_oss.py "$VER" "$DIST" --installers --latest
    fi

    if [ "$DO_GH" = 1 ]; then
        log '上传到 GitHub Release'
        REPO=$(gh repo view --json nameWithOwner --jq .nameWithOwner)
        assets=
        for f in "$DIST"/*.zip "$DIST"/*.tar.gz "$DIST"/SHA256SUMS; do
            [ -f "$f" ] && assets="$assets $f"
        done
        [ -n "$assets" ] || die "没有可上传的资产"

        if gh release view "$TAG" --repo "$REPO" >/dev/null 2>&1; then
            log "Release $TAG 已存在，覆盖同名资产"
            # shellcheck disable=SC2086
            gh release upload "$TAG" $assets --clobber --repo "$REPO"
        else
            log "创建 Release $TAG"
            # shellcheck disable=SC2086
            gh release create "$TAG" $assets --repo "$REPO" --title "$TAG" --generate-notes
        fi
        log "https://github.com/$REPO/releases/tag/$TAG"
    fi
fi

log '完成'
