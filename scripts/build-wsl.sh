#!/bin/sh
# 在 WSL/Linux 里编译 Linux 目标并打包成 tar.gz。
# 由 scripts/release.sh 调用，不建议直接跑。
#
# 用法：sh scripts/build-wsl.sh <version> <target>...
#
# 环境变量：
#   CARGO_TARGET_DIR  可选；指向 WSL 原生文件系统（如 ~/.cache/rswfn-target）
#                     可以明显加快构建，避免在 /mnt/e 上编译。

set -eu

[ $# -ge 2 ] || {
    printf '用法: %s <version> <target>...\n' "$0" >&2
    exit 2
}

cd "$(dirname "$0")/.."

VER=$1
shift
DIST=${RSWFN_DIST:-dist}

log() { printf '==> %s\n' "$*"; }
die() {
    printf '==> ERROR: %s\n' "$*" >&2
    exit 1
}

mkdir -p "$DIST"

for t in "$@"; do
    # 交叉编译需要显式指定链接器和 C 编译器（vendored OpenSSL 也要用它）
    case "$t" in
        x86_64-unknown-linux-gnu) ;;
        x86_64-unknown-linux-musl)
            export CARGO_TARGET_X86_64_UNKNOWN_LINUX_MUSL_LINKER=musl-gcc
            export CC_x86_64_unknown_linux_musl=musl-gcc
            ;;
        aarch64-unknown-linux-gnu)
            export CARGO_TARGET_AARCH64_UNKNOWN_LINUX_GNU_LINKER=aarch64-linux-gnu-gcc
            export CC_aarch64_unknown_linux_gnu=aarch64-linux-gnu-gcc
            ;;
        *)
            die "未知目标 $t（本脚本只处理 Linux 目标）"
            ;;
    esac

    command -v cargo >/dev/null 2>&1 || die 'WSL 里找不到 cargo，先装 Rust 工具链'
    rustup target list --installed 2>/dev/null | grep -qx "$t" || {
        die "WSL 缺少 target $t：rustup target add $t"
    }

    log "编译 $t"
    cargo build --release --locked --target "$t"

    target_dir=${CARGO_TARGET_DIR:-target}
    bin="$target_dir/$t/release/rswfn"
    [ -f "$bin" ] || die "找不到 $bin"

    stage="$DIST/rswfn-$VER-$t"
    asset="rswfn-$VER-$t.tar.gz"
    rm -rf "$stage"
    mkdir -p "$stage"
    cp "$bin" "$stage/"
    cp README.md LICENSE-MIT LICENSE-APACHE "$stage/"
    chmod 755 "$stage/rswfn"

    log "打包 $asset"
    rm -f "$DIST/$asset"
    # 与 release.yml 一致：tar 里带顶层目录 rswfn-<ver>-<target>/
    tar -czf "$DIST/$asset" -C "$DIST" "rswfn-$VER-$t"
    rm -rf "$stage"
done

log 'Linux 目标构建完成'
