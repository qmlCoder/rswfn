#!/bin/sh
# rswfn installer for Linux.
#
# Installs the prebuilt rswfn binary from GitHub Releases. It needs no Rust
# toolchain, no root privileges and never calls sudo.
#
# Typical use:
#   curl -fsSL https://raw.githubusercontent.com/qmlCoder/rswfn/main/install.sh | sh
#
# This script is POSIX sh (dash / busybox sh compatible). It never reads from
# stdin because stdin is the script itself in the `curl ... | sh` form.
# It never modifies shell rc files: it only prints the line to add to PATH.

set -u

DEFAULT_BASE_URL='https://github.com/qmlCoder/rswfn/releases/download'
LATEST_PAGE_URL='https://github.com/qmlCoder/rswfn/releases/latest'
LATEST_API_URL='https://api.github.com/repos/qmlCoder/rswfn/releases/latest'

VERSION=''
INSTALL_DIR=''
BASE_URL=''
TARGET=''
SKIP_VERIFY=0
UNINSTALL=0

BASE_MODE=''
BASE_PATH=''
WORK_DIR=''
ARCHIVE=''
SRC_BINARY=''

# ---------------------------------------------------------------------------
# Output helpers. Progress lines always start with "==>".
# ---------------------------------------------------------------------------
say() {
    printf '==> %s\n' "$*"
}

warn() {
    printf '==> WARNING: %s\n' "$*" >&2
}

error_line() {
    printf '==> ERROR: %s\n' "$*" >&2
}

die() {
    error_line "$1"
    exit 1
}

usage_error() {
    error_line "$1"
    usage >&2
    exit 2
}

usage() {
    cat <<'USAGE_EOF'
Usage: install.sh [options]

Install the rswfn command-line tool from GitHub Releases.

Options:
  -h, --help               Show this help and exit.
      --version <version>  Version to install, e.g. 0.1.0 or v0.1.0.
      --dir <path>         Install directory (default: $HOME/.local/bin).
      --base-url <url>     Download base URL or local release directory.
      --skip-verify        Do not verify the SHA256 checksum.
      --target <triple>    Override the detected target triple.
      --uninstall          Remove <dir>/rswfn and exit.

Environment variables (command-line options take precedence):
  RSWFN_VERSION, RSWFN_INSTALL_DIR, RSWFN_BASE_URL,
  RSWFN_SKIP_VERIFY, RSWFN_TARGET

Examples:
  curl -fsSL https://raw.githubusercontent.com/qmlCoder/rswfn/main/install.sh | sh
  sh install.sh --version 0.1.0 --target x86_64-unknown-linux-gnu
  sh install.sh --uninstall
USAGE_EOF
}

# ---------------------------------------------------------------------------
# Cleanup: the temporary directory is removed on exit, success or failure.
# ---------------------------------------------------------------------------
cleanup() {
    if [ -n "${WORK_DIR}" ] && [ -d "${WORK_DIR}" ]; then
        rm -rf "${WORK_DIR}" 2>/dev/null || true
        WORK_DIR=''
    fi
}

trap cleanup 0
trap 'cleanup; exit 1' 1 2 3 15

# ---------------------------------------------------------------------------
# Environment defaults (command line options are parsed afterwards).
# ---------------------------------------------------------------------------
if [ -n "${RSWFN_VERSION:-}" ]; then
    VERSION="${RSWFN_VERSION}"
fi
if [ -n "${RSWFN_INSTALL_DIR:-}" ]; then
    INSTALL_DIR="${RSWFN_INSTALL_DIR}"
fi
if [ -n "${RSWFN_BASE_URL:-}" ]; then
    BASE_URL="${RSWFN_BASE_URL}"
fi
if [ -n "${RSWFN_TARGET:-}" ]; then
    TARGET="${RSWFN_TARGET}"
fi
if [ -n "${RSWFN_SKIP_VERIFY:-}" ]; then
    case "${RSWFN_SKIP_VERIFY}" in
        1|true|TRUE|True|yes|YES|on|ON) SKIP_VERIFY=1 ;;
    esac
fi

# ---------------------------------------------------------------------------
# Command line parsing.
# ---------------------------------------------------------------------------
while [ "$#" -gt 0 ]; do
    case "$1" in
        -h|--help)
            usage
            exit 0
            ;;
        --version)
            [ "$#" -ge 2 ] || usage_error "option '--version' requires a value, e.g. --version 0.1.0"
            VERSION="$2"
            shift 2
            ;;
        --version=*)
            VERSION="${1#*=}"
            shift
            ;;
        --dir)
            [ "$#" -ge 2 ] || usage_error "option '--dir' requires a value, e.g. --dir \"\$HOME/.local/bin\""
            INSTALL_DIR="$2"
            shift 2
            ;;
        --dir=*)
            INSTALL_DIR="${1#*=}"
            shift
            ;;
        --base-url)
            [ "$#" -ge 2 ] || usage_error "option '--base-url' requires a value"
            BASE_URL="$2"
            shift 2
            ;;
        --base-url=*)
            BASE_URL="${1#*=}"
            shift
            ;;
        --target)
            [ "$#" -ge 2 ] || usage_error "option '--target' requires a value, e.g. --target x86_64-unknown-linux-gnu"
            TARGET="$2"
            shift 2
            ;;
        --target=*)
            TARGET="${1#*=}"
            shift
            ;;
        --skip-verify)
            SKIP_VERIFY=1
            shift
            ;;
        --uninstall)
            UNINSTALL=1
            shift
            ;;
        -*)
            usage_error "unknown option '$1'"
            ;;
        *)
            usage_error "unexpected argument '$1'"
            ;;
    esac
done

# ---------------------------------------------------------------------------
# Install directory.
# ---------------------------------------------------------------------------
prepare_install_dir() {
    if [ -n "${INSTALL_DIR}" ]; then
        mkdir -p "${INSTALL_DIR}" 2>/dev/null || die "cannot create install directory '${INSTALL_DIR}': check the path or pass --dir <path>"
        [ -d "${INSTALL_DIR}" ] || die "install directory '${INSTALL_DIR}' does not exist: pass --dir <path>"
        [ -w "${INSTALL_DIR}" ] || die "install directory '${INSTALL_DIR}' is not writable: pass --dir <writable path>"
        return 0
    fi
    home_value="${HOME:-}"
    if [ -n "${home_value}" ]; then
        candidate="${home_value}/.local/bin"
        if mkdir -p "${candidate}" 2>/dev/null && [ -w "${candidate}" ]; then
            INSTALL_DIR="${candidate}"
            return 0
        fi
    fi
    if [ -d '/usr/local/bin' ] && [ -w '/usr/local/bin' ]; then
        INSTALL_DIR='/usr/local/bin'
        return 0
    fi
    die "cannot find a writable install directory (tried \$HOME/.local/bin and /usr/local/bin); pass --dir <path>, e.g. --dir \"\$HOME/.local/bin\""
}

# ---------------------------------------------------------------------------
# Resolve the install directory for uninstall: never create directories.
# ---------------------------------------------------------------------------
resolve_install_dir_without_create() {
    if [ -n "${INSTALL_DIR}" ]; then
        return 0
    fi
    home_value="${HOME:-}"
    if [ -n "${home_value}" ]; then
        INSTALL_DIR="${home_value}/.local/bin"
        return 0
    fi
    INSTALL_DIR='/usr/local/bin'
}

# ---------------------------------------------------------------------------
# Uninstall: remove only <dir>/rswfn. No rc files, no sudo, no directories.
# ---------------------------------------------------------------------------
if [ "${UNINSTALL}" -eq 1 ]; then
    resolve_install_dir_without_create
    uninstall_bin="${INSTALL_DIR}/rswfn"
    if [ -d "${INSTALL_DIR}" ] && { [ -e "${uninstall_bin}" ] || [ -L "${uninstall_bin}" ]; }; then
        rm -f "${uninstall_bin}" || die "failed to remove '${uninstall_bin}': check permissions or run with --dir <path>"
        say "removed ${uninstall_bin}"
    else
        say "rswfn is not installed in ${INSTALL_DIR}; nothing to remove"
    fi
    exit 0
fi

# ---------------------------------------------------------------------------
# Version parsing.
# ---------------------------------------------------------------------------
normalize_version() {
    case "${VERSION}" in
        v*|V*) VERSION="${VERSION#?}" ;;
    esac
    [ -n "${VERSION}" ] || usage_error "version must not be empty: pass --version x.y.z"
    case "${VERSION}" in
        *[!0-9A-Za-z._-]*)
            usage_error "invalid version '${VERSION}': only characters [0-9A-Za-z._-] are allowed"
            ;;
    esac
}

# ---------------------------------------------------------------------------
# Download helpers: remote URL or local release directory.
# ---------------------------------------------------------------------------
normalize_base_url() {
    if [ -z "${BASE_URL}" ]; then
        BASE_URL="${DEFAULT_BASE_URL}"
    fi
    case "${BASE_URL}" in
        http://*|https://*)
            BASE_MODE='remote'
            while [ "${BASE_URL%/}" != "${BASE_URL}" ]; do
                BASE_URL="${BASE_URL%/}"
            done
            [ -n "${BASE_URL}" ] || usage_error "invalid base URL: '${BASE_URL}'"
            ;;
        file://*)
            BASE_MODE='local'
            BASE_PATH="${BASE_URL#file://}"
            case "${BASE_PATH}" in
                /[A-Za-z]:/*) BASE_PATH="${BASE_PATH#/}" ;;
            esac
            ;;
        *)
            BASE_MODE='local'
            BASE_PATH="${BASE_URL}"
            ;;
    esac
}

fetch_url() {
    url="$1"
    dest="$2"
    if command -v curl >/dev/null 2>&1; then
        curl -fsSL -o "${dest}" "${url}"
        return $?
    fi
    if command -v wget >/dev/null 2>&1; then
        wget -q -O "${dest}" "${url}"
        return $?
    fi
    error_line "neither 'curl' nor 'wget' is available: install curl (or wget) and retry"
    return 1
}

download_asset() {
    # $1 = asset file name, $2 = destination path
    asset_name="$1"
    asset_dest="$2"
    if [ "${BASE_MODE}" = 'remote' ]; then
        fetch_url "${BASE_URL}/v${VERSION}/${asset_name}" "${asset_dest}"
        return $?
    fi
    asset_src="${BASE_PATH}/v${VERSION}/${asset_name}"
    if [ ! -f "${asset_src}" ]; then
        return 1
    fi
    cp "${asset_src}" "${asset_dest}" 2>/dev/null
}

# ---------------------------------------------------------------------------
# Target detection. An explicit --target / RSWFN_TARGET skips uname probing.
# ---------------------------------------------------------------------------
detect_target() {
    if [ -n "${TARGET}" ]; then
        case "${TARGET}" in
            *[!0-9A-Za-z._-]*)
                usage_error "invalid target '${TARGET}': only characters [0-9A-Za-z._-] are allowed"
                ;;
        esac
        return 0
    fi
    os_name=$(uname -s 2>/dev/null) || os_name=''
    case "${os_name}" in
        Linux*) : ;;
        Darwin*)
            die "macOS is not supported by this installer: install rswfn with 'cargo install rswfn' instead"
            ;;
        '')
            die "could not detect the operating system with 'uname -s': pass --target <triple>"
            ;;
        *)
            die "unsupported operating system '${os_name}': this installer supports Linux only; pass --target <triple> to override"
            ;;
    esac
    cpu_arch=$(uname -m 2>/dev/null) || cpu_arch=''
    case "${cpu_arch}" in
        x86_64|amd64)
            TARGET='x86_64-unknown-linux-musl'
            ;;
        aarch64|arm64)
            TARGET='aarch64-unknown-linux-gnu'
            ;;
        '')
            die "could not detect the CPU architecture with 'uname -m': pass --target <triple>"
            ;;
        *)
            die "unsupported CPU architecture '${cpu_arch}': pass --target <triple>, e.g. --target x86_64-unknown-linux-musl"
            ;;
    esac
}

# ---------------------------------------------------------------------------
# Latest version from the GitHub releases API (only when not specified).
# ---------------------------------------------------------------------------
resolve_latest_from_redirect() {
    # Resolve the latest tag from the releases/latest 301/302 redirect: this does
    # not consume the GitHub API quota (the anonymous API allows only 60 requests
    # per IP per hour, which CI runners and corporate NATs hit easily with 403).
    latest_url=''
    if command -v curl >/dev/null 2>&1; then
        latest_url=$(curl -fsSL -o /dev/null -w '%{url_effective}' "${LATEST_PAGE_URL}" 2>/dev/null) || latest_url=''
    fi
    if [ -z "${latest_url}" ] && command -v wget >/dev/null 2>&1; then
        latest_url=$(wget -q -S --max-redirect=0 -O /dev/null "${LATEST_PAGE_URL}" 2>&1 \
            | sed -n 's/^[[:space:]]*[Ll]ocation:[[:space:]]*//p' | tr -d '\r' | head -n 1) || latest_url=''
    fi
    case "${latest_url}" in
        */releases/tag/*)
            tag_from_redirect="${latest_url##*/releases/tag/}"
            tag_from_redirect="${tag_from_redirect%%[?#]*}"
            printf '%s' "${tag_from_redirect}"
            ;;
        *)
            : ;;
    esac
}

resolve_version() {
    if [ -n "${VERSION}" ]; then
        normalize_version
        return 0
    fi
    if [ "${BASE_MODE}" = 'local' ]; then
        die "no version specified: pass --version x.y.z when using a local base directory '${BASE_PATH}' (layout: <base>/v<version>/<asset>)"
    fi
    say 'resolving the latest rswfn release from GitHub'
    tag_name=$(resolve_latest_from_redirect)
    if [ -z "${tag_name}" ]; then
        warn 'could not resolve the latest tag from the releases page redirect; falling back to the GitHub API'
        latest_json="${WORK_DIR}/latest-release.json"
        if ! fetch_url "${LATEST_API_URL}" "${latest_json}"; then
            die "could not determine the latest release (network unavailable, or the anonymous GitHub API rate limit was hit): pass --version x.y.z explicitly (for example: sh install.sh --version 0.1.1)"
        fi
        tag_name=$(sed -n 's/.*"tag_name"[[:space:]]*:[[:space:]]*"\([^"]*\)".*/\1/p' "${latest_json}" | head -n 1)
        if [ -z "${tag_name}" ]; then
            die "could not parse 'tag_name' from the GitHub API response: pass --version x.y.z explicitly"
        fi
    fi
    VERSION="${tag_name}"
    normalize_version
    say "using release v${VERSION}"
}

# ---------------------------------------------------------------------------
# Temporary directory.
# ---------------------------------------------------------------------------
create_work_dir() {
    WORK_DIR=$(mktemp -d 2>/dev/null) || WORK_DIR=''
    if [ -z "${WORK_DIR}" ] || [ ! -d "${WORK_DIR}" ]; then
        die "cannot create a temporary directory with 'mktemp -d': check that mktemp is installed and \$TMPDIR is writable"
    fi
}

# ---------------------------------------------------------------------------
# SHA256 checksum verification.
# ---------------------------------------------------------------------------
verify_archive() {
    if [ "${SKIP_VERIFY}" -eq 1 ]; then
        say 'skipping checksum verification (--skip-verify)'
        return 0
    fi
    sums_file="${WORK_DIR}/SHA256SUMS"
    if ! download_asset 'SHA256SUMS' "${sums_file}"; then
        warn "SHA256SUMS is not available for v${VERSION}: skipping checksum verification (older release?)"
        return 0
    fi
    sums_lf="${WORK_DIR}/SHA256SUMS.lf"
    tr -d '\r' < "${sums_file}" > "${sums_lf}" 2>/dev/null || cp "${sums_file}" "${sums_lf}"
    # SHA256SUMS lines are "<hash>  <file>"; a few toolchains (e.g. MSYS/Cygwin
    # sha256sum in binary mode) write "<hash> *<file>". Accept both forms.
    expected_hash=$(awk -v name="${ASSET}" '{ n = $2; sub(/^\*/, "", n); if (n == name) { print $1; exit } }' "${sums_lf}" 2>/dev/null)
    if [ -z "${expected_hash}" ]; then
        die "SHA256SUMS does not contain an entry for '${ASSET}': refusing to install an unverified archive (use --skip-verify only if you trust the download)"
    fi
    if command -v sha256sum >/dev/null 2>&1; then
        actual_hash=$(sha256sum "${ARCHIVE}" | awk '{ print $1 }')
    elif command -v shasum >/dev/null 2>&1; then
        actual_hash=$(shasum -a 256 "${ARCHIVE}" | awk '{ print $1 }')
    else
        warn "neither 'sha256sum' nor 'shasum' is available: skipping checksum verification"
        return 0
    fi
    if [ -z "${actual_hash}" ]; then
        die "failed to compute the SHA256 checksum of '${ARCHIVE}'"
    fi
    expected_hash=$(printf '%s' "${expected_hash}" | tr 'A-Z' 'a-z')
    actual_hash=$(printf '%s' "${actual_hash}" | tr 'A-Z' 'a-z')
    if [ "${expected_hash}" != "${actual_hash}" ]; then
        error_line "checksum mismatch for '${ASSET}'"
        error_line "  expected: ${expected_hash}"
        error_line "  actual:   ${actual_hash}"
        exit 1
    fi
    say "verified SHA256 checksum of ${ASSET}"
}

# ---------------------------------------------------------------------------
# Extract the archive and locate the binary.
# ---------------------------------------------------------------------------
extract_archive() {
    command -v tar >/dev/null 2>&1 || die "'tar' is required to extract '${ASSET}' but was not found: install tar and retry"
    extract_dir="${WORK_DIR}/extract"
    mkdir -p "${extract_dir}" || die "cannot create '${extract_dir}'"
    if ! tar -xzf "${ARCHIVE}" -C "${extract_dir}" 2>"${WORK_DIR}/tar.err"; then
        if command -v gzip >/dev/null 2>&1; then
            if ! gzip -dc "${ARCHIVE}" 2>/dev/null | tar -xf - -C "${extract_dir}" 2>>"${WORK_DIR}/tar.err"; then
                die "failed to extract '${ASSET}': $(cat "${WORK_DIR}/tar.err" 2>/dev/null)"
            fi
        else
            die "failed to extract '${ASSET}': $(cat "${WORK_DIR}/tar.err" 2>/dev/null)"
        fi
    fi
}

locate_binary() {
    extract_dir="${WORK_DIR}/extract"
    if [ -f "${extract_dir}/rswfn" ]; then
        SRC_BINARY="${extract_dir}/rswfn"
        return 0
    fi
    found_dir=''
    dir_count=0
    for entry in "${extract_dir}"/*; do
        [ -e "${entry}" ] || continue
        if [ -d "${entry}" ]; then
            dir_count=$((dir_count + 1))
            found_dir="${entry}"
        fi
    done
    if [ "${dir_count}" -eq 1 ] && [ -f "${found_dir}/rswfn" ]; then
        SRC_BINARY="${found_dir}/rswfn"
        return 0
    fi
    die "could not find the 'rswfn' binary in '${ASSET}' (unexpected archive layout); please report this at https://github.com/qmlCoder/rswfn/issues"
}

# ---------------------------------------------------------------------------
# Install the binary atomically inside the install directory.
# ---------------------------------------------------------------------------
install_binary() {
    mkdir -p "${INSTALL_DIR}" || die "cannot create install directory '${INSTALL_DIR}': pass --dir <path>"
    target_bin="${INSTALL_DIR}/rswfn"
    staged_bin="${INSTALL_DIR}/.rswfn.new.$$"
    rm -f "${staged_bin}" 2>/dev/null || true
    cp "${SRC_BINARY}" "${staged_bin}" || die "failed to copy the binary into '${INSTALL_DIR}': check permissions or pass --dir <path>"
    chmod 755 "${staged_bin}" || {
        rm -f "${staged_bin}" 2>/dev/null || true
        die "failed to make '${staged_bin}' executable"
    }
    mv -f "${staged_bin}" "${target_bin}" || {
        rm -f "${staged_bin}" 2>/dev/null || true
        die "failed to install the binary as '${target_bin}': check permissions or pass --dir <path>"
    }
}

# ---------------------------------------------------------------------------
# Post-install self test. Failures only warn: older binaries do not know
# --version, and that must not fail the installation.
# ---------------------------------------------------------------------------
self_test() {
    installed_bin="${INSTALL_DIR}/rswfn"
    if "${installed_bin}" --version </dev/null >"${WORK_DIR}/version.out" 2>&1; then
        version_output=$(cat "${WORK_DIR}/version.out" 2>/dev/null)
        if [ -n "${version_output}" ]; then
            say "self-test ok: ${version_output}"
            return 0
        fi
    fi
    warn "the installed binary did not report a version for '--version' (normal for very old releases)"
    if "${installed_bin}" --help </dev/null >"${WORK_DIR}/help.out" 2>&1; then
        if [ -s "${WORK_DIR}/help.out" ]; then
            printf '==> installed binary help output:\n'
            cat "${WORK_DIR}/help.out"
        fi
    fi
    return 0
}

# ---------------------------------------------------------------------------
# PATH hint. install.sh never modifies shell rc files.
# ---------------------------------------------------------------------------
path_hint() {
    case ":${PATH:-}:" in
        *":${INSTALL_DIR}:"*) return 0 ;;
    esac
    printf '\n'
    say 'rswfn is not in your PATH yet. Add this line to ~/.bashrc (or ~/.profile):'
    printf '    export PATH="%s:$PATH"\n' "${INSTALL_DIR}"
    say 'then reload your shell, e.g. run: . ~/.bashrc'
}

# ---------------------------------------------------------------------------
# Main flow.
# ---------------------------------------------------------------------------
normalize_base_url
detect_target
create_work_dir
resolve_version
prepare_install_dir

ASSET="rswfn-${VERSION}-${TARGET}.tar.gz"
ARCHIVE="${WORK_DIR}/${ASSET}"

if [ "${BASE_MODE}" = 'remote' ]; then
    say "downloading ${ASSET} for version ${VERSION}"
else
    say "copying ${ASSET} from local directory ${BASE_PATH}"
fi
if ! download_asset "${ASSET}" "${ARCHIVE}"; then
    if [ "${BASE_MODE}" = 'remote' ]; then
        error_line "failed to download '${ASSET}' from '${BASE_URL}/v${VERSION}/'"
        die "there is no release asset for target '${TARGET}' (HTTP 404?) or the network is unavailable: check --version '${VERSION}' and pass another --target, e.g. --target x86_64-unknown-linux-gnu"
    fi
    die "archive not found: '${BASE_PATH}/v${VERSION}/${ASSET}': check --version '${VERSION}', --target '${TARGET}' and the local base directory '${BASE_PATH}' (layout: <base>/v<version>/<asset>)"
fi

verify_archive
extract_archive
locate_binary
install_binary
self_test
path_hint

say "rswfn ${VERSION} installed to ${INSTALL_DIR}/rswfn"
say "run 'rswfn --help' to get started"
exit 0
