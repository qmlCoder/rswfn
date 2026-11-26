#!/usr/bin/env python3
"""把发布资产上传到阿里云 OSS。

布局与 install.sh / install.ps1 的镜像约定一致（见 docs/install.md 第 5 节）：

    oss://wmview/rswfn/v<version>/<asset>     版本目录（不可变）
    oss://wmview/rswfn/v<version>/SHA256SUMS
    oss://wmview/rswfn/install.sh             安装脚本（默认基址指向本镜像）
    oss://wmview/rswfn/install.ps1
    oss://wmview/rswfn/latest.txt             一行版本号，用于解析“最新版本”

有了后三个，就可以直接从 OSS 一行安装：

    curl -fsSL https://<bucket>.<endpoint>/rswfn/install.sh | sh
    irm https://<bucket>.<endpoint>/rswfn/install.ps1 | iex

环境变量：
    OSS_ACCESS_KEY_ID / OSS_ACCESS_KEY_SECRET  AccessKey（必填）
    OSS_ENDPOINT                               例如 oss-cn-hangzhou.aliyuncs.com
    OSS_BUCKET                                 默认 wmview
    OSS_PREFIX                                 默认 rswfn

取值顺序：进程环境变量 > 仓库根目录的 .env（模板见 .env.example）。
所以 CI 里靠 GitHub Secrets 提供的环境变量即可，本地则把凭据放在 .env。

用法：
    python3 scripts/upload_oss.py <version> [dist-dir] [--installers] [--latest]

    --installers  同时发布 install.sh / install.ps1（默认下载基址改写成本镜像）
    --latest      同时更新 latest.txt 为 <version>

不带 dist-dir 时只发布 installers / latest.txt，不碰版本目录。
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

DEFAULT_BUCKET = "wmview"
DEFAULT_PREFIX = "rswfn"

# 版本目录内容不会再变，可以长期缓存；根部文件（installers、latest.txt）会随发布更新
IMMUTABLE_CACHE = "public, max-age=31536000, immutable"
RECENT_CACHE = "public, max-age=60"

NAME_TYPES = {"SHA256SUMS": "text/plain; charset=utf-8", "latest.txt": "text/plain; charset=utf-8"}
SUFFIX_TYPES = {
    ".zip": "application/zip",
    ".gz": "application/gzip",
    ".sha256": "text/plain; charset=utf-8",
    ".sh": "text/plain; charset=utf-8",
    ".ps1": "text/plain; charset=utf-8",
    ".txt": "text/plain; charset=utf-8",
}

# 安装脚本里需要改成本镜像地址的字面量
SH_BASE_PREFIX = "DEFAULT_BASE_URL="
PS1_BASE_PREFIX = "$BaseUrl = "
GITHUB_RAW = "https://raw.githubusercontent.com/qmlCoder/rswfn/main"


def content_type(path: Path) -> str:
    if path.name in NAME_TYPES:
        return NAME_TYPES[path.name]
    return SUFFIX_TYPES.get(path.suffix, "application/octet-stream")


def load_dotenv() -> Path | None:
    """把 .env 里的键补进环境变量（已存在的不覆盖，所以 CI 的 secrets 优先）。"""
    repo_root = Path(__file__).resolve().parents[1]
    for candidate in (Path.cwd() / ".env", repo_root / ".env"):
        if not candidate.is_file():
            continue
        # utf-8-sig：容忍 Windows 记事本写入的 BOM
        for line in candidate.read_text(encoding="utf-8-sig").splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("export "):
                line = line[len("export ") :].lstrip()
            key, sep, value = line.partition("=")
            if not sep:
                continue
            key = key.strip()
            value = value.strip().strip('"').strip("'")
            if key and key not in os.environ:
                os.environ[key] = value
        return candidate
    return None


def rewrite_installers(base_url: str) -> dict[str, bytes]:
    """把安装脚本的默认下载基址改写成本镜像，返回 {对象名: 内容}。"""
    repo_root = Path(__file__).resolve().parents[1]

    sh_line = f"{SH_BASE_PREFIX}'{base_url}'"
    ps1_line = f"{PS1_BASE_PREFIX}'{base_url}'"

    def rewrite(name: str, is_default_line, replacement: str) -> list[str]:
        text = (repo_root / name).read_text(encoding="utf-8-sig")
        text = text.replace("\r\n", "\n")
        out = []
        for line in text.splitlines(keepends=True):
            if is_default_line(line):
                indent = line[: len(line) - len(line.lstrip())]
                out.append(indent + replacement + "\n")
            else:
                out.append(line)
        # 脚本里“先下载再执行”的提示也指向本镜像
        joined = "".join(out).replace(f"{GITHUB_RAW}/{name}", f"{base_url}/{name}")
        return joined

    sh_text = rewrite("install.sh", lambda line: line.startswith(SH_BASE_PREFIX), sh_line)
    # PowerShell 习惯 CRLF
    ps1_text = rewrite("install.ps1", lambda line: line.strip().startswith(PS1_BASE_PREFIX), ps1_line)
    ps1_text = ps1_text.replace("\n", "\r\n")

    return {"install.sh": sh_text.encode("utf-8"), "install.ps1": ps1_text.encode("utf-8")}


def main(argv: list[str]) -> int:
    version = None
    dist = None
    with_installers = False
    with_latest = False

    for arg in argv[1:]:
        if arg == "--installers":
            with_installers = True
        elif arg == "--latest":
            with_latest = True
        elif arg.startswith("-"):
            print(f"未知参数: {arg}", file=sys.stderr)
            return 2
        elif version is None:
            version = arg
        elif dist is None:
            dist = arg
        else:
            print(f"多余参数: {arg}", file=sys.stderr)
            return 2

    if version is None:
        print("用法: upload_oss.py <version> [dist-dir] [--installers] [--latest]", file=sys.stderr)
        return 2

    # 允许调用方直接传 tag 名（v0.1.3），统一去掉前缀
    version = version[1:] if version.startswith("v") else version

    env_file = load_dotenv()
    if env_file is not None:
        print(f"read defaults from {env_file}")

    missing = [
        name
        for name in ("OSS_ACCESS_KEY_ID", "OSS_ACCESS_KEY_SECRET", "OSS_ENDPOINT")
        if not os.environ.get(name)
    ]
    if missing:
        print("缺少环境变量: " + ", ".join(missing), file=sys.stderr)
        return 2

    files: list[Path] = []
    if dist is not None:
        dist_path = Path(dist)
        if not dist_path.is_dir():
            print(f"目录不存在: {dist_path}", file=sys.stderr)
            return 1
        files = sorted(path for path in dist_path.iterdir() if path.is_file())
        if not files:
            print(f"{dist_path} 里没有文件", file=sys.stderr)
            return 1

    if not files and not with_installers and not with_latest:
        print("没有要上传的内容：给出 dist 目录，或加 --installers / --latest", file=sys.stderr)
        return 1

    endpoint = os.environ["OSS_ENDPOINT"].strip()
    if not endpoint.startswith(("http://", "https://")):
        endpoint = "https://" + endpoint
    bucket_name = os.environ.get("OSS_BUCKET") or DEFAULT_BUCKET
    prefix = (os.environ.get("OSS_PREFIX") or DEFAULT_PREFIX).strip("/")
    host = endpoint.split("://", 1)[-1].rstrip("/")
    base_url = f"https://{bucket_name}.{host}/{prefix}"

    import oss2

    auth = oss2.Auth(os.environ["OSS_ACCESS_KEY_ID"], os.environ["OSS_ACCESS_KEY_SECRET"])
    bucket = oss2.Bucket(auth, endpoint, bucket_name)
    uploaded = 0

    def put(key: str, data: bytes, ctype: str, cache: str) -> None:
        nonlocal uploaded
        bucket.put_object(key, data, headers={"Content-Type": ctype, "Cache-Control": cache})
        print(f"uploaded oss://{bucket_name}/{key} ({len(data)} bytes)")
        uploaded += 1

    for path in files:
        put(
            f"{prefix}/v{version}/{path.name}",
            path.read_bytes(),
            content_type(path),
            IMMUTABLE_CACHE,
        )

    if with_installers:
        for name, data in rewrite_installers(base_url).items():
            put(f"{prefix}/{name}", data, content_type(Path(name)), RECENT_CACHE)

    if with_latest:
        put(
            f"{prefix}/latest.txt",
            f"{version}\n".encode("utf-8"),
            "text/plain; charset=utf-8",
            RECENT_CACHE,
        )

    print(f"done: {uploaded} object(s); versioned assets under oss://{bucket_name}/{prefix}/v{version}/")
    if with_installers:
        print(f"  Linux:   curl -fsSL {base_url}/install.sh | sh")
        print(f"  Windows: irm {base_url}/install.ps1 | iex")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
