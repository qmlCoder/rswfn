#!/usr/bin/env python3
"""把发布资产上传到阿里云 OSS。

布局与 install.sh / install.ps1 的镜像约定一致（见 docs/install.md 第 5 节）：

    oss://wmview/rswfn/v<version>/<asset>
    oss://wmview/rswfn/v<version>/SHA256SUMS

所以可以直接把 OSS 的公开访问地址当作下载基址：

    sh install.sh --version <version> \\
        --base-url https://<bucket>.<endpoint>/rswfn

环境变量：
    OSS_ACCESS_KEY_ID / OSS_ACCESS_KEY_SECRET  AccessKey（必填）
    OSS_ENDPOINT                               例如 oss-cn-hangzhou.aliyuncs.com
    OSS_BUCKET                                 默认 wmview
    OSS_PREFIX                                 默认 rswfn

用法：
    python3 .github/scripts/upload_oss.py <version> <dist-dir>
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

DEFAULT_BUCKET = "wmview"
DEFAULT_PREFIX = "rswfn"

NAME_TYPES = {"SHA256SUMS": "text/plain; charset=utf-8"}
SUFFIX_TYPES = {
    ".zip": "application/zip",
    ".gz": "application/gzip",
    ".sha256": "text/plain; charset=utf-8",
}


def content_type(path: Path) -> str:
    if path.name in NAME_TYPES:
        return NAME_TYPES[path.name]
    return SUFFIX_TYPES.get(path.suffix, "application/octet-stream")


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print("用法: upload_oss.py <version> <dist-dir>", file=sys.stderr)
        return 2

    raw_version, dist = argv[1], Path(argv[2])
    # 允许调用方直接传 tag 名（v0.1.3），统一去掉前缀
    version = raw_version[1:] if raw_version.startswith("v") else raw_version

    missing = [
        name
        for name in ("OSS_ACCESS_KEY_ID", "OSS_ACCESS_KEY_SECRET", "OSS_ENDPOINT")
        if not os.environ.get(name)
    ]
    if missing:
        print("缺少环境变量: " + ", ".join(missing), file=sys.stderr)
        return 2

    if not dist.is_dir():
        print(f"目录不存在: {dist}", file=sys.stderr)
        return 1
    files = sorted(path for path in dist.iterdir() if path.is_file())
    if not files:
        print(f"{dist} 里没有文件", file=sys.stderr)
        return 1

    endpoint = os.environ["OSS_ENDPOINT"].strip()
    if not endpoint.startswith(("http://", "https://")):
        endpoint = "https://" + endpoint
    bucket_name = os.environ.get("OSS_BUCKET") or DEFAULT_BUCKET
    prefix = (os.environ.get("OSS_PREFIX") or DEFAULT_PREFIX).strip("/")

    import oss2

    auth = oss2.Auth(os.environ["OSS_ACCESS_KEY_ID"], os.environ["OSS_ACCESS_KEY_SECRET"])
    bucket = oss2.Bucket(auth, endpoint, bucket_name)

    for path in files:
        key = f"{prefix}/v{version}/{path.name}"
        bucket.put_object_from_file(
            key,
            str(path),
            headers={
                "Content-Type": content_type(path),
                # 版本目录的内容不会再变，允许长期缓存
                "Cache-Control": "public, max-age=31536000, immutable",
            },
        )
        print(f"uploaded oss://{bucket_name}/{key} ({path.stat().st_size} bytes)")

    host = endpoint.split("://", 1)[-1].rstrip("/")
    print(f"done: {len(files)} object(s) under oss://{bucket_name}/{prefix}/v{version}/")
    print(f"install base url: https://{bucket_name}.{host}/{prefix}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
