#!/usr/bin/env python3
"""Write SHA-256 and size inventories for a reproducible artifact directory."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--json", type=Path)
    parser.add_argument("--sums", type=Path)
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    json_path = (args.json or root / "artifact_inventory.json").resolve()
    sums_path = (args.sums or root / "SHA256SUMS.txt").resolve()
    excluded = {json_path, sums_path}
    entries = []
    for path in sorted(root.rglob("*")):
        resolved = path.resolve()
        if not path.is_file() or resolved in excluded:
            continue
        entries.append(
            {
                "path": path.relative_to(root).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": _sha256(path),
            }
        )
    report = {
        "schema": "artifact-inventory-v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "root_name": root.name,
        "file_count": len(entries),
        "total_bytes": sum(entry["bytes"] for entry in entries),
        "files": entries,
    }
    json_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    sums_path.write_text(
        "".join(f"{entry['sha256']}  {entry['path']}\n" for entry in entries),
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "file_count": report["file_count"],
                "total_bytes": report["total_bytes"],
                "json": str(json_path),
                "sums": str(sums_path),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
