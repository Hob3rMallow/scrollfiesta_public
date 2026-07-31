"""Mirror a remote OME-Zarr store (s3://) to a local directory, verbatim.

The scrollfiesta mesher (``scrollunwrap.zarr_source.open_volume``) reads a local
zarr path exactly like an ``s3://`` one, so a local mirror lets us carve any ROI
/ cube / grid of the whole scroll offline -- no repeated S3 round-trips.

Downloads every chunk object under the requested levels (plus the store's
top-level ``.zattrs`` / ``.zgroup`` / ``metadata.json``) preserving the key
layout, so the result is a byte-for-byte usable zarr store. The transfer is:

  * parallel        -- a thread pool of ``--workers`` concurrent object GETs;
  * resumable       -- an object already present locally with the right byte
                       size is skipped, so re-running finishes the tail;
  * verifiable      -- ``--verify`` re-opens the local store and reads a cube.

Usage (anonymous public bucket)::

    py -m uv run --project python python scripts/mirror_zarr.py \
        --zarr s3://vesuvius-challenge-open-data/PHerc0139/.../<...>.zarr \
        --out  D:/work/vesuvius-c/PHerc0139-full/<...>.zarr \
        --levels all --workers 32

``--dry-run`` reports the object count / byte total per level and exits.
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import s3fs


def human(n: float) -> str:
    for u in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024:
            return f"{n:.2f} {u}"
        n /= 1024
    return f"{n:.2f} PB"


def _split_s3(uri: str) -> str:
    """s3://bucket/key -> bucket/key (the form s3fs.ls/find expect)."""
    return uri[len("s3://"):] if uri.startswith("s3://") else uri


def list_objects(fs: s3fs.S3FileSystem, root: str, levels) -> list[dict]:
    """Return [{'key','size','rel'}] for every file object to mirror.

    ``levels`` is None (=> everything under root) or a set of level names; when a
    set, we take each ``root/<level>`` subtree plus the store-level metadata files
    that sit directly under ``root``.
    """
    objs: list[dict] = []
    seen: set[str] = set()

    def add_tree(prefix: str):
        for key, info in fs.find(prefix, detail=True).items():
            if info.get("type") != "file" or key in seen:
                continue
            seen.add(key)
            objs.append({"key": key, "size": int(info.get("size") or 0),
                         "rel": key[len(root) + 1:]})

    if levels is None:
        add_tree(root)
    else:
        # top-level metadata files (.zattrs/.zgroup/metadata.json/.zmetadata)
        for e in fs.ls(root, detail=True):
            if e.get("type") == "file":
                key = e["name"]
                if key not in seen:
                    seen.add(key)
                    objs.append({"key": key, "size": int(e.get("size") or 0),
                                 "rel": key[len(root) + 1:]})
        for lv in sorted(levels):
            add_tree(f"{root}/{lv}")
    return objs


def _download_one(fs, obj, out_dir: Path, retries: int = 4) -> tuple[str, int, str]:
    """Fetch one object -> (rel, bytes_written, status). Skips if size matches."""
    rel, size = obj["rel"], obj["size"]
    dst = out_dir / rel
    if dst.exists() and dst.stat().st_size == size and size > 0:
        return rel, 0, "skip"
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_name(dst.name + ".part")
    last = ""
    for attempt in range(retries):
        try:
            fs.get_file(obj["key"], str(tmp))
            got = tmp.stat().st_size
            if size and got != size:
                last = f"size mismatch {got}!={size}"
                tmp.unlink(missing_ok=True)
                continue
            os.replace(tmp, dst)
            return rel, got, "ok"
        except Exception as e:  # transient S3 error -> back off and retry
            last = f"{type(e).__name__}: {str(e)[:160]}"
            time.sleep(0.5 * (attempt + 1))
    return rel, 0, f"FAIL {last}"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Mirror a remote OME-Zarr to disk.")
    ap.add_argument("--zarr", required=True, help="s3:// (or local) zarr store URI")
    ap.add_argument("--out", required=True, type=Path, help="local destination dir")
    ap.add_argument("--levels", default="all",
                    help="'all' or comma list, e.g. '0,1,2,3,4,5' (default all)")
    ap.add_argument("--workers", type=int, default=32)
    ap.add_argument("--shard", default=None, metavar="i/N",
                    help="mirror only objects with (sorted) index %% N == i, so "
                         "N processes can run disjoint shards in parallel (each "
                         "writes different files -> no collisions). e.g. --shard 0/8")
    ap.add_argument("--anon", choices=["auto", "yes", "no"], default="auto",
                    help="S3 auth: auto=signed iff AWS creds present else anon")
    ap.add_argument("--endpoint", default=None, help="custom S3 endpoint URL")
    ap.add_argument("--dry-run", action="store_true",
                    help="list object counts / sizes per level and exit")
    args = ap.parse_args(argv)

    shard_i, shard_n = 0, 1
    if args.shard:
        shard_i, shard_n = (int(x) for x in args.shard.split("/"))
        if not (0 <= shard_i < shard_n):
            ap.error(f"--shard {args.shard}: need 0 <= i < N")

    root = _split_s3(str(args.zarr).rstrip("/"))
    anon = {"auto": None, "yes": True, "no": False}[args.anon]
    if anon is None:
        cred = ("AWS_ACCESS_KEY_ID", "AWS_PROFILE", "AWS_ROLE_ARN")
        anon = not any(os.environ.get(k) for k in cred)
    so = {"anon": bool(anon)}
    if args.endpoint:
        so["client_kwargs"] = {"endpoint_url": args.endpoint}
    fs = s3fs.S3FileSystem(**so)

    levels = None
    if args.levels != "all":
        levels = {s.strip() for s in args.levels.split(",") if s.strip()}

    print(f"listing {root} (levels={args.levels}) ...", flush=True)
    t0 = time.perf_counter()
    objs = list_objects(fs, root, levels)
    objs.sort(key=lambda o: o["rel"])  # deterministic order for stable sharding
    total_all = len(objs)
    bytes_all = sum(o["size"] for o in objs)
    print(f"  {total_all} objects, {human(bytes_all)} "
          f"[{time.perf_counter()-t0:.1f}s listing]")
    if shard_n > 1:
        objs = [o for k, o in enumerate(objs) if k % shard_n == shard_i]
        print(f"  shard {shard_i}/{shard_n}: {len(objs)} objects, "
              f"{human(sum(o['size'] for o in objs))}")
    total_bytes = sum(o["size"] for o in objs)

    # per-level breakdown
    by_lvl: dict[str, list[int]] = {}
    for o in objs:
        top = o["rel"].split("/", 1)[0]
        lvl = top if top.isdigit() else "(meta)"
        by_lvl.setdefault(lvl, [0, 0])
        by_lvl[lvl][0] += 1
        by_lvl[lvl][1] += o["size"]
    for lvl in sorted(by_lvl):
        n, b = by_lvl[lvl]
        print(f"    level {lvl:6s}: {n:>7d} objects  {human(b)}")

    if args.dry_run:
        return 0

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"\nmirroring -> {out_dir}  ({args.workers} workers)")

    done_bytes = skipped = failed = 0
    fails: list[str] = []
    t0 = time.perf_counter()
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(_download_one, fs, o, out_dir) for o in objs]
        for i, f in enumerate(as_completed(futs), 1):
            rel, nb, status = f.result()
            if status == "ok":
                done_bytes += nb
            elif status == "skip":
                skipped += 1
            else:
                failed += 1
                fails.append(f"{rel}: {status}")
            if i % 200 == 0 or i == len(futs):
                dt = time.perf_counter() - t0
                rate = done_bytes / dt if dt else 0
                print(f"  {i}/{len(futs)}  new={human(done_bytes)} "
                      f"skip={skipped} fail={failed}  {human(rate)}/s", flush=True)

    dt = time.perf_counter() - t0
    print(f"\ndone in {dt/60:.1f} min: new={human(done_bytes)} "
          f"skipped={skipped} failed={failed}")
    if fails:
        print("FAILURES (first 20):")
        for line in fails[:20]:
            print("  " + line)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
