#!/usr/bin/env python3
"""Download a large HTTP file in parallel byte ranges and verify its MD5."""

import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.request import Request, urlopen


READ_SIZE = 1 << 20


def md5sum(path):
    digest = hashlib.md5()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(READ_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_range(url, start, end, part_path, retries):
    expected = end - start + 1
    for attempt in range(1, retries + 1):
        current = part_path.stat().st_size if part_path.exists() else 0
        if current == expected:
            return part_path
        if current > expected:
            part_path.unlink()
            current = 0
        headers = {
            "Range": f"bytes={start + current}-{end}",
            "Accept-Encoding": "identity",
            "User-Agent": "bmc-multipart-downloader/1.0",
        }
        try:
            curl = shutil.which("curl") or shutil.which("curl.exe")
            if curl:
                temp_path = part_path.with_suffix(part_path.suffix + ".curl")
                command = [
                    curl, "-fsSL", "--connect-timeout", "20",
                    "--speed-limit", "1024", "--speed-time", "20",
                    "--max-time", "120", "--retry", "5", "--retry-delay", "2",
                    "--retry-all-errors", "-r", f"{start + current}-{end}",
                    "-o", str(temp_path), url,
                ]
                try:
                    subprocess.run(command, check=True, timeout=7200)
                    with open(part_path, "ab") as handle, open(temp_path, "rb") as source:
                        for block in iter(lambda: source.read(READ_SIZE), b""):
                            handle.write(block)
                finally:
                    temp_path.unlink(missing_ok=True)
            else:
                request = Request(url, headers=headers)
                with urlopen(request, timeout=180) as response:
                    if current and response.status != 206:
                        raise RuntimeError("server ignored the resume range")
                    mode = "ab" if current else "wb"
                    with open(part_path, mode) as handle:
                        while True:
                            chunk = response.read(READ_SIZE)
                            if not chunk:
                                break
                            handle.write(chunk)
            if part_path.stat().st_size == expected:
                return part_path
        except Exception as exc:
            if attempt == retries:
                raise RuntimeError(f"range {start}-{end} failed: {exc}") from exc
            time.sleep(min(30, attempt * 2))
    raise RuntimeError(f"range {start}-{end} did not complete")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--size", type=int, required=True)
    parser.add_argument("--md5", help="Expected MD5; optional for SRA ODP objects")
    parser.add_argument("--parts", type=int, default=16)
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=0,
        help="Use fixed-size chunks instead of --parts when greater than zero",
    )
    parser.add_argument("--retries", type=int, default=20)
    args = parser.parse_args()

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    if output.exists() and output.stat().st_size == args.size:
        actual = md5sum(output)
        if not args.md5 or actual.lower() == args.md5.lower():
            print(f"[SKIP] verified existing file: {output}")
            return 0

    part_dir = output.parent / f".{output.name}.parts"
    part_dir.mkdir(parents=True, exist_ok=True)
    ranges = []
    if args.chunk_size > 0:
        index = 0
        for start in range(0, args.size, args.chunk_size):
            end = min(args.size - 1, start + args.chunk_size - 1)
            ranges.append((index, start, end, part_dir / f"part_{index:06d}"))
            index += 1
    else:
        part_size = (args.size + args.parts - 1) // args.parts
        for index in range(args.parts):
            start = index * part_size
            end = min(args.size - 1, start + part_size - 1)
            if start <= end:
                ranges.append((index, start, end, part_dir / f"part_{index:03d}"))

    print(f"[DOWNLOAD] {output} ({args.size / (1 << 30):.2f} GiB, {len(ranges)} parts)")
    completed = 0
    with ThreadPoolExecutor(max_workers=len(ranges)) as pool:
        futures = {
            pool.submit(download_range, args.url, start, end, part_path, args.retries): part_path
            for _, start, end, part_path in ranges
        }
        for future in as_completed(futures):
            future.result()
            completed += 1
            print(f"[PART] {completed}/{len(ranges)} complete", flush=True)

    merged = output.with_name(output.name + ".merge")
    with open(merged, "wb") as out_handle:
        for _, _, _, part_path in ranges:
            with open(part_path, "rb") as in_handle:
                for chunk in iter(lambda: in_handle.read(READ_SIZE), b""):
                    out_handle.write(chunk)

    if merged.stat().st_size != args.size:
        raise SystemExit(
            f"Merged size mismatch: {merged.stat().st_size} != {args.size}"
        )
    actual = md5sum(merged)
    if args.md5 and actual.lower() != args.md5.lower():
        raise SystemExit(f"MD5 mismatch: {actual} != {args.md5}")

    os.replace(merged, output)
    for part_path in part_dir.glob("part_*"):
        part_path.unlink(missing_ok=True)
    part_dir.rmdir()
    print(f"[OK] {output} md5={actual}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
