#!/usr/bin/env python3
from __future__ import annotations

import argparse
import ast
import multiprocessing
import os
import re
import subprocess
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXE = ROOT / "build" / ("cpp_com_link_gen.exe" if os.name == "nt" else "cpp_com_link_gen")
PD_HEADER_RE = re.compile(r"^//\s*PD_CODE:\s*(.*)$")


def one_line_error(exc: BaseException | str) -> str:
    if isinstance(exc, BaseException):
        text = f"{exc.__class__.__name__}: {exc}"
    else:
        text = str(exc)
    text = " ".join(text.replace("\r", " ").replace("\n", " ").split())
    if not text:
        text = "unknown error"
    return f"ERROR[{text}]"


def ensure_built() -> None:
    if EXE.is_file():
        return
    subprocess.check_call([sys.executable, str(ROOT / "build.py")], cwd=ROOT)


def parse_pd_code_from_file(path: Path) -> str:
    with path.open("r", encoding="utf-8", errors="replace") as fp:
        for line in fp:
            match = PD_HEADER_RE.match(line)
            if match:
                # Normalize the header through Python's parser so equivalent
                # list formatting is emitted consistently to the C++ parser.
                return repr(ast.literal_eval(match.group(1).strip()))
    raise ValueError(f"PD_CODE header not found: {path}")


def txt_sort_key(path: Path) -> tuple[int, int | str, str]:
    if path.stem.isdigit():
        return (0, int(path.stem), path.name)
    return (1, path.name, path.name)


def selected_txt_files(
    data_dir: Path,
    *,
    recursive: bool,
    numeric_only: bool,
    start_index: int | None,
    end_index: int | None,
    limit: int | None,
) -> list[Path]:
    iterator = data_dir.rglob("*.txt") if recursive else data_dir.glob("*.txt")
    paths = sorted((path for path in iterator if path.is_file()), key=txt_sort_key)
    if numeric_only:
        paths = [path for path in paths if path.stem.isdigit()]
    if start_index is not None:
        paths = [path for path in paths if path.stem.isdigit() and int(path.stem) >= start_index]
    if end_index is not None:
        paths = [path for path in paths if path.stem.isdigit() and int(path.stem) <= end_index]
    if limit is not None:
        paths = paths[:limit]
    return paths


def output_label(data_dir: Path, path: Path, recursive: bool) -> str:
    if recursive:
        return path.relative_to(data_dir).as_posix()
    return path.name


def compute_one(task: tuple[int, str, str, str, int]) -> dict[str, object]:
    index, path_text, label, exe_text, timeout = task
    path = Path(path_text)
    try:
        pd_code = parse_pd_code_from_file(path)
        result = subprocess.run(
            [exe_text, "kh", "--pd", pd_code],
            cwd=ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            check=False,
        )
        if result.returncode != 0:
            message = result.stderr.strip() or result.stdout.strip() or f"exit code {result.returncode}"
            homology = one_line_error(message)
            ok = False
        else:
            lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]
            homology = " || ".join(lines) if lines else one_line_error("empty cppkh output")
            ok = bool(lines)
        return {
            "index": index,
            "label": label,
            "path": path_text,
            "homology": homology,
            "ok": ok,
        }
    except Exception as exc:
        return {
            "index": index,
            "label": label,
            "path": path_text,
            "homology": one_line_error(exc),
            "ok": False,
        }


def export_cppkh_khovanov(
    data_dir: Path,
    output_path: Path,
    *,
    workers: int,
    chunksize: int,
    timeout: int,
    recursive: bool,
    numeric_only: bool,
    start_index: int | None,
    end_index: int | None,
    limit: int | None,
    progress_every: int,
) -> dict[str, object]:
    ensure_built()
    data_dir = data_dir.resolve()
    output_path = output_path.resolve()
    paths = selected_txt_files(
        data_dir,
        recursive=recursive,
        numeric_only=numeric_only,
        start_index=start_index,
        end_index=end_index,
        limit=limit,
    )
    if not paths:
        raise SystemExit(f"no selected .txt files under {data_dir}")

    worker_count = max(1, int(workers))
    tasks = [
        (index, str(path), output_label(data_dir, path, recursive), str(EXE), int(timeout))
        for index, path in enumerate(paths)
    ]
    results: list[dict[str, object] | None] = [None] * len(tasks)
    started_at = time.time()
    checked = 0
    errors = 0

    print(f"cppkh export starting: files={len(tasks)} workers={worker_count} output={output_path}")
    with multiprocessing.Pool(processes=worker_count) as pool:
        for item in pool.imap_unordered(compute_one, tasks, chunksize=max(1, int(chunksize))):
            checked += 1
            results[int(item["index"])] = item
            if not item["ok"]:
                errors += 1
            if progress_every and (
                checked == 1 or checked % progress_every == 0 or checked == len(tasks)
            ):
                elapsed = time.time() - started_at
                speed = checked / elapsed if elapsed > 0 else 0.0
                print(f"  computed {checked}/{len(tasks)} errors={errors} speed={speed:.2f}/s")

    missing = [index for index, item in enumerate(results) if item is None]
    if missing:
        raise RuntimeError(f"missing worker result indices: {missing[:10]}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = output_path.with_name(output_path.name + ".tmp")
    with temp_path.open("w", encoding="utf-8", newline="\n") as fp:
        for item in results:
            assert item is not None
            fp.write(f"{item['label']}: {item['homology']}\n")
    os.replace(temp_path, output_path)

    elapsed = time.time() - started_at
    print(f"cppkh export written: files={len(results)} errors={errors} output={output_path} elapsed={elapsed:.1f}s")
    return {
        "checked": len(results),
        "errors": errors,
        "output_path": str(output_path),
        "elapsed_seconds": elapsed,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Extract PD_CODE headers from generated txt files and export cppkh Khovanov values."
    )
    parser.add_argument("--dir", type=Path, required=True, help="directory containing generated .txt files")
    parser.add_argument("--output", type=Path, required=True, help="output text file")
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    parser.add_argument("--chunksize", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=600, help="seconds per file")
    parser.add_argument("--recursive", action="store_true")
    parser.add_argument("--numeric-only", action="store_true")
    parser.add_argument("--start-index", type=int)
    parser.add_argument("--end-index", type=int)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--progress-every", type=int, default=25)
    args = parser.parse_args()

    export_cppkh_khovanov(
        args.dir,
        args.output,
        workers=args.workers,
        chunksize=args.chunksize,
        timeout=args.timeout,
        recursive=args.recursive,
        numeric_only=args.numeric_only,
        start_index=args.start_index,
        end_index=args.end_index,
        limit=args.limit,
        progress_every=args.progress_every,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
