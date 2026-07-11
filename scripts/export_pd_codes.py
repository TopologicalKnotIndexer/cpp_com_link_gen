#!/usr/bin/env python3
from __future__ import annotations

import argparse
import ast
import os
import re
import time
from pathlib import Path


PD_HEADER_RE = re.compile(r"^//\s*PD_CODE:\s*(.*)$")


def one_line_error(exc: BaseException | str) -> str:
    if isinstance(exc, BaseException):
        text = f"{exc.__class__.__name__}: {exc}"
    else:
        text = str(exc)
    text = " ".join(text.replace("\r", " ").replace("\n", " ").split())
    return f"ERROR[{text or 'unknown error'}]"


def parse_pd_code_from_file(path: Path) -> str:
    with path.open("r", encoding="utf-8", errors="replace") as fp:
        for line in fp:
            match = PD_HEADER_RE.match(line)
            if match:
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


def export_pd_codes(
    data_dir: Path,
    output_path: Path,
    *,
    recursive: bool,
    numeric_only: bool,
    start_index: int | None,
    end_index: int | None,
    limit: int | None,
    progress_every: int,
) -> dict[str, object]:
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

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = output_path.with_name(output_path.name + ".tmp")

    started_at = time.time()
    errors = 0
    with temp_path.open("w", encoding="utf-8", newline="\n") as fp:
        for index, path in enumerate(paths, start=1):
            try:
                pd_code = parse_pd_code_from_file(path)
            except Exception as exc:
                pd_code = one_line_error(exc)
                errors += 1
            fp.write(f"{output_label(data_dir, path, recursive)}: {pd_code}\n")

            if progress_every and (
                index == 1 or index % progress_every == 0 or index == len(paths)
            ):
                elapsed = time.time() - started_at
                speed = index / elapsed if elapsed > 0 else 0.0
                print(f"  extracted {index}/{len(paths)} errors={errors} speed={speed:.2f}/s")

    os.replace(temp_path, output_path)
    elapsed = time.time() - started_at
    print(f"PD_CODE export written: files={len(paths)} errors={errors} output={output_path} elapsed={elapsed:.1f}s")
    return {
        "checked": len(paths),
        "errors": errors,
        "output_path": str(output_path),
        "elapsed_seconds": elapsed,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Extract PD_CODE headers from generated txt files."
    )
    parser.add_argument("--dir", type=Path, required=True, help="directory containing generated .txt files")
    parser.add_argument("--output", type=Path, required=True, help="output text file")
    parser.add_argument("--recursive", action="store_true")
    parser.add_argument("--numeric-only", action="store_true")
    parser.add_argument("--start-index", type=int)
    parser.add_argument("--end-index", type=int)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--progress-every", type=int, default=1000)
    args = parser.parse_args()

    export_pd_codes(
        args.dir,
        args.output,
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
