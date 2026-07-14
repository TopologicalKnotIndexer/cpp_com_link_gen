#!/usr/bin/env python3
"""Cross-platform build script for cpp_com_link_gen.

The script intentionally avoids requiring CMake for the wrapper project.  It
uses CMake-style compiler conventions where possible: set CXX to choose a
compiler, and pass --debug for an unoptimized build.
"""

from __future__ import annotations

import argparse
import os
import platform
import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parent
BUILD = ROOT / "build"
THIRD_PARTY_CPPKH = ROOT / "third_party" / "cppkh"
THIRD_PARTY_PD_DIAGRAM = ROOT / "third_party" / "pd_code_to_diagram" / "cpp_src"


def exe_name(name: str) -> str:
    return name + (".exe" if os.name == "nt" else "")


def run(cmd: list[str], cwd: Path | None = None) -> None:
    print("+ " + " ".join(cmd))
    subprocess.check_call(cmd, cwd=str(cwd) if cwd else None)


def require_cppkh() -> None:
    src = THIRD_PARTY_CPPKH / "src" / "main.cpp"
    if src.is_file():
        return
    raise SystemExit(
        f"Tracked dependency source is missing: {src}. "
        "Restore the ordinary third_party/cppkh files from this repository."
    )


def find_compiler() -> str:
    cxx = os.environ.get("CXX")
    if cxx:
        return cxx

    candidates = ["clang++", "g++"]
    if os.name == "nt":
        candidates = ["cl", "clang++", "g++"]
    for candidate in candidates:
        if shutil.which(candidate):
            return candidate
    raise SystemExit("No C++ compiler found. Set CXX to your compiler path.")


def compiler_kind(cxx: str) -> str:
    base = Path(cxx.split()[0]).name.lower()
    if base in {"cl", "cl.exe"}:
        return "msvc"
    return "posix"


def quote_define_path(path: Path) -> str:
    value = str(path).replace("\\", "\\\\").replace('"', '\\"')
    return f'CPP_COM_LINK_GEN_SOURCE_DIR="{value}"'


def build(args: argparse.Namespace) -> Path:
    require_cppkh()
    BUILD.mkdir(parents=True, exist_ok=True)

    cxx = find_compiler()
    kind = compiler_kind(cxx)
    out = BUILD / exe_name("cpp_com_link_gen")
    sources = [
        ROOT / "src" / "main.cpp",
        THIRD_PARTY_CPPKH / "src" / "main.cpp",
    ]

    if kind == "msvc":
        flags = [
            "/std:c++17",
            "/EHsc",
            "/DCPPKH_SHARED_LIBRARY",
            "/DKH_THREAD_BACKEND_WIN32",
            "/D" + quote_define_path(ROOT),
            "/I" + str(THIRD_PARTY_PD_DIAGRAM),
        ]
        flags += ["/Od", "/Zi"] if args.debug else ["/O2", "/DNDEBUG"]
        cmd = [cxx, *flags, *(str(src) for src in sources), f"/Fe:{out}"]
    else:
        flags = [
            "-std=c++17",
            "-DCPPKH_SHARED_LIBRARY",
            f"-D{quote_define_path(ROOT)}",
            "-I",
            str(THIRD_PARTY_PD_DIAGRAM),
            "-Wall",
            "-Wextra",
        ]
        system = platform.system().lower()
        if system == "windows":
            flags.append("-DKH_THREAD_BACKEND_WIN32")
        else:
            flags.append("-DKH_THREAD_BACKEND_PTHREAD")
        flags += ["-O0", "-g"] if args.debug else ["-O3", "-DNDEBUG"]
        libs: list[str] = []
        if system in {"linux", "darwin"}:
            libs.append("-pthread")
        cmd = [cxx, *flags, *(str(src) for src in sources), "-o", str(out), *libs]

    run(cmd)
    print(f"Built {out}")
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description="Build cpp_com_link_gen")
    parser.add_argument("--debug", action="store_true", help="build without optimizations")
    args = parser.parse_args()
    build(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
