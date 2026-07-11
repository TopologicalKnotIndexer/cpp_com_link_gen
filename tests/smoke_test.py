#!/usr/bin/env python3
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EXE = ROOT / "build" / ("cpp_com_link_gen.exe" if os.name == "nt" else "cpp_com_link_gen")


def run(args: list[str], cwd: Path = ROOT) -> str:
    result = subprocess.run(
        [str(EXE), *args],
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    return result.stdout


def main() -> int:
    if not EXE.is_file():
        subprocess.check_call([sys.executable, str(ROOT / "build.py")], cwd=ROOT)

    trefoil = "[[1,5,2,4],[3,1,4,6],[5,3,6,2]]"
    kh = run(["kh", "--pd", trefoil]).splitlines()
    assert len(kh) == 1
    assert any("q^1*t^0*Z[0]" in line for line in kh)

    with tempfile.TemporaryDirectory(prefix="cpp_com_link_gen_test_") as tmp:
        tmp_path = Path(tmp)
        data_root = tmp_path / "data"
        shutil.copytree(ROOT / "data" / "prime_link_knot_10", data_root / "prime_link_knot_10")

        out = run([
            "generate",
            "--total-crs",
            "3",
            "--max-prime-cnt",
            "1",
            "--jobs",
            "2",
            "--data-root",
            str(data_root),
        ])
        generated = Path(out.strip().splitlines()[-1].strip('"'))
        files = sorted(generated.glob("*.txt"))
        assert len(files) == 4

        run(["process-one", str(files[0])])
        content = files[0].read_text(encoding="utf-8")
        assert "// KHOVANOV:" in content
        assert "// PD_CODE:" in content

    print("smoke test passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
