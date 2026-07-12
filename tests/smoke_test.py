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


def run_failure(args: list[str], cwd: Path = ROOT) -> str:
    result = subprocess.run(
        [str(EXE), *args],
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    assert result.returncode != 0
    return result.stderr


def main() -> int:
    if not EXE.is_file():
        subprocess.check_call([sys.executable, str(ROOT / "build.py")], cwd=ROOT)

    trefoil = "[[1,5,2,4],[3,1,4,6],[5,3,6,2]]"
    kh = run(["kh", "--pd", trefoil]).splitlines()
    assert len(kh) == 1
    assert any("q^1*t^0*Z[0]" in line for line in kh)
    relabeled_trefoil = "[[10,40,30,60],[20,10,60,50],[40,20,50,30]]"
    assert run(["kh", "--pd", relabeled_trefoil]).splitlines() == kh
    assert "twice" in run_failure(["kh", "--pd", "[[1,2,3,4]]"]).lower()

    hopf = "[[2,3,1,4],[4,1,3,2]]"
    hopf_kh = run(["kh", "--pd", hopf]).splitlines()
    assert len(hopf_kh) == 1
    hopf_all = run(["kh-all-orientations", "--pd", hopf]).splitlines()
    assert len(hopf_all) == 2

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

        run(["khovanov", "--dir", str(generated), "--jobs", "2"])
        run(["postprocess", "--dir", str(generated), "--jobs", "2"])
        cluster_dir = data_root / "cluster"
        re_cluster_dir = data_root / "re_cluster"
        assert cluster_dir.is_dir()
        assert re_cluster_dir.is_dir()
        assert list(cluster_dir.rglob("*.svg"))

        run(["process-one", str(files[0])])
        content = files[0].read_text(encoding="utf-8")
        assert "// KHOVANOV:" in content
        assert "// PD_CODE:" in content

        manual_dir = tmp_path / "manual_kh"
        manual_dir.mkdir()
        manual_pd = "[[2, 3, 1, 4], [4, 1, 3, 2]]"
        for index, z_text in enumerate(("Z[2,0]", "Z[0,2]"), start=1):
            (manual_dir / f"{index:07d}.txt").write_text(
                f"// KHOVANOV: q^0*t^0*{z_text}\n"
                f"// PD_CODE: {manual_pd}\n"
                "[L2a1]\n",
                encoding="utf-8",
            )
        classify_out = run([
            "classify",
            "--dir",
            str(manual_dir),
            "--cluster-dir",
            str(tmp_path / "manual_cluster"),
            "--dry-run",
        ])
        assert "ignore_pd_code: 1" in classify_out
        assert "khovanov_classes: 1" in classify_out

    print("smoke test passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
