# cpp_com_link_gen

C++ replacement for the Python `com-link-khovanov-10` pipeline.

It keeps the same generated file format:

```text
// KHOVANOV: ...
// PD_CODE: [[...], ...]
<link representation>
```

The generator uses the bundled `prime-link-knot-10` PD tables, parses the
`link-rep` text format directly, fixes connected-sum/orientation handling in
C++, and computes Khovanov homology through
[`GGN-2015/cppkh`](https://github.com/GGN-2015/cppkh).

## Build

Use the Python build script:

```bash
python build.py
```

Set `CXX` to choose a compiler:

```bash
CXX=clang++ python build.py
```

On Windows with the requested conda environment:

```powershell
conda run -n math_env python build.py
```

The executable is written to `build/cpp_com_link_gen` or
`build/cpp_com_link_gen.exe`.

## Usage

Generate the default 10-crossing dataset:

```bash
build/cpp_com_link_gen generate --total-crs 10 --max-prime-cnt 3 --jobs 8
```

Compute Khovanov homology for the generated files:

```bash
build/cpp_com_link_gen khovanov --dir data/com_link_gen_10-v0.1.0-com_link_gen-10-3 --jobs 8
```

Run the downstream classification pipeline after Khovanov headers are present:

```bash
build/cpp_com_link_gen postprocess --dir data/com_link_gen_10-v0.1.0-com_link_gen-10-3 --jobs 8
```

This replaces the old `get_khovanov.py`, `make_all_diagram.py`, and
`re_cluster.py` scripts. It writes:

- `data/cluster/<md5>/`: one folder per distinct Khovanov key, with duplicate
  `PD_CODE` entries ignored.
- `khovanov.txt`: the normalized Khovanov values for the class.
- `*.svg`: self-contained C++ generated PD/component diagrams for each retained
  representative.
- `data/re_cluster/NNN/<md5>/`: the same classes regrouped by the number of
  retained link representatives in the class.

The individual downstream stages are also available:

```bash
build/cpp_com_link_gen classify --dir data/com_link_gen_10-v0.1.0-com_link_gen-10-3
build/cpp_com_link_gen diagrams --cluster-dir data/cluster --jobs 8 --force
build/cpp_com_link_gen re-cluster --cluster-dir data/cluster
```

Generate and process in one command:

```bash
build/cpp_com_link_gen all --total-crs 10 --max-prime-cnt 3 --jobs 8
```

Compatibility helpers:

```bash
build/cpp_com_link_gen process-one path/to/0000001.txt
build/cpp_com_link_gen legacy --dir path/to/data --mod 16 --res 0
build/cpp_com_link_gen pd --file path/to/link_rep.txt
build/cpp_com_link_gen kh --pd "[[1,5,2,4],[3,1,4,6],[5,3,6,2]]"
```

If run without arguments, the executable prompts for `process_count>>>` and
processes the default generated `10,3` directory, matching the old Python entry
point.

## Notes

- `total_crs` follows the original package behavior: it means "up to this many
  crossings", not exactly this many crossings.
- Connected sum is implemented by cutting oriented component arcs and gluing
  them before canonical renumbering. This avoids the old local-crossing
  replacement bug-prone path and keeps merged-cluster component representatives
  updated for the whole cluster.
- Component orientation enumeration computes all `2^n` component orientations.
  For each orientation the program passes explicit crossing signs to `cppkh`,
  then removes duplicate homology strings. The final distinct count is checked
  against the `2^(n-1)` theoretical upper bound.
- Downstream clustering is implemented in C++ and does not call the old Python
  scripts. The diagram stage writes SVG files instead of the old Pillow PNG
  files so the pipeline stays self-contained and cross-platform.

## Smoke Test

```bash
python tests/smoke_test.py
```

The test uses a temporary data root and does not write generated datasets into
the repository.

## Algorithm Validation

Run the stronger validation suite with:

```bash
python tests/algorithm_validation.py
```

This test intentionally does not treat the old Python implementation as the
source of truth. It checks known Khovanov outputs for the unknot, trefoil, and
Hopf link; independently recomputes small generation counts from the bundled
prime-link table; validates generated PD codes as 2-regular component graphs;
checks connected-sum crossing/component-count invariants; runs Khovanov on the
generated sample set; and enforces the `2^(n-1)` distinct-Khovanov upper bound.

For an additional mature-library cross-check, install Spherogram in the Python
environment and require it during validation:

```powershell
conda run -n math_env python -m pip install spherogram
conda run -n math_env python tests\algorithm_validation.py --require-spherogram
```

When Spherogram is available, the validation suite constructs independent
Spherogram `Link` objects from generated PD codes and verifies component counts
and PD round trips.

## SageMath Khovanov Cross-Check

There is also a SageMath-only checker for recomputing every component
orientation with Sage's `Link.khovanov_homology()` implementation:

```sage
load("sage/check_oriented_khovanov.sage")
check_khovanov_file("data/com_link_gen_10-v0.1.0-com_link_gen-10-3/0000001.txt")
```

For a parallel directory sample:

```sage
summary = check_khovanov_directory_parallel(
    "data/com_link_gen_10-v0.1.0-com_link_gen-10-3",
    limit=20,
    workers=8,
    progress_every=1,
    json_report_path="sage_khovanov_report.json",
)
summary["ok"]
```

For a split full run, use numeric filename ranges. This is useful when checking
the full 7000+ retained PD-code workload across several machines:

```sage
check_khovanov_directory_parallel(
    "data/com_link_gen_10-v0.1.0-com_link_gen-10-3",
    start_index=1,
    end_index=1000,
    workers=8,
    progress_every=25,
    json_report_path="sage_khovanov_0001_1000.json",
)
```

The serial version is still available as `check_khovanov_directory(...)`, but
the parallel version is recommended for real validation. It uses process-based
parallelism, prefers the `fork` start method when Sage provides it, and checks
one generated txt file per worker task.

The Sage checker intentionally enumerates all `2^n` component orientations. For
each orientation it builds an oriented Gauss code, computes integral Khovanov
homology in Sage, converts Sage's abelian-group output into the same canonical
`Z[...]` form used by `cppkh`, and compares the distinct set with the
`KHOVANOV` headers in the generated txt file. It also checks the `2^(n-1)`
distinct-result upper bound.

## References

- `cppkh`: <https://github.com/GGN-2015/cppkh>
- Spherogram link diagrams: <https://github.com/3-manifolds/Spherogram>
- SageMath knot/link functionality and PD-code conventions:
  <https://doc.sagemath.org/html/en/reference/knots/>
- Knot Atlas planar diagrams: <https://katlas.org/wiki/Planar_Diagrams>
