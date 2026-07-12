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
build/cpp_com_link_gen kh-all-orientations --pd "[[2,3,1,4],[4,1,3,2]]"
```

Python export helpers:

```bash
python scripts/export_pd_codes.py \
  --dir data/com_link_gen_10-v0.1.0-com_link_gen-10-3 \
  --output build/pd_codes_10_3.txt \
  --numeric-only
```

This writes one ordered line per selected file:

```text
0000001.txt: [[...], ...]
```

To recompute exactly one Khovanov value from each file's `PD_CODE` through the
C++/`cppkh` path and write one ordered line per file:

```bash
python scripts/export_cppkh_pd_khovanov.py \
  --dir data/com_link_gen_10-v0.1.0-com_link_gen-10-3 \
  --output build/cppkh_pd_khovanov_10_3.txt \
  --numeric-only \
  --workers 8
```

Both scripts support `--start-index`, `--end-index`, `--limit`, and
`--recursive`. The cppkh export writes directly to the target file and flushes
ordered result lines as they become available. If a file cannot be parsed or
computed, its output line contains a single-line `ERROR[...]` value in place of
the PD code or homology.

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
- Every C++ path that calls `cppkh` first normalizes the `PD_CODE`, verifies
  contiguous component numbering, checks that every strand pair is adjacent in
  the component cycle, and then passes explicit crossing signs into `cppkh`.
- Khovanov strings are normalized so every `Z[...]` invariant-factor list is
  sorted in ascending numeric order before output or clustering.
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

## SageMath Khovanov Polynomial Export

The Sage helper can export one Sage-computed Khovanov polynomial for every
selected generated `.txt` file. It extracts only `PD_CODE`, lets Sage construct
the link directly from that PD code, computes
`Link(pd).khovanov_polynomial(var1="q", var2="t")`, and writes ordered lines of
the form `filename: polynomial`. Polynomial terms are sorted by ascending
`t` exponent, then ascending `q` exponent. The output path is converted to an
absolute path and printed at startup. Relative paths are resolved from Sage's
current working directory, so use an absolute output path if you want the file
in this repository's `build` directory.

```sage
load("sage/check_oriented_khovanov.sage")

write_sage_pd_khovanov_directory(
    "data/com_link_gen_10-v0.1.0-com_link_gen-10-3",
    "build/sage_pd_khovanov.txt",
    numeric_only=True,
    workers=8,
    chunksize=1,
    progress_every=25,
)
```

The output file keeps the selected file order, even though computation is
parallel. It is created immediately and flushed as ordered results become
available. If a file fails to parse or Sage fails to compute its Khovanov
polynomial, that file still gets a line and the polynomial field is replaced by a
single-line `ERROR[...]` value. Progress is printed every `progress_every`
completed output lines and includes cache hits, so the progress bar tracks the
final output file rather than only newly computed unique PD codes.

For speed, the exporter computes duplicate `PD_CODE` values only once by
default and writes the result to every matching file line. It also keeps a
persistent success cache at `output_path + ".cache.json"` by default, so reruns
skip previously computed `PD_CODE` values; pass `cache_path=None` to disable
that cache. It schedules uncached jobs by descending crossing count to reduce
parallel tail latency. It uses Python `flush()` for live file visibility but
does not call `fsync()` by default; pass `fsync_every=100` or another positive
line count only if you need periodic disk syncs. Keep `chunksize=1` unless a
benchmark shows your selected range has very uniform cost.

```sage
write_sage_pd_khovanov_directory(
    "data/com_link_gen_10-v0.1.0-com_link_gen-10-3",
    "build/sage_pd_khovanov_0001_1000.txt",
    start_index=1,
    end_index=1000,
    numeric_only=True,
    workers=8,
    chunksize=1,
    progress_every=25,
)
```

## References

- `cppkh`: <https://github.com/GGN-2015/cppkh>
- Spherogram link diagrams: <https://github.com/3-manifolds/Spherogram>
- SageMath knot/link functionality and PD-code conventions:
  <https://doc.sagemath.org/html/en/reference/knots/>
- Knot Atlas planar diagrams: <https://katlas.org/wiki/Planar_Diagrams>
