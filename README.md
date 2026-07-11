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

## Smoke Test

```bash
python tests/smoke_test.py
```

The test uses a temporary data root and does not write generated datasets into
the repository.

## References

- `cppkh`: <https://github.com/GGN-2015/cppkh>
- SageMath knot/link functionality and PD-code conventions:
  <https://doc.sagemath.org/html/en/reference/knots/>
- Knot Atlas planar diagrams: <https://katlas.org/wiki/Planar_Diagrams>
