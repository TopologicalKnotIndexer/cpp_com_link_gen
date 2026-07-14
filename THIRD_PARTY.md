# Third-Party Components

This project vendors or consumes the following MIT-licensed components and
data sources.

- `TopologicalKnotIndexer/cppkh`, from
  <https://github.com/TopologicalKnotIndexer/cppkh>, is used as the
  Khovanov homology backend. The vendored source is synchronized to upstream
  commit `37b3cc345b3b80844239708edce2848062dfcde1`. It carries a local C API
  extension for evaluating multiple explicit crossing-sign variants and a fix
  that makes relative `build.py --out` paths resolve from the caller's working
  directory.
- `pd-code-to-diagram` 0.1.8, from <https://pypi.org/project/pd-code-to-diagram/>,
  is vendored under `third_party/pd_code_to_diagram` and used to lay out PD_CODE
  diagrams before SVG rendering.
- The prime knot/link PD tables under `data/prime_link_knot_10` come from
  [`TopologicalKnotIndexer/pd_code_for_prime_link_and_knot`](https://github.com/TopologicalKnotIndexer/pd_code_for_prime_link_and_knot).
  The tables were introduced in commit `4dfa6ce` and are tracked directly so a
  fresh checkout can generate datasets without downloading another project.
  The source project's MIT license is copied to
  `data/prime_link_knot_10/LICENSE`.
- The original Python packages used for behavior comparison are
  `com-link-gen-10`, `link-rep`, `link-rep-to-pd-code`,
  `pd-code-connected-sum`, `pd-code-components`, `pd-code-sanity`,
  `pd-code-pre-nxt`, `pd-code-reverse-component`, and
  `group-diagram-combination`.

The project license remains MIT; preserve the license notices of vendored
third-party sources when redistributing binaries or source archives.

All required dependency sources are committed as ordinary files. The build
does not use Git submodules and does not clone a missing dependency.
