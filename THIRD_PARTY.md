# Third-Party Components

This project vendors or consumes the following MIT-licensed components and
data sources.

- `GGN-2015/cppkh`, from <https://github.com/GGN-2015/cppkh>, is used as the
  Khovanov homology backend.
- `pd-code-to-diagram` 0.1.8, from <https://pypi.org/project/pd-code-to-diagram/>,
  is vendored under `third_party/pd_code_to_diagram` and used to lay out PD_CODE
  diagrams before SVG rendering.
- The prime knot/link PD tables under `data/prime_link_knot_10` come from the
  `prime-link-knot-10` Python package version `0.0.5`; its MIT license is
  copied to `data/prime_link_knot_10/LICENSE`.
- The original Python packages used for behavior comparison are
  `com-link-gen-10`, `link-rep`, `link-rep-to-pd-code`,
  `pd-code-connected-sum`, `pd-code-components`, `pd-code-sanity`,
  `pd-code-pre-nxt`, `pd-code-reverse-component`, and
  `group-diagram-combination`.

The project license remains MIT; preserve the license notices of vendored
third-party sources when redistributing binaries or source archives.
