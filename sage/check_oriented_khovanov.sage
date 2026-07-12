"""
SageMath cross-checks for cpp_com_link_gen Khovanov orientation output.

Load this file inside Sage:

    sage: load("sage/check_oriented_khovanov.sage")

Then run, for example:

    sage: write_sage_pd_khovanov_directory(
    ....:     "data/com_link_gen_10-v0.1.0-com_link_gen-10-3",
    ....:     "sage_pd_khovanov.txt",
    ....:     numeric_only=True,
    ....:     workers=8,
    ....: )

The main export helper extracts PD_CODE from each selected txt file, lets Sage
compute one Khovanov polynomial directly from that PD code with variables
``q`` and ``t``, and writes ordered lines of the form ``filename: polynomial``.
The output file is created immediately and flushed as ordered results become
available.  Per-file errors are written as single-line ``ERROR[...]`` values in
the polynomial field.
"""

import ast
import json
import multiprocessing
import os
import re
import time
import traceback
from collections import Counter, defaultdict

try:
    from sage.all import Link, Knots, ZZ
except Exception:
    # When loaded by Sage, these names are already available in the global
    # namespace.  Keeping this fallback lets `load(...)` work naturally.
    pass


PD_HEADER_RE = re.compile(r"^//\s*PD_CODE:\s*(.*)$")
KH_HEADER_RE = re.compile(r"^//\s*KHOVANOV:\s*(.*)$")
CPPKH_TERM_RE = re.compile(r"q\^(-?\d+)\*t\^(-?\d+)\*Z\[([^\]]*)\]")


def _trim_cppkh_spaces(text):
    return text.strip().replace(" ", "")


def parse_generated_khovanov_file(path):
    """Return ``(pd_code, khovanov_headers)`` from a generated txt file."""
    pd_code = None
    khovanov = []
    with open(path, "r", encoding="utf-8", errors="replace") as fp:
        for line in fp:
            pd_match = PD_HEADER_RE.match(line)
            if pd_match:
                pd_code = ast.literal_eval(pd_match.group(1).strip())
                continue

            kh_match = KH_HEADER_RE.match(line)
            if kh_match:
                khovanov.append(kh_match.group(1).strip())

    if pd_code is None:
        raise ValueError("PD_CODE header not found: {}".format(path))
    if not khovanov:
        raise ValueError("KHOVANOV header not found: {}".format(path))
    return pd_code, khovanov


def parse_pd_code_from_file(path):
    """Return only the ``PD_CODE`` header from a generated txt file."""
    with open(path, "r", encoding="utf-8", errors="replace") as fp:
        for line in fp:
            pd_match = PD_HEADER_RE.match(line)
            if pd_match:
                return ast.literal_eval(pd_match.group(1).strip())
    raise ValueError("PD_CODE header not found: {}".format(path))


def parse_cppkh_homology(text):
    """Parse one cppkh homology line into a canonical tuple."""
    terms = []
    for match in CPPKH_TERM_RE.finditer(_trim_cppkh_spaces(text)):
        q_degree = int(match.group(1))
        t_degree = int(match.group(2))
        inv_text = match.group(3).strip()
        invariants = tuple(
            int(item) for item in inv_text.split(",") if item.strip() != ""
        )
        terms.append((q_degree, t_degree, invariants))

    if not terms and _trim_cppkh_spaces(text):
        raise ValueError("could not parse cppkh homology line: {!r}".format(text))
    return tuple(sorted(terms))


def canonical_homology_to_cppkh_text(canonical):
    """Render a canonical homology tuple in cppkh-like notation."""
    chunks = []
    for q_degree, t_degree, invariants in canonical:
        inv = ",".join(str(x) for x in invariants)
        chunks.append("q^{}*t^{}*Z[{}]".format(q_degree, t_degree, inv))
    return " + ".join(chunks)


def _add_component_edge(graph, a, b):
    graph[a].add(b)
    graph[b].add(a)


def components_from_pd(pd_code):
    """Return sorted PD label components using the (a,c), (b,d) strand pairs."""
    if not pd_code:
        return []

    counts = Counter(label for crossing in pd_code for label in crossing)
    bad = {label: count for label, count in counts.items() if count != 2}
    if bad:
        raise ValueError("invalid PD label multiplicities: {}".format(bad))

    graph = defaultdict(set)
    for a, b, c, d in pd_code:
        _add_component_edge(graph, a, c)
        _add_component_edge(graph, b, d)

    seen = set()
    components = []
    for start in sorted(graph):
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        component = []
        while stack:
            now = stack.pop()
            component.append(now)
            for nxt in sorted(graph[now]):
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        components.append(sorted(component))
    return sorted(components)


def label_component_map(components):
    out = {}
    for index, component in enumerate(components):
        for label in component:
            out[label] = index
    return out


def orientation_bound(component_count):
    return 1 if component_count <= 1 else 2 ** (component_count - 1)


def _sage_base_crossing_signs(pd_code):
    """Use Sage's own PD convention to get crossing signs in PD order."""
    link = Link(pd_code)
    return [int(x) for x in link.orientation()]


def oriented_gauss_code_for_mask(pd_code, mask):
    """
    Build Sage oriented Gauss input for one component-orientation mask.

    The PD convention is [incoming under, incoming over, outgoing under,
    outgoing over].  Reversing a component swaps incoming/outgoing labels on
    every strand of that component.  The crossing sign changes exactly when
    the under and over strands belong to different components and one, but not
    both, of those components is reversed.
    """
    components = components_from_pd(pd_code)
    component_of = label_component_map(components)
    base_signs = _sage_base_crossing_signs(pd_code)
    crossing_signs = list(base_signs)
    event_by_incoming_label = {}

    for crossing_index, crossing in enumerate(pd_code):
        a, b, c, d = crossing
        under_component = component_of[a]
        over_component = component_of[b]
        if component_of[c] != under_component or component_of[d] != over_component:
            raise ValueError("PD crossing is inconsistent with component strands: {}".format(crossing))

        reverse_under = ((mask >> under_component) & 1) != 0
        reverse_over = ((mask >> over_component) & 1) != 0
        if under_component != over_component and reverse_under != reverse_over:
            crossing_signs[crossing_index] = -crossing_signs[crossing_index]

        crossing_id = crossing_index + 1
        under_in, under_out = (c, a) if reverse_under else (a, c)
        over_in, over_out = (d, b) if reverse_over else (b, d)

        if under_in in event_by_incoming_label or over_in in event_by_incoming_label:
            raise ValueError("duplicate incoming PD label while building Gauss code")
        event_by_incoming_label[under_in] = (-crossing_id, under_out)
        event_by_incoming_label[over_in] = (crossing_id, over_out)

    gauss_components = []
    for component in components:
        labels = set(component)
        start = min(component)
        now = start
        seen = set()
        gauss = []
        while now not in seen:
            seen.add(now)
            if now not in event_by_incoming_label:
                raise ValueError("missing incoming event for PD label {}".format(now))
            signed_crossing, nxt = event_by_incoming_label[now]
            gauss.append(signed_crossing)
            now = nxt
        if now != start or seen != labels:
            raise ValueError(
                "component traversal did not close correctly: start={}, now={}, seen={}, labels={}".format(
                    start, now, sorted(seen), sorted(labels)
                )
            )
        gauss_components.append(gauss)

    return [gauss_components, crossing_signs]


def _parse_module_repr(text):
    text = text.replace(" ", "").replace("\u00d7", "x")
    if text in ("", "0"):
        return ()

    invariants = []
    for part in re.split(r"x|\*", text):
        if not part:
            continue
        if part == "Z":
            invariants.append(0)
            continue
        match = re.match(r"Z\^(\d+)$", part)
        if match:
            invariants.extend([0] * int(match.group(1)))
            continue
        match = re.match(r"C(\d+)(?:\^(\d+))?$", part)
        if match:
            order = int(match.group(1))
            multiplicity = int(match.group(2) or "1")
            invariants.extend([order] * multiplicity)
            continue
        raise ValueError("cannot parse Sage module representation: {!r}".format(text))
    return tuple(invariants)


def sage_module_to_invariants(module):
    """Convert a Sage abelian group/module to cppkh-style invariant factors."""
    if module == 0:
        return ()

    for method_name in ("invariants", "elementary_divisors"):
        method = getattr(module, method_name, None)
        if method is None:
            continue
        try:
            values = [int(x) for x in method()]
        except Exception:
            continue
        if values:
            return tuple(0 if x == 0 else abs(x) for x in values)

    return _parse_module_repr(str(module))


def _sage_khovanov_homology(link, implementation=None):
    if implementation is None:
        return link.khovanov_homology(ring=ZZ)
    return link.khovanov_homology(ring=ZZ, implementation=implementation)


def sage_khovanov_for_pd(pd_code, implementation=None):
    """Compute canonical Sage Khovanov homology directly from a PD code."""
    link = Link(pd_code)
    homology = _sage_khovanov_homology(link, implementation=implementation)
    return sage_homology_to_canonical(homology)


def kh_poly_string_q_t_ascending(P):
    """Render a Sage Khovanov polynomial sorted by ascending q, then t degree."""
    q, t = P.parent().gens()
    terms = sorted(P.dict().items(), key=lambda kv: (kv[0][0], kv[0][1]))

    pieces = []
    for (i, j), c in terms:
        term = c * (q ** int(i)) * (t ** int(j))
        pieces.append(str(term))

    if not pieces:
        return "0"
    return " + ".join(pieces).replace("+ -", "- ")


def sage_khovanov_polynomial_for_pd(pd_code):
    """Compute Sage's q,t Khovanov polynomial directly from a PD code."""
    link = Link(pd_code)
    polynomial = link.khovanov_polynomial(var1="q", var2="t")
    return kh_poly_string_q_t_ascending(polynomial)


def sage_homology_to_canonical(homology):
    """Convert Sage's nested Khovanov dictionary to a canonical tuple."""
    terms = []
    for q_degree, by_t_degree in homology.items():
        for t_degree, module in by_t_degree.items():
            invariants = sage_module_to_invariants(module)
            if invariants:
                terms.append((int(q_degree), int(t_degree), tuple(invariants)))
    return tuple(sorted(terms))


def sage_khovanov_for_orientation(pd_code, mask, implementation=None):
    """Compute canonical Sage Khovanov homology for one orientation mask."""
    if not pd_code:
        link = Knots().one()
    else:
        link = Link(oriented_gauss_code_for_mask(pd_code, mask))
    homology = _sage_khovanov_homology(link, implementation=implementation)
    return sage_homology_to_canonical(homology)


def sage_khovanov_all_orientations(pd_code, implementation=None):
    """Return ``[(mask, canonical_homology), ...]`` for all 2^n orientations."""
    component_count = len(components_from_pd(pd_code))
    orientation_count = 1 if component_count == 0 else 2 ** component_count
    return [
        (mask, sage_khovanov_for_orientation(pd_code, mask, implementation=implementation))
        for mask in range(orientation_count)
    ]


def expected_homology_set_from_file(path):
    pd_code, khovanov_headers = parse_generated_khovanov_file(path)
    return pd_code, set(parse_cppkh_homology(line) for line in khovanov_headers)


def check_khovanov_file(path, implementation=None, verbose=True):
    """
    Check one generated txt file against Sage for every component orientation.

    Returns a dictionary with fields:
      ``ok``: bool
      ``component_count``: number of link components from the PD code
      ``expected_count``: number of KHOVANOV headers in the file
      ``sage_distinct_count``: number of distinct Sage results across 2^n masks
      ``missing_from_sage`` / ``extra_from_sage``: canonical homology tuples
    """
    pd_code, expected = expected_homology_set_from_file(path)
    components = components_from_pd(pd_code)
    sage_by_mask = sage_khovanov_all_orientations(pd_code, implementation=implementation)
    actual = set(value for _, value in sage_by_mask)
    bound = orientation_bound(len(components))

    result = {
        "path": path,
        "ok": expected == actual and len(actual) <= bound,
        "component_count": len(components),
        "orientation_count": 1 if not components else 2 ** len(components),
        "bound": bound,
        "expected_count": len(expected),
        "sage_distinct_count": len(actual),
        "missing_from_sage": sorted(expected - actual),
        "extra_from_sage": sorted(actual - expected),
        "mask_results": sage_by_mask,
    }

    if verbose:
        status = "OK" if result["ok"] else "FAIL"
        print(
            "{} {} components={} expected={} sage_distinct={} bound={}".format(
                status, path, result["component_count"], result["expected_count"],
                result["sage_distinct_count"], result["bound"]
            )
        )
        if not result["ok"]:
            if result["missing_from_sage"]:
                print("  in file but not Sage:")
                for item in result["missing_from_sage"]:
                    print("    {}".format(canonical_homology_to_cppkh_text(item)))
            if result["extra_from_sage"]:
                print("  in Sage but not file:")
                for item in result["extra_from_sage"]:
                    print("    {}".format(canonical_homology_to_cppkh_text(item)))
    return result


def _homology_set_text(values, limit=None):
    values = sorted(values)
    if limit is not None:
        values = values[:int(limit)]
    return [canonical_homology_to_cppkh_text(value) for value in values]


def check_sage_membership_file(path, mask=0, implementation=None, verbose=True, raise_on_failure=True):
    """
    Compute one Sage Khovanov homology and check membership in file headers.

    This is a deliberately weaker but faster validation than
    ``check_khovanov_file``.  It asks Sage to construct the link directly from
    the file's ``PD_CODE`` and verifies that Sage's result appears among the
    generated ``KHOVANOV`` values.  Any parse error or Sage error is allowed to
    propagate.

    Only ``mask=0`` is supported here.  The older ``oriented_gauss_code_for_mask``
    helper is intentionally not used for this check, because reconstructing a
    Sage oriented Gauss code from a PD code is convention-sensitive and can
    create a different diagram than the original PD code.
    """
    mask = int(mask)
    if mask != 0:
        raise NotImplementedError(
            "check_sage_membership_file currently supports only mask=0; "
            "it validates Sage's direct PD-code result against the generated KHOVANOV set"
        )
    pd_code, expected = expected_homology_set_from_file(path)
    components = components_from_pd(pd_code)
    orientation_count = 1 if not components else 2 ** len(components)

    sage_value = sage_khovanov_for_pd(pd_code, implementation=implementation)
    ok = sage_value in expected
    result = {
        "path": path,
        "ok": ok,
        "mask": mask,
        "component_count": len(components),
        "orientation_count": orientation_count,
        "expected_count": len(expected),
        "sage_homology": sage_value,
        "sage_homology_text": canonical_homology_to_cppkh_text(sage_value),
        "file_homology_text": _homology_set_text(expected),
    }

    if verbose:
        status = "OK" if ok else "FAIL"
        print(
            "{} {} mask={} components={} orientations={} file_kh={}".format(
                status, path, mask, len(components), orientation_count, len(expected)
            )
        )
        if not ok:
            print("  Sage homology:")
            print("    {}".format(result["sage_homology_text"]))
            print("  File KHOVANOV values:")
            for item in result["file_homology_text"]:
                print("    {}".format(item))

    if not ok and raise_on_failure:
        raise AssertionError(
            "{}: Sage Khovanov result for mask={} is not present in KHOVANOV headers\n"
            "Sage: {}\n"
            "File values:\n{}".format(
                path,
                mask,
                result["sage_homology_text"],
                "\n".join("  " + item for item in result["file_homology_text"]),
            )
        )
    return result


def diagnose_khovanov_file(path, mask=0, implementation=None):
    """
    Run one file and one orientation mask without suppressing exceptions.

    This is useful when a parallel run reports many failures quickly.  The
    printed phase markers show whether parsing, PD component traversal, Sage
    Gauss-code construction, or Sage Khovanov computation is failing.
    """
    print("diagnose: parsing {}".format(path))
    pd_code, expected = expected_homology_set_from_file(path)
    print("diagnose: crossings={} expected_headers={}".format(len(pd_code), len(expected)))

    components = components_from_pd(pd_code)
    print(
        "diagnose: components={} orientations={} bound={}".format(
            len(components),
            1 if not components else 2 ** len(components),
            orientation_bound(len(components)),
        )
    )

    print("diagnose: building oriented Gauss code for mask={}".format(int(mask)))
    gauss_code = oriented_gauss_code_for_mask(pd_code, int(mask))
    print(
        "diagnose: gauss_components={} crossing_signs={}".format(
            len(gauss_code[0]), len(gauss_code[1])
        )
    )

    print("diagnose: computing Sage Khovanov homology")
    canonical = sage_khovanov_for_orientation(pd_code, int(mask), implementation=implementation)
    print("diagnose: computed_terms={}".format(len(canonical)))
    return {
        "path": path,
        "component_count": len(components),
        "mask": int(mask),
        "homology": canonical,
        "homology_text": canonical_homology_to_cppkh_text(canonical),
    }


def _numeric_txt_files(data_dir):
    paths = []
    for filename in os.listdir(data_dir):
        if not filename.endswith(".txt"):
            continue
        stem = filename[:-4]
        if stem.isdigit():
            paths.append(os.path.join(data_dir, filename))
    return sorted(paths, key=lambda p: int(os.path.basename(p)[:-4]))


def _txt_sort_key(path):
    stem = os.path.basename(path)[:-4]
    if stem.isdigit():
        return (0, int(stem), stem)
    return (1, stem)


def _txt_files(data_dir, recursive=False):
    paths = []
    if recursive:
        for root, _, filenames in os.walk(data_dir):
            for filename in filenames:
                if filename.endswith(".txt"):
                    paths.append(os.path.join(root, filename))
    else:
        for filename in os.listdir(data_dir):
            path = os.path.join(data_dir, filename)
            if filename.endswith(".txt") and os.path.isfile(path):
                paths.append(path)
    return sorted(paths, key=_txt_sort_key)


def _selected_numeric_txt_files(data_dir, limit=None, start_index=None, end_index=None):
    paths = _numeric_txt_files(data_dir)
    if start_index is not None:
        paths = [p for p in paths if int(os.path.basename(p)[:-4]) >= int(start_index)]
    if end_index is not None:
        paths = [p for p in paths if int(os.path.basename(p)[:-4]) <= int(end_index)]
    if limit is not None:
        paths = paths[:int(limit)]
    return paths


def _selected_txt_files(
    data_dir,
    limit=None,
    start_index=None,
    end_index=None,
    recursive=False,
    numeric_only=False,
):
    paths = _numeric_txt_files(data_dir) if numeric_only else _txt_files(data_dir, recursive=recursive)
    if start_index is not None:
        paths = [
            p for p in paths
            if os.path.basename(p)[:-4].isdigit()
            and int(os.path.basename(p)[:-4]) >= int(start_index)
        ]
    if end_index is not None:
        paths = [
            p for p in paths
            if os.path.basename(p)[:-4].isdigit()
            and int(os.path.basename(p)[:-4]) <= int(end_index)
        ]
    if limit is not None:
        paths = paths[:int(limit)]
    return paths


def check_sage_membership_directory(
    data_dir,
    mask=0,
    implementation=None,
    limit=None,
    start_index=None,
    end_index=None,
    recursive=False,
    numeric_only=False,
    progress_every=25,
    json_report_path=None,
):
    """
    Check every selected txt file by one Sage Khovanov membership test.

    For each file, this extracts ``PD_CODE``, lets Sage build the link directly
    from that PD code, computes one Sage Khovanov homology, and checks that this
    value is present in the file's existing ``KHOVANOV`` headers.  The function
    stops immediately on the first parse error, Sage error, or membership
    failure.  Only ``mask=0`` is supported.
    """
    paths = _selected_txt_files(
        data_dir,
        limit=limit,
        start_index=start_index,
        end_index=end_index,
        recursive=recursive,
        numeric_only=numeric_only,
    )
    if not paths:
        raise ValueError("no txt files selected under {}".format(data_dir))

    print(
        "Sage single-membership check starting: files={} mask={} recursive={} numeric_only={}".format(
            len(paths), int(mask), bool(recursive), bool(numeric_only)
        )
    )

    checked = 0
    started_at = time.time()
    results = []
    for path in paths:
        checked += 1
        verbose = progress_every and (
            checked == 1 or checked % int(progress_every) == 0 or checked == len(paths)
        )
        result = check_sage_membership_file(
            path,
            mask=mask,
            implementation=implementation,
            verbose=verbose,
            raise_on_failure=True,
        )
        results.append(result)

    elapsed = time.time() - started_at
    summary = {
        "data_dir": data_dir,
        "implementation": implementation,
        "mode": "single-membership",
        "mask": int(mask),
        "checked": int(checked),
        "ok": True,
        "elapsed_seconds": float(elapsed),
        "recursive": bool(recursive),
        "numeric_only": bool(numeric_only),
        "results": [
            {
                "path": item["path"],
                "mask": item["mask"],
                "component_count": item["component_count"],
                "orientation_count": item["orientation_count"],
                "expected_count": item["expected_count"],
                "sage_homology": item["sage_homology_text"],
            }
            for item in results
        ],
    }
    print(
        "Sage single-membership check passed: checked={} elapsed={:.1f}s".format(
            checked, elapsed
        )
    )
    _write_json_report(summary, json_report_path)
    return summary


def _sage_pd_khovanov_export_worker(task):
    index, path, label, implementation = task
    del implementation
    try:
        pd_code = parse_pd_code_from_file(path)
        polynomial = sage_khovanov_polynomial_for_pd(pd_code)
        return {
            "ok": True,
            "index": int(index),
            "path": path,
            "label": label,
            "homology": polynomial,
            "error": None,
        }
    except Exception as exc:
        return {
            "ok": False,
            "index": int(index),
            "path": path,
            "label": label,
            "homology": _one_line_error(exc),
            "error": _one_line_error(exc),
            "error_type": exc.__class__.__name__,
            "error_message": str(exc),
            "traceback": traceback.format_exc(),
        }


def _one_line_error(exc):
    text = "{}: {}".format(exc.__class__.__name__, str(exc))
    text = " ".join(text.replace("\r", " ").replace("\n", " ").split())
    if not text:
        text = exc.__class__.__name__
    return "ERROR[{}]".format(text)


def _flush_output_file(fp):
    fp.flush()
    try:
        os.fsync(fp.fileno())
    except OSError:
        pass


def _output_label_for_path(data_dir, path, recursive=False):
    if recursive:
        return os.path.relpath(path, data_dir).replace(os.sep, "/")
    return os.path.basename(path)


def write_sage_pd_khovanov_directory(
    data_dir,
    output_path,
    implementation=None,
    limit=None,
    start_index=None,
    end_index=None,
    recursive=False,
    numeric_only=False,
    workers=None,
    chunksize=1,
    start_method=None,
    progress_every=25,
):
    """
    Compute one Sage Khovanov polynomial per txt file and write ordered lines.

    This function does not inspect existing ``KHOVANOV`` headers and does not
    compare results.  It extracts only ``PD_CODE`` from each selected ``.txt``
    file, computes ``Link(pd_code).khovanov_polynomial(var1="q", var2="t")``
    in parallel, and writes one line per file:

        filename: q^...*t^...

    Polynomial terms are sorted by ascending ``q`` exponent, then ascending
    ``t`` exponent.  The output order is the selected file order, not worker
    completion order.
    The output file is opened immediately and flushed as soon as the next
    ordered result line is available, so it can be watched while Sage runs.
    If parsing or Sage computation fails for a file, the corresponding output
    line is still written, with the polynomial field replaced by a single-line
    ``ERROR[...]`` value.
    """
    data_dir = os.path.abspath(str(data_dir))
    output_path = os.path.abspath(str(output_path))
    paths = _selected_txt_files(
        data_dir,
        limit=limit,
        start_index=start_index,
        end_index=end_index,
        recursive=recursive,
        numeric_only=numeric_only,
    )
    if not paths:
        raise ValueError("no txt files selected under {}".format(data_dir))

    worker_count = _default_worker_count(workers)
    ctx = _multiprocessing_context(start_method)
    tasks = [
        (
            index,
            path,
            _output_label_for_path(data_dir, path, recursive=recursive),
            implementation,
        )
        for index, path in enumerate(paths)
    ]

    print(
        "Sage PD Khovanov export starting: files={} workers={} start_method={} output={}".format(
            len(tasks), worker_count, ctx.get_start_method(), output_path
        )
    )

    started_at = time.time()
    checked = 0
    written = 0
    error_count = 0
    results = [None] * len(tasks)
    pool = ctx.Pool(processes=worker_count)
    pool_closed = False
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    try:
        with open(output_path, "w", encoding="utf-8", newline="\n") as fp:
            _flush_output_file(fp)
            iterator = pool.imap_unordered(
                _sage_pd_khovanov_export_worker,
                tasks,
                chunksize=max(1, int(chunksize)),
            )
            for item in iterator:
                checked += 1
                results[int(item["index"])] = item
                if item.get("error") is not None:
                    error_count += 1

                while written < len(results) and results[written] is not None:
                    ready = results[written]
                    fp.write("{}: {}\n".format(ready["label"], ready["homology"]))
                    written += 1
                _flush_output_file(fp)

                if progress_every and (
                    checked == 1 or checked % int(progress_every) == 0 or checked == len(tasks)
                ):
                    elapsed = time.time() - started_at
                    speed = checked / elapsed if elapsed > 0 else 0.0
                    print(
                        "  computed {}/{} written={} errors={} speed={:.2f}/s".format(
                            checked, len(tasks), written, error_count, speed
                        )
                    )
            else:
                pool.close()
                pool_closed = True
    except Exception:
        if not pool_closed:
            pool.terminate()
            pool_closed = True
        raise
    finally:
        if not pool_closed:
            pool.terminate()
        pool.join()

    missing = [index for index, item in enumerate(results) if item is None]
    if missing:
        raise RuntimeError("missing worker result indices: {}".format(missing[:10]))

    if written != len(results):
        raise RuntimeError("only wrote {}/{} result lines".format(written, len(results)))

    elapsed = time.time() - started_at
    print(
        "Sage PD Khovanov export written: files={} errors={} output={} elapsed={:.1f}s".format(
            len(results), error_count, output_path, elapsed
        )
    )
    return {
        "data_dir": data_dir,
        "output_path": output_path,
        "implementation": implementation,
        "checked": int(len(results)),
        "errors": int(error_count),
        "elapsed_seconds": float(elapsed),
        "workers": int(worker_count),
        "recursive": bool(recursive),
        "numeric_only": bool(numeric_only),
    }


def _canonical_list_to_text(values, limit=None):
    values = list(values or [])
    if limit is not None:
        values = values[:int(limit)]
    return [canonical_homology_to_cppkh_text(value) for value in values]


def _failure_reasons(item):
    reasons = []
    if item.get("error"):
        error_type = item.get("error_type") or "exception"
        error_message = item.get("error_message") or item.get("error")
        reasons.append("{}: {}".format(error_type, error_message))

    expected_count = item.get("expected_count")
    sage_distinct_count = item.get("sage_distinct_count")
    bound = item.get("bound")
    if expected_count is not None and sage_distinct_count is not None:
        if int(expected_count) != int(sage_distinct_count):
            reasons.append(
                "distinct-count mismatch: file={} sage={}".format(
                    int(expected_count), int(sage_distinct_count)
                )
            )

    if bound is not None and sage_distinct_count is not None:
        if int(sage_distinct_count) > int(bound):
            reasons.append(
                "Sage distinct count exceeds 2^(n-1) bound: sage={} bound={}".format(
                    int(sage_distinct_count), int(bound)
                )
            )

    missing_count = len(item.get("missing_from_sage", []))
    extra_count = len(item.get("extra_from_sage", []))
    if missing_count:
        reasons.append(
            "{} KHOVANOV value(s) are in the generated file but not in Sage".format(
                missing_count
            )
        )
    if extra_count:
        reasons.append(
            "{} Sage value(s) are missing from the generated file".format(extra_count)
        )

    if not reasons:
        reasons.append("ok flag is false, but no detailed mismatch was recorded")
    return reasons


def _compact_mask_results(item, limit=None):
    mask_results = item.get("mask_results", [])
    if limit is not None:
        mask_results = mask_results[:int(limit)]
    return [
        {
            "mask": int(mask),
            "homology": canonical_homology_to_cppkh_text(value),
        }
        for mask, value in mask_results
    ]


def _compact_failure(item):
    compact = {
        "path": item["path"],
        "reasons": _failure_reasons(item),
        "component_count": item.get("component_count"),
        "orientation_count": item.get("orientation_count"),
        "expected_count": item.get("expected_count"),
        "sage_distinct_count": item.get("sage_distinct_count"),
        "bound": item.get("bound"),
        "missing_from_sage": _canonical_list_to_text(item.get("missing_from_sage", [])),
        "extra_from_sage": _canonical_list_to_text(item.get("extra_from_sage", [])),
    }
    if item.get("mask_results"):
        compact["sage_mask_results"] = _compact_mask_results(item)
    for key in ("error", "error_type", "error_message", "traceback"):
        if item.get(key):
            compact[key] = item[key]
    return compact


def _json_safe(value):
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, set):
        return [_json_safe(item) for item in sorted(value, key=str)]
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    try:
        return int(value)
    except Exception:
        pass
    try:
        return float(value)
    except Exception:
        pass
    return str(value)


def _write_json_report(summary, json_report_path):
    if json_report_path is not None:
        with open(json_report_path, "w", encoding="utf-8") as fp:
            json.dump(_json_safe(summary), fp, ensure_ascii=False, indent=2)


def _summary_from_failures(data_dir, implementation, checked, failures, mode, workers=None, elapsed_seconds=None):
    summary = {
        "data_dir": data_dir,
        "implementation": implementation,
        "mode": mode,
        "workers": None if workers is None else int(workers),
        "checked": int(checked),
        "failed": int(len(failures)),
        "ok": bool(len(failures) == 0),
        "elapsed_seconds": None if elapsed_seconds is None else float(elapsed_seconds),
        "failures": [_compact_failure(item) for item in failures],
    }
    return summary


def _failure_summary_line(item):
    return "; ".join(_failure_reasons(item)) + " | counts: components={} orientations={} expected={} sage={} bound={}".format(
        item.get("component_count"),
        item.get("orientation_count"),
        item.get("expected_count"),
        item.get("sage_distinct_count"),
        item.get("bound"),
    )


def _print_text_block(label, values, limit):
    if not values:
        return
    print("      {}:".format(label))
    for text in _canonical_list_to_text(values, limit=limit):
        print("        {}".format(text))
    if len(values) > int(limit):
        print("        ... {} more".format(len(values) - int(limit)))


def _print_traceback_block(item, traceback_lines):
    tb = item.get("traceback")
    if not tb:
        return
    lines = tb.strip().splitlines()
    print("      traceback last {} line(s):".format(int(traceback_lines)))
    for line in lines[-int(traceback_lines):]:
        print("        {}".format(line))


def _print_mask_results_block(item, limit):
    mask_results = item.get("mask_results", [])
    if not mask_results:
        return
    print("      first Sage mask result(s):")
    for entry in _compact_mask_results(item, limit=limit):
        print("        mask={}: {}".format(entry["mask"], entry["homology"]))
    if len(mask_results) > int(limit):
        print("        ... {} more mask result(s)".format(len(mask_results) - int(limit)))


def _print_failure_sample(failures, limit=3, homology_limit=2, mask_limit=4, traceback_lines=8):
    if not failures:
        return
    print("First failure samples with reasons:")
    for index, item in enumerate(failures[:int(limit)], start=1):
        print("  [{}] {}".format(index, item.get("path")))
        for reason in _failure_reasons(item):
            print("      reason: {}".format(reason))
        print(
            "      counts: components={} orientations={} expected={} sage={} bound={}".format(
                item.get("component_count"),
                item.get("orientation_count"),
                item.get("expected_count"),
                item.get("sage_distinct_count"),
                item.get("bound"),
            )
        )
        _print_text_block("file-only KHOVANOV values", item.get("missing_from_sage", []), homology_limit)
        _print_text_block("sage-only KHOVANOV values", item.get("extra_from_sage", []), homology_limit)
        _print_mask_results_block(item, mask_limit)
        _print_traceback_block(item, traceback_lines)


def check_khovanov_directory(
    data_dir,
    implementation=None,
    limit=None,
    start_index=None,
    end_index=None,
    progress_every=1,
    stop_on_first_failure=False,
    json_report_path=None,
    failure_sample_limit=3,
    failure_homology_limit=2,
    failure_mask_limit=4,
    failure_traceback_lines=8,
):
    """
    Serial checker for generated numbered txt files in a directory.

    Use ``limit`` for a quick sample, or ``start_index`` / ``end_index`` for a
    numeric filename range.  For the full 10-crossing data set, prefer
    ``check_khovanov_directory_parallel``.
    """
    paths = _selected_numeric_txt_files(data_dir, limit, start_index, end_index)

    failures = []
    checked = 0
    started_at = time.time()
    for path in paths:
        checked += 1
        verbose = progress_every and (checked == 1 or checked % int(progress_every) == 0)
        result = check_khovanov_file(path, implementation=implementation, verbose=verbose)
        if not result["ok"]:
            failures.append(result)
            if stop_on_first_failure:
                break

    elapsed = time.time() - started_at
    summary = _summary_from_failures(
        data_dir, implementation, checked, failures, "serial", workers=1, elapsed_seconds=elapsed
    )

    print(
        "Sage Khovanov directory check: checked={} failed={} elapsed={:.1f}s".format(
            checked, len(failures), elapsed
        )
    )
    _print_failure_sample(
        failures,
        limit=failure_sample_limit,
        homology_limit=failure_homology_limit,
        mask_limit=failure_mask_limit,
        traceback_lines=failure_traceback_lines,
    )
    _write_json_report(summary, json_report_path)
    return summary


def _check_khovanov_worker(task):
    path, implementation = task
    try:
        result = check_khovanov_file(path, implementation=implementation, verbose=False)
        return {
            "path": path,
            "ok": result["ok"],
            "result": result,
            "error": None,
        }
    except Exception as exc:
        return {
            "path": path,
            "ok": False,
            "result": {
                "path": path,
                "component_count": None,
                "orientation_count": None,
                "expected_count": None,
                "sage_distinct_count": None,
                "bound": None,
                "missing_from_sage": [],
                "extra_from_sage": [],
            },
            "error": repr(exc),
            "error_type": exc.__class__.__name__,
            "error_message": str(exc),
            "traceback": traceback.format_exc(),
        }


def _default_worker_count(worker_count):
    if worker_count is not None:
        return max(1, int(worker_count))
    cpu_count = multiprocessing.cpu_count()
    return max(1, cpu_count - 1)


def _multiprocessing_context(start_method):
    if start_method is not None:
        return multiprocessing.get_context(start_method)
    methods = multiprocessing.get_all_start_methods()
    if "fork" in methods:
        return multiprocessing.get_context("fork")
    return multiprocessing.get_context(methods[0])


def check_khovanov_directory_parallel(
    data_dir,
    implementation=None,
    limit=None,
    start_index=None,
    end_index=None,
    workers=None,
    chunksize=1,
    start_method=None,
    progress_every=10,
    stop_on_first_failure=False,
    json_report_path=None,
    failure_sample_limit=3,
    failure_homology_limit=2,
    failure_mask_limit=4,
    failure_traceback_lines=8,
):
    """
    Process-based parallel checker for generated numbered txt files.

    This is the recommended checker for the full data set.  Each process checks
    one file at a time, and within that file enumerates all ``2^n`` component
    orientations.  The default start method prefers ``fork`` when available,
    which is the most reliable mode after loading this script inside Sage.

    Example:

        check_khovanov_directory_parallel(
            "data/com_link_gen_10-v0.1.0-com_link_gen-10-3",
            workers=8,
            progress_every=25,
            json_report_path="sage_khovanov_parallel.json",
        )
    """
    paths = _selected_numeric_txt_files(data_dir, limit, start_index, end_index)
    worker_count = _default_worker_count(workers)
    if not paths:
        summary = _summary_from_failures(
            data_dir, implementation, 0, [], "parallel", workers=worker_count, elapsed_seconds=0.0
        )
        _write_json_report(summary, json_report_path)
        print("Sage Khovanov parallel check: checked=0 failed=0")
        return summary

    ctx = _multiprocessing_context(start_method)
    tasks = [(path, implementation) for path in paths]
    failures = []
    checked = 0
    started_at = time.time()
    print(
        "Sage Khovanov parallel check starting: files={} workers={} start_method={}".format(
            len(tasks), worker_count, ctx.get_start_method()
        )
    )

    pool = ctx.Pool(processes=worker_count)
    pool_closed = False
    try:
        iterator = pool.imap_unordered(_check_khovanov_worker, tasks, chunksize=max(1, int(chunksize)))
        for item in iterator:
            checked += 1
            if not item["ok"]:
                failure = item["result"]
                if item.get("error") is not None:
                    failure = dict(failure)
                    for key in ("error", "error_type", "error_message", "traceback"):
                        if item.get(key):
                            failure[key] = item[key]
                failures.append(failure)
                if stop_on_first_failure:
                    pool.terminate()
                    break

            if progress_every and (checked == 1 or checked % int(progress_every) == 0 or checked == len(tasks)):
                elapsed = time.time() - started_at
                speed = checked / elapsed if elapsed > 0 else 0.0
                remaining = (len(tasks) - checked) / speed if speed > 0 else None
                left_text = "{:.1f}s".format(remaining) if remaining is not None else "unknown"
                print(
                    "  checked {}/{} failed={} speed={:.2f}/s left={}".format(
                        checked, len(tasks), len(failures), speed, left_text
                    )
            )
        else:
            pool.close()
            pool_closed = True
    except Exception:
        pool.terminate()
        pool_closed = True
        raise
    finally:
        if not pool_closed:
            pool.terminate()
        pool.join()

    elapsed = time.time() - started_at
    summary = _summary_from_failures(
        data_dir,
        implementation,
        checked,
        failures,
        "parallel",
        workers=worker_count,
        elapsed_seconds=elapsed,
    )
    print(
        "Sage Khovanov parallel check: checked={} failed={} elapsed={:.1f}s".format(
            checked, len(failures), elapsed
        )
    )
    _print_failure_sample(
        failures,
        limit=failure_sample_limit,
        homology_limit=failure_homology_limit,
        mask_limit=failure_mask_limit,
        traceback_lines=failure_traceback_lines,
    )
    _write_json_report(summary, json_report_path)
    return summary


print("Loaded Sage Khovanov orientation checker.")
print("Export Sage PD Khovanov polynomials: write_sage_pd_khovanov_directory(..., output_path, workers=8)")
