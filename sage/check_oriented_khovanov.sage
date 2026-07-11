"""
SageMath cross-checks for cpp_com_link_gen Khovanov orientation output.

Load this file inside Sage:

    sage: load("sage/check_oriented_khovanov.sage")

Then run, for example:

    sage: check_khovanov_file("data/com_link_gen_10-v0.1.0-com_link_gen-10-3/0000001.txt")
    sage: check_khovanov_directory("data/com_link_gen_10-v0.1.0-com_link_gen-10-3", limit=20)

The checker enumerates all 2^n component orientations.  It builds oriented
Gauss codes explicitly, asks Sage to compute integral Khovanov homology, and
compares the distinct results with the KHOVANOV headers written by the C++
pipeline.
"""

import ast
import json
import multiprocessing
import os
import re
import time
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
    text = text.replace(" ", "")
    if text in ("", "0"):
        return ()

    invariants = []
    for part in re.split(r"x|×|\*", text):
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


def sage_homology_to_canonical(homology):
    """Convert Sage's nested Khovanov dictionary to a canonical tuple."""
    terms = []
    for q_degree, by_t_degree in homology.items():
        for t_degree, module in by_t_degree.items():
            invariants = sage_module_to_invariants(module)
            if invariants:
                terms.append((int(q_degree), int(t_degree), tuple(invariants)))
    return tuple(sorted(terms))


def sage_khovanov_for_orientation(pd_code, mask, implementation="native"):
    """Compute canonical Sage Khovanov homology for one orientation mask."""
    if not pd_code:
        link = Knots().one()
    else:
        link = Link(oriented_gauss_code_for_mask(pd_code, mask))
    homology = link.khovanov_homology(ring=ZZ, implementation=implementation)
    return sage_homology_to_canonical(homology)


def sage_khovanov_all_orientations(pd_code, implementation="native"):
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


def check_khovanov_file(path, implementation="native", verbose=True):
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


def _numeric_txt_files(data_dir):
    paths = []
    for filename in os.listdir(data_dir):
        if not filename.endswith(".txt"):
            continue
        stem = filename[:-4]
        if stem.isdigit():
            paths.append(os.path.join(data_dir, filename))
    return sorted(paths, key=lambda p: int(os.path.basename(p)[:-4]))


def _selected_numeric_txt_files(data_dir, limit=None, start_index=None, end_index=None):
    paths = _numeric_txt_files(data_dir)
    if start_index is not None:
        paths = [p for p in paths if int(os.path.basename(p)[:-4]) >= int(start_index)]
    if end_index is not None:
        paths = [p for p in paths if int(os.path.basename(p)[:-4]) <= int(end_index)]
    if limit is not None:
        paths = paths[:int(limit)]
    return paths


def _compact_failure(item):
    compact = {
        "path": item["path"],
        "component_count": item["component_count"],
        "expected_count": item["expected_count"],
        "sage_distinct_count": item["sage_distinct_count"],
        "bound": item["bound"],
        "missing_from_sage": [
            canonical_homology_to_cppkh_text(x) for x in item["missing_from_sage"]
        ],
        "extra_from_sage": [
            canonical_homology_to_cppkh_text(x) for x in item["extra_from_sage"]
        ],
    }
    if "error" in item:
        compact["error"] = item["error"]
    return compact


def _write_json_report(summary, json_report_path):
    if json_report_path is not None:
        with open(json_report_path, "w", encoding="utf-8") as fp:
            json.dump(summary, fp, ensure_ascii=False, indent=2)


def _summary_from_failures(data_dir, implementation, checked, failures, mode, workers=None, elapsed_seconds=None):
    summary = {
        "data_dir": data_dir,
        "implementation": implementation,
        "mode": mode,
        "workers": workers,
        "checked": checked,
        "failed": len(failures),
        "ok": len(failures) == 0,
        "elapsed_seconds": elapsed_seconds,
        "failures": [_compact_failure(item) for item in failures],
    }
    return summary


def check_khovanov_directory(
    data_dir,
    implementation="native",
    limit=None,
    start_index=None,
    end_index=None,
    progress_every=1,
    stop_on_first_failure=False,
    json_report_path=None,
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
    _write_json_report(summary, json_report_path)

    print(
        "Sage Khovanov directory check: checked={} failed={} elapsed={:.1f}s".format(
            checked, len(failures), elapsed
        )
    )
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
                "expected_count": None,
                "sage_distinct_count": None,
                "bound": None,
                "missing_from_sage": [],
                "extra_from_sage": [],
            },
            "error": repr(exc),
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
    implementation="native",
    limit=None,
    start_index=None,
    end_index=None,
    workers=None,
    chunksize=1,
    start_method=None,
    progress_every=10,
    stop_on_first_failure=False,
    json_report_path=None,
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
                if item["error"] is not None:
                    failure = dict(failure)
                    failure["error"] = item["error"]
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
    _write_json_report(summary, json_report_path)
    print(
        "Sage Khovanov parallel check: checked={} failed={} elapsed={:.1f}s".format(
            checked, len(failures), elapsed
        )
    )
    return summary


print("Loaded Sage Khovanov orientation checker.")
print("Try: check_khovanov_file('data/com_link_gen_10-v0.1.0-com_link_gen-10-3/0000001.txt')")
print("For full data, use: check_khovanov_directory_parallel(..., workers=8)")
