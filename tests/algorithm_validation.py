#!/usr/bin/env python3
from __future__ import annotations

import argparse
import ast
import importlib.util
import itertools
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
EXE = ROOT / "build" / ("cpp_com_link_gen.exe" if os.name == "nt" else "cpp_com_link_gen")
DATA_ROOT = ROOT / "data"

NAME_RE = re.compile(r"^(m?)([LK])(\d+)([an])(\d+)$")
PD_HEADER_RE = re.compile(r"^//\s*PD_CODE:\s*(.*)$")
KH_HEADER_RE = re.compile(r"^//\s*KHOVANOV:\s*(.*)$")
METHOD_RE = re.compile(r"^L\[(\d+)\s*,\s*(\d+)\]\s*#\s*L\[(\d+)\s*,\s*(\d+)\]\s*$")

KNOWN_KH_SINGLE = {
    "unknot": (
        "[]",
        "q^-1*t^0*Z[0] + q^1*t^0*Z[0]",
    ),
    "right_trefoil": (
        "[[1,5,2,4],[3,1,4,6],[5,3,6,2]]",
        "q^1*t^0*Z[0] + q^3*t^0*Z[0] + q^5*t^2*Z[0] + q^7*t^3*Z[2] + q^9*t^3*Z[0]",
    ),
    "hopf_link": (
        "[[2,3,1,4],[4,1,3,2]]",
        "q^-6*t^-2*Z[0] + q^-4*t^-2*Z[0] + q^-2*t^0*Z[0] + q^0*t^0*Z[0]",
    ),
}

KNOWN_KH_ALL_ORIENTATIONS = {
    "hopf_link": (
        "[[2,3,1,4],[4,1,3,2]]",
        {
            "q^-6*t^-2*Z[0] + q^-4*t^-2*Z[0] + q^-2*t^0*Z[0] + q^0*t^0*Z[0]",
            "q^0*t^0*Z[0] + q^2*t^0*Z[0] + q^4*t^2*Z[0] + q^6*t^2*Z[0]",
        },
    ),
}


def run(args: list[str], *, timeout: int = 120, cwd: Path = ROOT) -> str:
    result = subprocess.run(
        [str(EXE), *args],
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=timeout,
        check=True,
    )
    return result.stdout


def ensure_built() -> None:
    if not EXE.is_file():
        subprocess.check_call([sys.executable, str(ROOT / "build.py")], cwd=ROOT)


def parse_pd(text: str) -> list[list[int]]:
    value = ast.literal_eval(text.strip())
    if not isinstance(value, list):
        raise AssertionError(f"PD is not a list: {text!r}")
    pd: list[list[int]] = []
    for crossing in value:
        if (
            not isinstance(crossing, list)
            or len(crossing) != 4
            or not all(isinstance(x, int) for x in crossing)
        ):
            raise AssertionError(f"bad crossing in PD: {crossing!r}")
        pd.append(crossing)
    return pd


def add_multiedge(adj: dict[int, list[int]], a: int, b: int) -> None:
    if a == b:
        adj[a].extend([b, b])
    else:
        adj[a].append(b)
        adj[b].append(a)


def component_count(pd: list[list[int]], *, require_contiguous: bool = False) -> int:
    if not pd:
        return 0

    labels = [label for crossing in pd for label in crossing]
    if any(label <= 0 for label in labels):
        raise AssertionError(f"PD labels must be positive: {pd!r}")

    counts = Counter(labels)
    bad_counts = {label: count for label, count in counts.items() if count != 2}
    if bad_counts:
        raise AssertionError(f"each PD label must occur exactly twice: {bad_counts!r}")

    if require_contiguous:
        expected = set(range(1, 2 * len(pd) + 1))
        actual = set(labels)
        if actual != expected:
            raise AssertionError(f"generated PD labels are not contiguous: expected {expected}, got {actual}")

    adj: dict[int, list[int]] = defaultdict(list)
    for a, b, c, d in pd:
        add_multiedge(adj, a, c)
        add_multiedge(adj, b, d)

    bad_degrees = {label: len(adj[label]) for label in counts if len(adj[label]) != 2}
    if bad_degrees:
        raise AssertionError(f"PD component graph must be 2-regular: {bad_degrees!r}")

    seen: set[int] = set()
    components = 0
    for start in sorted(counts):
        if start in seen:
            continue
        components += 1
        stack = [start]
        seen.add(start)
        while stack:
            now = stack.pop()
            for nxt in adj[now]:
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
    return components


def orientation_bound(component_total: int) -> int:
    return 1 if component_total <= 1 else 1 << (component_total - 1)


def parse_link_name(name: str) -> tuple[int, bool, bool, int, bool]:
    match = NAME_RE.match(name)
    if not match:
        raise AssertionError(f"invalid link name: {name!r}")
    mirror = match.group(1) == "m"
    is_link = match.group(2) == "L"
    crossing = int(match.group(3))
    non_alternating = match.group(4) == "n"
    index = int(match.group(5))
    return crossing, is_link, non_alternating, index, mirror


def parse_name_list(line: str) -> list[str]:
    text = line.strip()
    if not (text.startswith("[") and text.endswith("]")):
        raise AssertionError(f"invalid link set line: {line!r}")
    body = text[1:-1].strip()
    if not body:
        return []
    return [item.strip() for item in body.split(",") if item.strip()]


def parse_link_rep(text: str) -> tuple[dict[str, list[list[int]]], list[str], list[tuple[tuple[int, int], tuple[int, int]]]]:
    definitions: dict[str, list[list[int]]] = {}
    link_set: list[str] | None = None
    methods: list[tuple[tuple[int, int], tuple[int, int]]] = []

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("//"):
            continue
        if ":" in line:
            name, pd_text = line.split(":", 1)
            definitions[name.strip()] = parse_pd(pd_text)
        elif "#" in line:
            match = METHOD_RE.match(line)
            if not match:
                raise AssertionError(f"invalid method line: {line!r}")
            a_link, a_comp, b_link, b_comp = map(int, match.groups())
            methods.append(((a_link, a_comp), (b_link, b_comp)))
        elif line.startswith("["):
            if link_set is not None:
                raise AssertionError("link representation contains more than one link set")
            link_set = parse_name_list(line)
        else:
            raise AssertionError(f"unrecognized link representation line: {line!r}")

    if link_set is None:
        raise AssertionError("link representation is missing its link set")
    return definitions, link_set, methods


def parse_generated_file(path: Path) -> tuple[list[list[int]], list[str], str]:
    pd: list[list[int]] | None = None
    kh: list[str] = []
    text = path.read_text(encoding="utf-8")
    for line in text.splitlines():
        pd_match = PD_HEADER_RE.match(line)
        if pd_match:
            pd = parse_pd(pd_match.group(1))
            continue
        kh_match = KH_HEADER_RE.match(line)
        if kh_match:
            kh.append(kh_match.group(1).strip())
    if pd is None:
        raise AssertionError(f"{path} has no PD_CODE header")
    return pd, kh, text


def load_amphicheiral(data_root: Path) -> set[str]:
    out: set[str] = set()
    for path in (data_root / "prime_link_knot_10" / "amphicheiral").glob("*.txt"):
        for line in path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line:
                out.add(line)
    return out


def load_prime_pd_table(data_root: Path) -> dict[str, list[list[int]]]:
    amphicheiral = load_amphicheiral(data_root)
    table: dict[str, list[list[int]]] = {}
    for path in (data_root / "prime_link_knot_10" / "pd_code").glob("*.txt"):
        for line in path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            name, pd_text = line.split(":", 1)
            name = name.strip()
            if name.startswith("m") and name[1:] in amphicheiral:
                continue
            pd = parse_pd(pd_text)
            component_count(pd)
            table[name] = pd
    return table


def all_prime_combinations(names: list[str], total_crossing: int) -> list[list[str]]:
    combinations: list[list[str]] = []
    current: list[str] = []

    def dfs(start: int, current_crossing: int) -> None:
        if current_crossing > total_crossing:
            return
        if current:
            combinations.append(list(current))
        if current_crossing >= total_crossing:
            return
        for i in range(start, len(names)):
            crossing = parse_link_name(names[i])[0]
            if current_crossing + crossing > total_crossing:
                continue
            current.append(names[i])
            dfs(i, current_crossing + crossing)
            current.pop()

    dfs(0, 0)
    return combinations


def edge_set_connects(edges: Iterable[tuple[tuple[int, int], tuple[int, int]]], groups: int) -> bool:
    if groups == 1:
        return True
    graph: dict[int, list[int]] = {i: [] for i in range(1, groups + 1)}
    for (a, _), (b, _) in edges:
        if a == b:
            continue
        graph[a].append(b)
        graph[b].append(a)
    seen = {1}
    stack = [1]
    while stack:
        now = stack.pop()
        for nxt in graph[now]:
            if nxt not in seen:
                seen.add(nxt)
                stack.append(nxt)
    return len(seen) == groups


def count_edge_sets(component_counts: list[int]) -> int:
    groups = len(component_counts)
    choose_count = groups - 1
    if choose_count == 0:
        return 1
    nodes = [
        (group_index + 1, component_index + 1)
        for group_index, count in enumerate(component_counts)
        for component_index in range(count)
    ]
    possible_edges = [
        (a, b)
        for a, b in itertools.combinations(nodes, 2)
        if a[0] != b[0]
    ]
    return sum(
        1
        for edges in itertools.combinations(possible_edges, choose_count)
        if edge_set_connects(edges, groups)
    )


def expected_generation_count(data_root: Path, total_crossing: int, max_prime_count: int) -> int:
    table = load_prime_pd_table(data_root)
    names = sorted(table, key=parse_link_name)
    total = 0
    for solution in all_prime_combinations(names, total_crossing):
        if len(solution) > max_prime_count:
            continue
        counts = [component_count(table[name]) for name in solution]
        total += count_edge_sets(counts)
    return total


def assert_generated_file_invariants(path: Path) -> tuple[int, int, int]:
    pd, kh, text = parse_generated_file(path)
    final_components = component_count(pd, require_contiguous=True)
    definitions, link_set, methods = parse_link_rep(text)

    if len(methods) != max(0, len(link_set) - 1):
        raise AssertionError(f"{path.name}: generated method count is not a tree")

    source_crossings = 0
    source_components = 0
    component_counts: list[int] = []
    for name in link_set:
        if name not in definitions:
            raise AssertionError(f"{path.name}: undefined link name {name!r}")
        prime_pd = definitions[name]
        source_crossings += len(prime_pd)
        count = component_count(prime_pd)
        component_counts.append(count)
        source_components += count

    if len(pd) != source_crossings:
        raise AssertionError(f"{path.name}: crossing count changed during connected sum")
    if final_components != source_components - len(methods):
        raise AssertionError(
            f"{path.name}: component count mismatch after connected sums "
            f"({final_components} != {source_components} - {len(methods)})"
        )

    for (a_link, a_comp), (b_link, b_comp) in methods:
        if not (1 <= a_link <= len(link_set) and 1 <= b_link <= len(link_set)):
            raise AssertionError(f"{path.name}: method references an invalid link index")
        if a_link == b_link:
            raise AssertionError(f"{path.name}: method connects a link to itself")
        if not (1 <= a_comp <= component_counts[a_link - 1]):
            raise AssertionError(f"{path.name}: method references an invalid component")
        if not (1 <= b_comp <= component_counts[b_link - 1]):
            raise AssertionError(f"{path.name}: method references an invalid component")

    if methods and not edge_set_connects(methods, len(link_set)):
        raise AssertionError(f"{path.name}: method graph is disconnected")

    if kh:
        if len(kh) != len(set(kh)):
            raise AssertionError(f"{path.name}: duplicate KHOVANOV headers were not deduplicated")
        bound = orientation_bound(final_components)
        if len(kh) > bound:
            raise AssertionError(f"{path.name}: {len(kh)} KH values exceeds orientation bound {bound}")

    return len(pd), final_components, len(kh)


def assert_known_khovanov_values() -> None:
    for name, (pd_text, expected) in KNOWN_KH_SINGLE.items():
        actual = run(["kh", "--pd", pd_text], timeout=120).splitlines()
        if actual != [expected]:
            raise AssertionError(
                f"{name}: unexpected single-PD Khovanov output\n"
                f"expected={[expected]}\nactual={actual}"
            )
    for name, (pd_text, expected) in KNOWN_KH_ALL_ORIENTATIONS.items():
        actual = set(run(["kh-all-orientations", "--pd", pd_text], timeout=120).splitlines())
        if actual != expected:
            raise AssertionError(
                f"{name}: unexpected all-orientations output\n"
                f"expected={expected}\nactual={actual}"
            )
    print("known Khovanov examples passed")


def assert_generation_counts() -> None:
    expected = {
        (3, 1): 4,
        (4, 2): 19,
        (5, 2): 33,
    }
    for params, fixed_count in expected.items():
        actual = expected_generation_count(DATA_ROOT, *params)
        if actual != fixed_count:
            raise AssertionError(f"independent generation count for {params} changed: {actual} != {fixed_count}")
    print("independent generation counts passed")


def run_generated_dataset_validation(
    total_crossing: int,
    max_prime_count: int,
    jobs: int,
) -> list[tuple[list[list[int]], int]]:
    with tempfile.TemporaryDirectory(prefix="cpp_com_link_gen_validation_") as tmp:
        tmp_path = Path(tmp)
        temp_data = tmp_path / "data"
        shutil.copytree(DATA_ROOT / "prime_link_knot_10", temp_data / "prime_link_knot_10")

        expected_count = expected_generation_count(temp_data, total_crossing, max_prime_count)
        generated_stdout = run(
            [
                "generate",
                "--total-crs",
                str(total_crossing),
                "--max-prime-cnt",
                str(max_prime_count),
                "--jobs",
                str(jobs),
                "--data-root",
                str(temp_data),
            ],
            timeout=180,
        )
        generated_dir = Path(generated_stdout.strip().splitlines()[-1].strip('"'))
        files = sorted(generated_dir.glob("*.txt"))
        if len(files) != expected_count:
            raise AssertionError(f"generated {len(files)} files, expected {expected_count}")

        for i, path in enumerate(files, start=1):
            if path.name != f"{i:07d}.txt":
                raise AssertionError(f"non-contiguous generated filename: {path.name}")
            assert_generated_file_invariants(path)

        run(["khovanov", "--dir", str(generated_dir), "--jobs", str(jobs)], timeout=240)

        kh_files = 0
        max_kh = 0
        max_components = 0
        samples: list[tuple[list[list[int]], int]] = []
        for path in files:
            pd, _, _ = parse_generated_file(path)
            _, components, kh_count = assert_generated_file_invariants(path)
            samples.append((pd, components))
            kh_files += int(kh_count > 0)
            max_kh = max(max_kh, kh_count)
            max_components = max(max_components, components)
        if kh_files != len(files):
            raise AssertionError(f"only {kh_files}/{len(files)} files received KHOVANOV headers")

        print(
            f"generated dataset validation passed: {len(files)} files, "
            f"max_components={max_components}, max_distinct_kh={max_kh}"
        )
        return samples


def run_optional_spherogram_check(
    samples: list[tuple[list[list[int]], int]],
    require_spherogram: bool,
) -> None:
    if importlib.util.find_spec("spherogram") is None:
        message = "Spherogram not found; optional mature-link-library cross-check skipped"
        if require_spherogram:
            raise AssertionError(message)
        print(message)
        return

    from spherogram import Link

    checked = 0
    for pd, expected_components in samples:
        if not pd:
            continue
        link = Link(pd, check_planarity=True)
        actual_components = len(link.link_components)
        if actual_components != expected_components:
            raise AssertionError(
                f"Spherogram component mismatch: {actual_components} != {expected_components} for {pd!r}"
            )
        if len(link.crossings) != len(pd):
            raise AssertionError(f"Spherogram crossing-count mismatch for {pd!r}")

        round_trip_pd = [[label + 1 for label in crossing] for crossing in link.PD_code()]
        if len(round_trip_pd) != len(pd):
            raise AssertionError(f"Spherogram PD round-trip crossing-count mismatch for {pd!r}")
        round_trip_components = component_count(round_trip_pd)
        if round_trip_components != expected_components:
            raise AssertionError(
                "Spherogram PD round-trip component mismatch: "
                f"{round_trip_components} != {expected_components} for {pd!r}"
            )
        checked += 1

    print(f"Spherogram cross-check passed for {checked} PD codes")


def main() -> int:
    parser = argparse.ArgumentParser(description="Run stronger algorithm-level validation checks.")
    parser.add_argument("--total-crs", type=int, default=5)
    parser.add_argument("--max-prime-cnt", type=int, default=2)
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--require-spherogram", action="store_true")
    args = parser.parse_args()

    ensure_built()
    assert_known_khovanov_values()
    assert_generation_counts()

    samples = run_generated_dataset_validation(args.total_crs, args.max_prime_cnt, args.jobs)
    for pd_text, _ in KNOWN_KH_SINGLE.values():
        pd = parse_pd(pd_text)
        samples.append((pd, component_count(pd)))
    run_optional_spherogram_check(samples, args.require_spherogram)

    print("algorithm validation passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
