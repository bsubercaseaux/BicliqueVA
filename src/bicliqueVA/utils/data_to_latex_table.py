#!/usr/bin/env python3
"""
Convert `experimental_data.json` to a LaTeX table.

Usage:
    python3 src/bicliqueVA/utils/data_to_latex_table.py experimental_data.json > table.tex
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any


PREFERRED_ALGO_ORDER = [
    "bva",
    "factor",
    "cpp_bp",
    # "bva_over_bp",
    "factor_over_bp",
    # "depth_reencode",
]

ALGO_LABELS = {
    "original": "orig.",
    "bva": "BVA",
    "factor": "factor",
    "cpp_bp": "BiVA",
    # "bva_over_bp": "BiVA+BVA",
    "factor_over_bp": "BiVA+factor",
    # "depth_reencode": "Depth",
}


def is_timeout(entry: Any) -> bool:
    return isinstance(entry, str)


def load_results(path: Path) -> list[dict[str, Any]]:
    with path.open() as f:
        data = json.load(f)
    return data["experimental_results"]


def collect_algorithms(results: list[dict[str, Any]]) -> list[str]:
    present: set[str] = set()
    for row in results:
        present.update(row.get("outcomes", {}).keys())
    present.discard("original")

    ordered = [algo for algo in PREFERRED_ALGO_ORDER if algo in present]
    extras = sorted(present - set(ordered))
    return ordered
    # return ordered + extras


def fmt_cls(entry: Any) -> str:
    if entry is None:
        return "--"
    if is_timeout(entry):
        return r"\timeout"
    return str(entry["cls"])


def fmt_time_ms(entry: Any) -> str:
    if entry is None:
        return "--"
    if is_timeout(entry):
        return r"\timeout"
    return f"{entry['time']:.1f}"


def fmt_aux_vars(entry: Any, original_vars: int) -> str:
    if entry is None:
        return "--"
    if is_timeout(entry):
        return r"\timeout"
    return str(max(0, int(entry["vrs"]) - original_vars))


def generate_table(results: list[dict[str, Any]]) -> str:
    algos = collect_algorithms(results)
    clause_algos = ["original"] + algos
    aux_algos = algos
    time_algos = algos

    n_clause = len(clause_algos)
    n_aux = len(aux_algos)
    n_time = len(time_algos)
    colspec = "r" + "r" * (n_clause + n_aux + n_time)

    start_clause = 2
    end_clause = 1 + n_clause
    start_aux = end_clause + 1
    end_aux = end_clause + n_aux
    start_time = end_aux + 1
    end_time = end_aux + n_time

    lines: list[str] = []
    lines.append(r"\newcommand{\timeout}{\textsc{TO}}")
    lines.append("")
    lines.append(r"\begin{table}")
    lines.append(r"\centering")
    lines.append(r"\small")
    lines.append(rf"\begin{{tabular}}{{{colspec}}}")
    lines.append(r"\toprule")
    lines.append(
        rf"$n$ & \multicolumn{{{n_clause}}}{{c}}{{clauses}} & "
        rf"\multicolumn{{{n_aux}}}{{c}}{{aux vars}} & "
        rf"\multicolumn{{{n_time}}}{{c}}{{time (ms)}} \\"
    )
    lines.append(
        rf"\cmidrule(lr){{{start_clause}-{end_clause}}}"
        rf"\cmidrule(lr){{{start_aux}-{end_aux}}}"
        rf"\cmidrule(lr){{{start_time}-{end_time}}}"
    )

    clause_headers = " & ".join(ALGO_LABELS.get(a, a) for a in clause_algos)
    aux_headers = " & ".join(ALGO_LABELS.get(a, a) for a in aux_algos)
    time_headers = " & ".join(ALGO_LABELS.get(a, a) for a in time_algos)
    lines.append(f" & {clause_headers} & {aux_headers} & {time_headers} \\\\")
    lines.append(r"\midrule")

    for row in sorted(results, key=lambda r: (r["n_vertices"], r["p_edge"])):
        outcomes = row["outcomes"]
        original_vars = int(outcomes["original"]["vrs"])
        cls_cells = [fmt_cls(outcomes.get(a)) for a in clause_algos]
        aux_cells = [fmt_aux_vars(outcomes.get(a), original_vars) for a in aux_algos]
        time_cells = [fmt_time_ms(outcomes.get(a)) for a in time_algos]
        lines.append(
            f"{row['n_vertices']} & "
            f"{' & '.join(cls_cells + aux_cells + time_cells)} \\\\"
        )

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\caption{Experimental results by algorithm (clauses, auxiliary variables, and runtime).}")
    lines.append(r"\label{tab:experimental-results}")
    lines.append(r"\end{table}")
    return "\n".join(lines)


def main(argv: list[str]) -> None:
    if len(argv) != 2:
        print(
            "Usage: python3 src/bicliqueVA/utils/data_to_latex_table.py experimental_data.json",
            file=sys.stderr,
        )
        raise SystemExit(1)

    results = load_results(Path(argv[1]))
    print(generate_table(results))


if __name__ == "__main__":
    main(sys.argv)
