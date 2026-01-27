#!/usr/bin/env python3
"""
Usage:
    python make_table.py experimental_data.json > table.tex
"""

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple


ALGO_ORDER: List[str] = ["original", "bva", "factor", "cpp_bp"]
ALGO_LABELS: Dict[str, str] = {
    "original": "orig.",
    "bva": "BVA",
    "factor": "Factor",
    "cpp_bp": "CPP-BP",
}


def load_results(path: Path) -> List[Dict[str, Any]]:
    with path.open() as f:
        data = json.load(f)
    return data["experimental_results"]


def is_timeout(entry: Any) -> bool:
    """Return True if this outcome is a timeout string."""
    return isinstance(entry, str)


def fmt_int(x: int) -> str:
    """Format integers for LaTeX (no thousands separator by default)."""
    return str(x)


def fmt_time(x: float) -> str:
    """Format time in seconds (1 decimal place)."""
    return f"{x:.1f}"


def generate_table(results: List[Dict[str, Any]]) -> str:
    # Sort by n_vertices, then p_edge
    results_sorted = sorted(results, key=lambda r: (r["n_vertices"], r["p_edge"]))

    lines: List[str] = []

    # Optional: macro for timeout
    lines.append(r"\newcommand{\timeout}{\textsc{TO}}")
    lines.append("")

    lines.append(r"\begin{table}")
    lines.append(r"\centering")
    lines.append(r"\small")
    # n, p, 4x clauses, 3x times
    lines.append(r"\begin{tabular}{rr rrrr rrr}")
    lines.append(r"\toprule")
    lines.append(
        r"$n$ & $p$ & \multicolumn{4}{c}{clauses} & "
        r"\multicolumn{3}{c}{time (s)} \\"
    )
    lines.append(r"\cmidrule(lr){3-6}\cmidrule(lr){7-9}")
    lines.append(
        r" & & orig. & BVA & Factor & CPP-BP & "
        r"BVA & Factor & CPP-BP \\"
    )
    lines.append(r"\midrule")

    # Build table body
    for row in results_sorted:
        n = row["n_vertices"]
        p = row["p_edge"]
        outcomes = row["outcomes"]

        # Original (no time)
        orig_cls = fmt_int(outcomes["original"]["cls"])

        # For each heuristic algo, get cls and time (or timeout)
        row_cls: Dict[str, str] = {}
        row_time: Dict[str, str] = {}

        for algo in ["bva", "factor", "cpp_bp"]:
            entry = outcomes.get(algo, None)
            if entry is None:
                row_cls[algo] = r"--"
                row_time[algo] = r"--"
            elif is_timeout(entry):
                row_cls[algo] = r"\timeout"
                row_time[algo] = r"\timeout"
            else:
                row_cls[algo] = fmt_int(entry["cls"])
                row_time[algo] = fmt_time(entry["time"])

        line = (
            f"{n} & {p:.1f} & "
            f"{orig_cls} & "
            f"{row_cls['bva']} & {row_cls['factor']} & {row_cls['cpp_bp']} & "
            f"{row_time['bva']} & {row_time['factor']} & {row_time['cpp_bp']} \\\\"
        )
        lines.append(line)

    lines.append(r"\bottomrule")
    lines.append(r"\end{tabular}")
    lines.append(r"\caption{Experimental results for CNF reductions.}")
    lines.append(r"\label{tab:cnf-reduction}")
    lines.append(r"\end{table}")

    return "\n".join(lines)


def main(argv: List[str]) -> None:
    if len(argv) != 2:
        print("Usage: python make_table.py experimental_data.json", file=sys.stderr)
        sys.exit(1)

    path = Path(argv[1])
    results = load_results(path)
    table = generate_table(results)
    print(table)


if __name__ == "__main__":
    main(sys.argv)
