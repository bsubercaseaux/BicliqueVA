#!/usr/bin/env python3
"""
Create a stacked bar PGFPlots figure from `experimental_data.json`.

Each bar corresponds to one method and has two stacked parts:
  1) encoding time (`time`)
  2) solving time with IS constraint (`is_time`)

Usage:
    python3 src/bicliqueVA/utils/data_to_stacked_time_bar.py experimental_data.json > stacked_time_bar.tex
    python3 src/bicliqueVA/utils/data_to_stacked_time_bar.py experimental_data.json --n-vertices 2000 --p-edge 0.5 --out-file plot.tex
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any


PREFERRED_ALGO_ORDER = [
    "original",
    "bva",
    "factor",
    "cpp_bp",
    "bva_over_bp",
    "factor_over_bp",
    "depth_reencode",
]

ALGO_LABELS = {
    "original": "orig.",
    "bva": "BVA",
    "factor": "factor",
    "cpp_bp": "BiVA",
    "bva_over_bp": "BiVA+BVA",
    "factor_over_bp": "BiVA+factor",
    "depth_reencode": "Depth",
}


def load_results(path: Path) -> list[dict[str, Any]]:
    with path.open() as f:
        data = json.load(f)
    return data["experimental_results"]


def as_float_if_numeric(value: Any) -> float | None:
    if isinstance(value, (int, float)):
        return float(value)
    return None


def choose_row(
    results: list[dict[str, Any]],
    n_vertices: int | None,
    p_edge: float | None,
) -> dict[str, Any]:
    if len(results) == 0:
        raise ValueError("No experimental results found.")

    if n_vertices is None:
        n_vertices = max(int(r["n_vertices"]) for r in results)

    candidates = [r for r in results if int(r["n_vertices"]) == n_vertices]
    if len(candidates) == 0:
        raise ValueError(f"No row found with n_vertices={n_vertices}.")

    if p_edge is None:
        p_edge = sorted(float(r["p_edge"]) for r in candidates)[0]

    filtered = [r for r in candidates if abs(float(r["p_edge"]) - float(p_edge)) < 1e-12]
    if len(filtered) == 0:
        raise ValueError(f"No row found with n_vertices={n_vertices}, p_edge={p_edge}.")
    return filtered[0]


def collect_algorithms(outcomes: dict[str, Any]) -> list[str]:
    present = [k for k in outcomes.keys() if k not in {"is_k", "trials"}]
    present_set = set(present)
    ordered = [a for a in PREFERRED_ALGO_ORDER if a in present_set]
    extras = sorted(present_set - set(ordered))
    return ordered + extras


def extract_times_for_row(row: dict[str, Any]) -> tuple[list[str], list[float], list[float], list[str]]:
    outcomes = row["outcomes"]
    algos = collect_algorithms(outcomes)

    selected_algos: list[str] = []
    encoding_times: list[float] = []
    solving_times: list[float] = []
    skipped: list[str] = []

    for algo in algos:
        entry = outcomes.get(algo)
        if not isinstance(entry, dict):
            skipped.append(f"{algo} (non-dict/timeout)")
            continue

        solve_time = as_float_if_numeric(entry.get("is_time"))
        if solve_time is None:
            skipped.append(f"{algo} (missing/timeout is_time)")
            continue

        if algo == "original":
            encode_time = 0.0
        else:
            encode_time = as_float_if_numeric(entry.get("time"))
            if encode_time is None:
                skipped.append(f"{algo} (missing/timeout time)")
                continue

        selected_algos.append(algo)
        encoding_times.append(encode_time)
        solving_times.append(solve_time)

    if len(selected_algos) == 0:
        raise ValueError("No methods with both encoding and solving times available for this row.")
    return selected_algos, encoding_times, solving_times, skipped


def generate_stacked_bar_tex(
    row: dict[str, Any],
    algos: list[str],
    encoding_times: list[float],
    solving_times: list[float],
) -> str:
    n = int(row["n_vertices"])
    p = float(row["p_edge"])
    outcomes = row["outcomes"]
    k = outcomes.get("is_k")
    trials = outcomes.get("trials")

    xcoords = ",".join(algos)
    xticklabels = ",".join(ALGO_LABELS.get(a, a) for a in algos)
    enc_coords = " ".join(f"({algo},{t:.3f})" for algo, t in zip(algos, encoding_times))
    sat_coords = " ".join(f"({algo},{t:.3f})" for algo, t in zip(algos, solving_times))

    title = f"Encoding + Solving Time (n={n}, p={p:g}"
    if isinstance(k, (int, float)):
        title += f", k={int(k)}"
    if isinstance(trials, int):
        title += f", trials={trials}"
    title += ")"

    lines: list[str] = []
    lines.append(r"\begin{tikzpicture}")
    lines.append(r"\begin{axis}[")
    lines.append(r"    ybar stacked,")
    lines.append(r"    bar width=13pt,")
    lines.append(r"    width=0.98\linewidth,")
    lines.append(r"    height=0.58\linewidth,")
    lines.append(r"    ylabel={Time (ms)},")
    lines.append(rf"    title={{{title}}},")
    lines.append(rf"    symbolic x coords={{{xcoords}}},")
    lines.append(r"    xtick=data,")
    lines.append(rf"    xticklabels={{{xticklabels}}},")
    lines.append(r"    x tick label style={rotate=20, anchor=east},")
    lines.append(r"    legend style={font=\small, at={(0.5,-0.2)}, anchor=north, legend columns=2},")
    lines.append(r"    grid=major,")
    lines.append(r"]")
    lines.append(rf"\addplot+[draw=black, fill=blue!55] coordinates {{{enc_coords}}};")
    lines.append(r"\addlegendentry{Encoding}")
    lines.append(rf"\addplot+[draw=black, fill=orange!75] coordinates {{{sat_coords}}};")
    lines.append(r"\addlegendentry{Solving}")
    lines.append(r"\end{axis}")
    lines.append(r"\end{tikzpicture}")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate stacked (encoding+solving) time bar plot in PGFPlots.")
    parser.add_argument("input_json", type=Path, help="Path to experimental_data.json")
    parser.add_argument("--n-vertices", type=int, default=None, help="Select n value (default: max available).")
    parser.add_argument("--p-edge", type=float, default=None, help="Select p value (default: smallest for selected n).")
    parser.add_argument("--out-file", type=Path, default=None, help="Write output to file instead of stdout.")
    args = parser.parse_args()

    results = load_results(args.input_json)
    row = choose_row(results, n_vertices=args.n_vertices, p_edge=args.p_edge)
    algos, encoding_times, solving_times, skipped = extract_times_for_row(row)
    tex = generate_stacked_bar_tex(row, algos, encoding_times, solving_times)

    if args.out_file is None:
        print(tex, end="")
    else:
        args.out_file.parent.mkdir(parents=True, exist_ok=True)
        args.out_file.write_text(tex)
        print(f"Wrote {args.out_file}")

    if len(skipped) > 0:
        print("Skipped methods: " + ", ".join(skipped), file=sys.stderr)


if __name__ == "__main__":
    main()
