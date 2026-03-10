#!/usr/bin/env python3
"""
Export PGFPlots/TikZ plots from `experimental_data.json`.

This script writes three LaTeX files:
  - clauses_plot.tex
  - time_plot.tex
  - aux_vars_plot.tex

Usage:
    python3 src/bicliqueVA/utils/data_to_latex_plots.py experimental_data.json --out-dir latex_plots
    python3 src/bicliqueVA/utils/data_to_latex_plots.py experimental_data.json --p-edge 0.5
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


PREFERRED_ALGO_ORDER = [
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

MARKERS = [
    "*",
    "square*",
    "triangle*",
    "diamond*",
    "pentagon*",
    "otimes*",
    "oplus*",
    "star",
    "x",
    "+",
]


def load_results(path: Path) -> list[dict[str, Any]]:
    with path.open() as f:
        data = json.load(f)
    return data["experimental_results"]


def is_timeout(entry: Any) -> bool:
    return isinstance(entry, str)


def collect_algorithms(results: list[dict[str, Any]]) -> list[str]:
    present: set[str] = set()
    for row in results:
        present.update(row.get("outcomes", {}).keys())
    present.discard("original")

    ordered = [algo for algo in PREFERRED_ALGO_ORDER if algo in present]
    extras = sorted(present - set(ordered))
    return ordered + extras


def filter_rows_by_p(results: list[dict[str, Any]], p_edge: float | None) -> tuple[list[dict[str, Any]], float | None]:
    p_values = sorted({float(row["p_edge"]) for row in results})
    if len(p_values) == 0:
        return [], None

    chosen_p = p_edge
    if chosen_p is None:
        if len(p_values) == 1:
            chosen_p = p_values[0]
        else:
            chosen_p = p_values[0]
            print(
                f"[info] Multiple p_edge values found ({p_values}); defaulting to p={chosen_p}. "
                f"Use --p-edge to select a different one."
            )

    filtered = [row for row in results if abs(float(row["p_edge"]) - float(chosen_p)) < 1e-12]
    filtered.sort(key=lambda r: r["n_vertices"])
    return filtered, chosen_p


def marker_for_algo(algo: str, index: int) -> str:
    return MARKERS[index % len(MARKERS)]


def coordinates_for_metric(
    rows: list[dict[str, Any]],
    algo: str,
    metric: str,
) -> list[tuple[int, float]]:
    pts: list[tuple[int, float]] = []
    for row in rows:
        outcomes = row.get("outcomes", {})
        original = outcomes.get("original")
        entry = outcomes.get(algo)
        if entry is None or is_timeout(entry):
            continue
        n = int(row["n_vertices"])
        if metric == "clauses":
            y = float(entry["cls"])
        elif metric == "time":
            y = float(entry["time"])
        elif metric == "aux":
            if original is None or is_timeout(original):
                continue
            y = float(max(0, int(entry["vrs"]) - int(original["vrs"])))
        else:
            raise ValueError(f"Unknown metric: {metric}")
        pts.append((n, y))
    return pts


def generate_plot_tex(
    rows: list[dict[str, Any]],
    algos: list[str],
    metric: str,
    y_label: str,
    title: str,
    marker_map: dict[str, str],
) -> str:
    def fmt_y(val: float) -> str:
        if metric == "time":
            return f"{val:.1f}"
        return str(int(round(val)))

    lines: list[str] = []
    lines.append(r"\begin{tikzpicture}")
    lines.append(r"\begin{axis}[")
    lines.append(r"    width=0.95\linewidth,")
    lines.append(r"    height=0.55\linewidth,")
    lines.append(r"    xlabel={$n$},")
    lines.append(rf"    ylabel={{{y_label}}},")
    lines.append(rf"    title={{{title}}},")
    lines.append(r"    grid=major,")
    lines.append(r"    legend style={font=\small, at={(0.5,-0.2)}, anchor=north, legend columns=3},")
    lines.append(r"    mark size=2.5pt,")
    lines.append(r"]")

    for algo in algos:
        pts = coordinates_for_metric(rows, algo=algo, metric=metric)
        if len(pts) == 0:
            continue
        marker = marker_map[algo]
        coords = " ".join(f"({x},{fmt_y(y)})" for x, y in pts)
        lines.append(rf"\addplot+[thick, mark={marker}] coordinates {{{coords}}};")
        lines.append(rf"\addlegendentry{{{ALGO_LABELS.get(algo, algo)}}}")

    lines.append(r"\end{axis}")
    lines.append(r"\end{tikzpicture}")
    return "\n".join(lines) + "\n"


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def main() -> None:
    parser = argparse.ArgumentParser(description="Export LaTeX/PGFPlots from experimental_data.json")
    parser.add_argument("input_json", type=Path, help="Path to experimental_data.json")
    parser.add_argument("--out-dir", type=Path, default=Path("."), help="Directory for output .tex plots")
    parser.add_argument("--p-edge", type=float, default=None, help="Filter by a specific p_edge value")
    args = parser.parse_args()

    results = load_results(args.input_json)
    rows, chosen_p = filter_rows_by_p(results, args.p_edge)
    if len(rows) == 0:
        raise SystemExit("No rows found for the selected p_edge.")

    algos = collect_algorithms(rows)
    clause_algos = ["original"] + algos
    metric_algos = algos

    ordered_for_markers = clause_algos
    marker_map = {
        algo: marker_for_algo(algo, i)
        for i, algo in enumerate(ordered_for_markers)
    }

    p_suffix = f" (p={chosen_p:g})" if chosen_p is not None else ""
    clauses_tex = generate_plot_tex(
        rows=rows,
        algos=clause_algos,
        metric="clauses",
        y_label="Clauses",
        title=f"Clauses vs n{p_suffix}",
        marker_map=marker_map,
    )
    time_tex = generate_plot_tex(
        rows=rows,
        algos=metric_algos,
        metric="time",
        y_label="Time (ms)",
        title=f"Time vs n{p_suffix}",
        marker_map=marker_map,
    )
    aux_tex = generate_plot_tex(
        rows=rows,
        algos=metric_algos,
        metric="aux",
        y_label="Auxiliary Variables",
        title=f"Aux Variables vs n{p_suffix}",
        marker_map=marker_map,
    )

    out_dir = args.out_dir
    write_text(out_dir / "clauses_plot.tex", clauses_tex)
    write_text(out_dir / "time_plot.tex", time_tex)
    write_text(out_dir / "aux_vars_plot.tex", aux_tex)

    print(f"Wrote {out_dir / 'clauses_plot.tex'}")
    print(f"Wrote {out_dir / 'time_plot.tex'}")
    print(f"Wrote {out_dir / 'aux_vars_plot.tex'}")


if __name__ == "__main__":
    main()
