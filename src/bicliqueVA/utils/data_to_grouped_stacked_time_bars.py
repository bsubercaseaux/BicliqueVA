#!/usr/bin/env python3
"""
Generate grouped stacked bar PGFPlots from `experimental_data.json`.

For each n-value (x-axis group), this plot draws one bar per method.
Each bar has two stacked segments:
  - encoding time (`time`, or 0 for `original`)
  - solving time (`is_time`)

Usage:
    python3 src/bicliqueVA/utils/data_to_grouped_stacked_time_bars.py experimental_data.json > grouped_stacked_time_bars.tex
    python3 src/bicliqueVA/utils/data_to_grouped_stacked_time_bars.py experimental_data.json --p-edge 0.5 --out-file plot.tex
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
    "original": "original",
    "bva": "BVA",
    "factor": "factor",
    "cpp_bp": "BiVA",
    "bva_over_bp": "BiVA+BVA",
    "factor_over_bp": "BiVA+factor",
    "depth_reencode": "Depth",
}

METHOD_COLORS = {
    "original": "gray",
    "bva": "blue!50!white",
    "factor": "red!70!white",
    "cpp_bp": "green!60!black",
    "bva_over_bp": "teal!70!black",
    "factor_over_bp": "lipicsYellow",
    "depth_reencode": "brown!80!black",
}
DEFAULT_METHOD_COLOR = "black!70"


def load_results(path: Path) -> list[dict[str, Any]]:
    with path.open() as f:
        data = json.load(f)
    return data["experimental_results"]


def is_numeric(x: Any) -> bool:
    return isinstance(x, (int, float))


def choose_rows(results: list[dict[str, Any]], p_edge: float | None) -> tuple[list[dict[str, Any]], float]:
    if len(results) == 0:
        raise ValueError("No experimental results found.")

    p_values = sorted({float(r["p_edge"]) for r in results})
    if len(p_values) == 0:
        raise ValueError("No p_edge values found in results.")

    chosen_p = p_edge
    if chosen_p is None:
        chosen_p = p_values[0]
        if len(p_values) > 1:
            print(
                f"[info] Multiple p_edge values found ({p_values}); defaulting to p={chosen_p}. "
                f"Use --p-edge to select another.",
                file=sys.stderr,
            )

    rows = [r for r in results if abs(float(r["p_edge"]) - float(chosen_p)) < 1e-12]
    rows.sort(key=lambda r: int(r["n_vertices"]))
    if len(rows) == 0:
        raise ValueError(f"No rows found for p_edge={chosen_p}.")
    return rows, float(chosen_p)


def collect_methods(rows: list[dict[str, Any]]) -> list[str]:
    present: set[str] = set()
    for row in rows:
        outcomes = row.get("outcomes", {})
        for k in outcomes.keys():
            if k not in {"is_k", "trials"}:
                present.add(k)

    ordered = [m for m in PREFERRED_ALGO_ORDER if m in present]
    extras = sorted(present - set(ordered))
    return ordered + extras


def extract_method_series(rows: list[dict[str, Any]], method: str) -> tuple[list[tuple[int, float]], list[tuple[int, float]], list[int]]:
    enc_pts: list[tuple[int, float]] = []
    sat_pts: list[tuple[int, float]] = []
    missing_ns: list[int] = []

    for row in rows:
        n = int(row["n_vertices"])
        entry = row["outcomes"].get(method)
        if not isinstance(entry, dict):
            missing_ns.append(n)
            continue

        sat_time = entry.get("is_time")
        if not is_numeric(sat_time):
            missing_ns.append(n)
            continue

        if method == "original":
            enc_time = 0.0
        else:
            enc_time = entry.get("time")
            if not is_numeric(enc_time):
                missing_ns.append(n)
                continue

        enc_pts.append((n, float(enc_time)))
        sat_pts.append((n, float(sat_time)))

    return enc_pts, sat_pts, missing_ns


def method_bar_shifts_pt(n_methods: int, bar_width_pt: float, gap_pt: float) -> list[float]:
    if n_methods <= 1:
        return [0.0]
    step = bar_width_pt + gap_pt
    return [(i - (n_methods - 1) / 2.0) * step for i in range(n_methods)]


def generate_plot_tex(rows: list[dict[str, Any]], methods: list[str], p_edge: float) -> str:
    n_values = [int(r["n_vertices"]) for r in rows]
    trials_vals = sorted({r["outcomes"].get("trials") for r in rows if isinstance(r["outcomes"], dict)})
    k_vals = sorted({r["outcomes"].get("is_k") for r in rows if isinstance(r["outcomes"], dict)})

    title_suffix = f"p={p_edge:g}"
    if len(k_vals) == 1 and k_vals[0] is not None:
        title_suffix += f", k={k_vals[0]}"
    if len(trials_vals) == 1 and trials_vals[0] is not None:
        title_suffix += f", trials={trials_vals[0]}"

    lines: list[str] = []
    lines.append(r"\begin{tikzpicture}")
    lines.append(r"\begin{axis}[")
    lines.append(r"    width=0.99\linewidth,")
    lines.append(r"    height=0.60\linewidth,")
    lines.append(r"    ybar=1pt,")
    lines.append(r"    bar width=3pt,")
    lines.append(r"    xlabel={$n$},")
    lines.append(r"    ymin=0,")
    lines.append(r"    ylabel={Time (s)},")
    lines.append(r"  y filter/.code={\pgfmathparse{#1/1000}\pgfmathresult},")
    lines.append(rf"    symbolic x coords={{{','.join(str(n) for n in n_values)}}},")
    lines.append(r"    xtick=data,")
    lines.append(r"    enlarge x limits=0.10,")
    lines.append(r"    ymajorgrids,")
    lines.append(r"    legend style={font=\small, at={(0.2,0.98)}, anchor=north, legend columns=1},")
    lines.append(r"]")

    series_by_method: dict[str, tuple[list[tuple[int, float]], list[tuple[int, float]], list[int]]] = {}
    plottable_methods: list[str] = []
    for method in methods:
        series = extract_method_series(rows, method)
        series_by_method[method] = series
        if len(series[0]) > 0:
            plottable_methods.append(method)

    if len(plottable_methods) == 0:
        raise ValueError("No plottable methods found (missing numeric time/is_time).")

    shifts = method_bar_shifts_pt(len(plottable_methods), bar_width_pt=4.0, gap_pt=0.8)
    legend_items: list[tuple[str, str]] = []
    for active_method_count, method in enumerate(plottable_methods):
        color = METHOD_COLORS.get(method, DEFAULT_METHOD_COLOR)
        shift_pt = shifts[active_method_count]
        enc_pts, sat_pts, missing_ns = series_by_method[method]

        sat_by_n = {n: v for (n, v) in sat_pts}
        enc_by_n = {n: v for (n, v) in enc_pts}
        total_coords = []
        enc_coords = []
        for n in n_values:
            if n in enc_by_n and n in sat_by_n:
                enc = enc_by_n[n]
                sat = sat_by_n[n]
                total_coords.append(f"({n},{enc + sat:.3f})")
                enc_coords.append(f"({n},{enc:.3f})")

        lines.append(
            rf"\addplot+[ybar,  bar shift={shift_pt:.2f}pt, draw={color}, fill={color}!85] "
            rf"coordinates {{{' '.join(total_coords)}}};"
        )
        lines.append(rf"\addlegendentry{{{ALGO_LABELS.get(method, method)} (solving)}}")
        lines.append(
            rf"\addplot+[ybar, bar shift={shift_pt:.2f}pt, draw={color}, fill={color}!45] "
            rf"coordinates {{{' '.join(enc_coords)}}};"
        )
        lines.append(rf"\addlegendentry{{{ALGO_LABELS.get(method, method)} (reencoding)}}")
        legend_items.append((ALGO_LABELS.get(method, method), color))
        if len(missing_ns) > 0:
            ns_text = ",".join(str(x) for x in sorted(set(missing_ns)))
            print(f"[warn] Missing data for method '{method}' at n={ns_text}", file=sys.stderr)

    # for label, color in legend_items:
    #     lines.append(rf"\addlegendimage{{area legend, draw={color}, fill={color}!70}}")
    #     lines.append(rf"\addlegendentry{{{label}}}")

    lines.append(r"\end{axis}")
    lines.append(r"\end{tikzpicture}")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate grouped stacked (encoding+solving) time bars in PGFPlots.")
    parser.add_argument("input_json", type=Path, help="Path to experimental_data.json")
    parser.add_argument("--p-edge", type=float, default=None, help="Select p value (default: smallest available).")
    parser.add_argument("--out-file", type=Path, default=None, help="Write output to file instead of stdout.")
    args = parser.parse_args()

    results = load_results(args.input_json)
    rows, p_edge = choose_rows(results, args.p_edge)
    methods = collect_methods(rows)
    print(methods)
    tex = generate_plot_tex(rows, methods, p_edge)

    if args.out_file is None:
        print(tex, end="")
    else:
        args.out_file.parent.mkdir(parents=True, exist_ok=True)
        args.out_file.write_text(tex)
        print(f"Wrote {args.out_file}")


if __name__ == "__main__":
    main()
