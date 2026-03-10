from ..utils.formula_gen import random_2cnf
from ..partition_algs.externals import run_factor_on_formula, run_bva_on_formula, run_sbva_on_formula, run_cpp_bp_on_formula
from ..partition_algs.biclique_partition import run_biclique_partition_on_formula
from ..partition_algs.custom_bva import run_custom_bva_on_formula
from ..partition_algs.depth_reencoder import run_depth_reencode
from ..partition_algs.biclique_finder import run_finder_partition_on_formula
import random
import argparse

from rich.console import Console
from rich.table import Table

def get_cnf_stats(formula_filename):
    with open(formula_filename, 'r') as f:
        for line in f:
            if line.startswith('p cnf'):
                tokens = line.split()
                return int(tokens[2]), int(tokens[3])
    return 0, 0

def main():
    console = Console()
    parser = argparse.ArgumentParser(description='Run experiments on a single formula.')
    parser.add_argument('-f', '--formula', type=str, required=True, help='Path to the formula file')
    args = parser.parse_args()
    formula_filename = args.formula

    init_vars, init_cls = get_cnf_stats(formula_filename)

    console.print(f"[bold blue]Running experiments on:[/bold blue] {formula_filename}")
    console.print(f"[bold white]Initial Stats:[/bold white] {init_vars} variables, {init_cls} clauses\n")

    table = Table(title="Reencoding Results")
    table.add_column("Method", style="bold cyan", no_wrap=True)
    table.add_column("Vars", justify="right", style="cyan")
    table.add_column("Aux Vars (%)", justify="right", style="magenta")
    table.add_column("Clauses", justify="right", style="cyan")
    table.add_column("Reduction (%)", justify="right", style="magenta")
    table.add_column("Time (s)", justify="right", style="green")


    table.add_row(
        "Initial",
        str(init_vars),
        "",
        str(init_cls),
        "",
        ""
    )


    table.add_section()

    def format_perc(final, initial, is_reduction=False):
        if initial == 0: return "0.0%"
        perc = (final - initial) / initial * 100
        sign = "+" if perc >= 0 else ""
        color = "green" if (is_reduction and perc < 0) or (not is_reduction and perc < 0) else "red"
        # For aux vars, positive is "more vars", for reduction, negative is "better"
        if is_reduction:
             color = "green" if perc < 0 else "red"
        else:
             color = "red" if perc > 0 else "green"
        
        return f"[{color}]{sign}{perc:.1f}%[/{color}]"

    with console.status("[bold green]Running BVA...[/bold green]"):
        output_bva = args.formula.replace(".cnf", "_bva.cnf")
        bva_vars, bva_cls, bva_time = run_bva_on_formula(formula_filename, output_bva)
        table.add_row(
            "BVA", 
            str(bva_vars), 
            format_perc(bva_vars, init_vars),
            str(bva_cls), 
            format_perc(bva_cls, init_cls, is_reduction=True),
            f"{bva_time/1000:.2f}"
        )
    
    with console.status("[bold green]Running SBVA...[/bold green]"):
        output_sbva = args.formula.replace(".cnf", "_sbva.cnf")
        sbva_vars, sbva_cls, sbva_time = run_sbva_on_formula(formula_filename, output_sbva)
        table.add_row(
            "SBVA", 
            str(sbva_vars), 
            format_perc(sbva_vars, init_vars),
            str(sbva_cls), 
            format_perc(sbva_cls, init_cls, is_reduction=True),
            f"{sbva_time/1000:.2f}"
        )

    with console.status("[bold green]Running Factor...[/bold green]"):
        output_factor = args.formula.replace(".cnf", "_factor.cnf")
        factor_vars, factor_cls, factor_time = run_factor_on_formula(formula_filename, output_factor)
        table.add_row(
            "factor", 
            str(factor_vars), 
            format_perc(factor_vars, init_vars),
            str(factor_cls), 
            format_perc(factor_cls, init_cls, is_reduction=True),
            f"{factor_time/1000:.2f}"
        )

    # with console.status("[bold green]Running Biclique Partition (Python)...[/bold green]"):
    #     output_bp = args.formula.replace(".cnf", "_bp.cnf")
    #     bp_vars, bp_cls, bp_time = run_biclique_partition_on_formula(formula_filename, output_formula_filename=output_bp)
    #     table.add_row(
    #         "BicliqueVA", 
    #         str(bp_vars), 
    #         format_perc(bp_vars, init_vars),
    #         str(bp_cls), 
    #         format_perc(bp_cls, init_cls, is_reduction=True),
    #         f"{bp_time/1000:.2f}"
    #     )

    with console.status("[bold green]Running Biclique Partition (cpp)...[/bold green]"):
        output_bp = args.formula.replace(".cnf", "_bp.cnf")
        bp_vars, bp_cls, bp_time = run_cpp_bp_on_formula(formula_filename, output_file_path=output_bp)
        table.add_row(
            "BicliqueVA", 
            str(bp_vars), 
            format_perc(bp_vars, init_vars),
            str(bp_cls), 
            format_perc(bp_cls, init_cls, is_reduction=True),
            f"{bp_time/1000:.2f}"
        )

    console.print(table)


    # custom_bva_vars, custom_bva_cls, custom_bva_time = run_custom_bva_on_formula(formula_filename, output_file_path="tmp/custom_bva_out.cnf")
    # print(f"Custom BVA results: vars={custom_bva_vars}, clauses={custom_bva_cls}, time={custom_bva_time:.2f} ms")
    # bp_vars, bp_cls, bp_time = run_biclique_partition_on_formula(formula_filename)
    # print(f"Biclique Partition (python3) results: vars={bp_vars}, clauses={bp_cls}, time={bp_time:.2f} ms")
    
    # finder_vars, finder_cls, finder_time = run_finder_partition_on_formula(formula_filename)
    # print(f"Finder Partition (python3) results: vars={finder_vars}, clauses={finder_cls}, time={finder_time:.2f} ms")
    

    # finder_bp_vars, finder_bp_cls, finder_bp_time = run_finder_partition_on_formula(f"tmp/result_bp.cnf")
    # print(f"Finder on BP results: vars={finder_bp_vars}, clauses={finder_bp_cls}, time={finder_bp_time}")
    # depth reencoding
    # depth_reencode_vars, depth_reencode_cls, depth_reencode_time = run_depth_reencode(formula_filename)
    # print(f"Depth reencoding results: vars={depth_reencode_vars}, clauses={depth_reencode_cls}, time={depth_reencode_time:.2f} ms")

if __name__ == "__main__":
    main()
