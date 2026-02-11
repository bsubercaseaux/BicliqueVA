from ..utils.formula_gen import random_2cnf
from ..partition_algs.externals import run_factor_on_formula, run_bva_on_formula, run_cpp_bp_on_formula
from ..partition_algs.biclique_partition import run_biclique_partition_on_formula
from ..partition_algs.custom_bva import run_custom_bva_on_formula
from ..partition_algs.depth_reencoder import run_depth_reencode
from ..partition_algs.biclique_finder import run_finder_partition_on_formula
import random
import argparse

def main():
    parser = argparse.ArgumentParser(description='Run experiments on a single formula.')
    parser.add_argument('-f', '--formula', type=str, help='Path to the formula file')
    args = parser.parse_args()
    formula_filename = args.formula

    output_bva = args.formula.replace(".cnf", "_bva.cnf")
    bva_vars, bva_cls, bva_time = run_bva_on_formula(formula_filename, output_bva)
    print(f"BVA results: vars={bva_vars}, clauses={bva_cls}, time={bva_time:.2f}")

    print('*'*50)
    output_factor = args.formula.replace(".cnf", "_factor.cnf")
    factor_vars, factor_cls, factor_time = run_factor_on_formula(formula_filename, output_factor)
    print(f"Factor results: vars={factor_vars}, clauses={factor_cls}, time={factor_time:.2f} ms")
    print('*'*50)
    output_bp = args.formula.replace(".cnf", "_bp.cnf")
    bp_vars, bp_cls, bp_time = run_biclique_partition_on_formula(formula_filename, output_formula_filename=output_bp)
    print(f"Biclique Partition (python3) results: vars={bp_vars}, clauses={bp_cls}, time={bp_time:.2f} ms")
    # for t in [1]: # [1,2,4,8]:
    #     cpp_bp_vars, cpp_bp_cls, cpp_bp_time = run_cpp_bp_on_formula(formula_filename, n_threads=t)
    #     print(f"Biclique Partition (cpp), {t} threads, results: vars={cpp_bp_vars}, clauses={cpp_bp_cls}, time={cpp_bp_time:.2f} ms")

    print('*'*50)


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
