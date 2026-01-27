from utils.formula_gen import random_2cnf
from partition_algs.externals import run_factor_on_formula, run_bva_on_formula, run_cpp_bp_on_formula
from partition_algs.biclique_partition import run_biclique_partition_on_formula

if __name__ == "__main__":

    DATA = {}
    TIMEOUT = 30
    N_THREADS = 8
    
    n_vals = list(range(1000, 5001, 500))
    p_vals = [0.1, 0.3, 0.5, 0.7, 0.9]
    JSON_DATA = { 'experimental_results': []}
    algs = {
        'bva': run_bva_on_formula,
        'factor': run_factor_on_formula,
        'cpp_bp': run_cpp_bp_on_formula
    }

    additional_params = {
        'cpp_bp': {'n_threads': N_THREADS}
    }

    for n in n_vals:
        DATA[n] = {}
        for p in p_vals:
            data_n_p = {}
            formula = random_2cnf(n, p)
            print("-"*30)
            print(f"Random CNF generated. vars={formula.n_vars()}, clauses={formula.n_clauses()}")
            formula_filename = f"tmp/random_2cnf_{n}_{p}.cnf"
            formula.serialize(formula_filename)
            data_n_p["original"] = {"vrs": formula.n_vars(), "cls": formula.n_clauses()}

            for alg_name in algs.keys():
                additional_params_alg = additional_params.get(alg_name, dict())
                try:
                    alg_vars, alg_cls, alg_time = algs[alg_name](formula_filename, timeout=TIMEOUT, **additional_params_alg)
                    print(f"{alg_name} results: vars={alg_vars}, clauses={alg_cls}, time={alg_time}")
                    data_n_p[alg_name] = {"vrs": alg_vars, "cls": alg_cls, "time": alg_time}
                except TimeoutError:
                    print(f"{alg_name} timed out! (timeout set to {TIMEOUT})")
                    data_n_p[alg_name] = f'TIMEOUT ({TIMEOUT})'


            json_result = {'n_vertices': n, 'p_edge': p, 'outcomes': data_n_p}
            DATA[n][p] = data_n_p
            JSON_DATA['experimental_results'].append(json_result)
    with open("experimental_data.json", "w") as results_file:
        results_file.write(str(JSON_DATA))
