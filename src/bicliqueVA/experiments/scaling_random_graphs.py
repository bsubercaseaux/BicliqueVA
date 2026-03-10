from bicliqueVA.partition_algs.externals import run_sbva_on_formula
import json
import math
import time
import subprocess
import random
random.seed(43)
import argparse
from collections import Counter
from ..utils.formula_gen import random_2cnf
from pysat.formula import CNF
from pysat.card import CardEnc, EncType
from ..partition_algs.externals import (
    run_factor_on_formula,
    run_bva_on_formula,
    run_cpp_bp_on_formula,
    run_bva_over_bp,
    run_factor_over_bp,
    run_sbva_on_formula
)
from ..partition_algs.biclique_partition import run_biclique_partition_on_formula
from ..partition_algs.depth_reencoder import run_depth_reencode


def add_at_least_k_to_cnf(input_formula_filename, output_formula_filename, base_n_vars, k):
    cnf = CNF(from_file=input_formula_filename)
    top_id = max(cnf.nv, base_n_vars)
    base_vars = list(range(1, base_n_vars + 1))
    card = CardEnc.atleast(
        lits=base_vars,
        bound=k,
        top_id=top_id,
        encoding=EncType.mtotalizer,
    )
    cnf.extend(card.clauses)
    cnf.nv = max(top_id, card.nv)
    cnf.to_file(output_formula_filename)


def run_sat_with_is_constraint(formula_filename, base_n_vars, k, timeout=None):
    constrained_formula_filename = formula_filename.replace(".cnf", f"_is{k}.cnf")
    add_at_least_k_to_cnf(
        input_formula_filename=formula_filename,
        output_formula_filename=constrained_formula_filename,
        base_n_vars=base_n_vars,
        k=k,
    )

    start_time = time.perf_counter_ns()
    try:
        result = subprocess.run(
            ["kissat", constrained_formula_filename],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            timeout=timeout,
        )
        sat_time_ms = (time.perf_counter_ns() - start_time) / 1e6
    except subprocess.TimeoutExpired:
        raise TimeoutError

    if result.returncode == 10:
        sat_status = "SAT"
    elif result.returncode == 20:
        sat_status = "UNSAT"
    else:
        sat_status = f"UNKNOWN({result.returncode})"

    return sat_time_ms, sat_status


def aggregate_trial_entries(entries, timeout):
    successful = [entry for entry in entries if isinstance(entry, dict)]
    n_trials = len(entries)
    n_success = len(successful)
    n_timeout = n_trials - n_success

    if n_success == 0:
        return f"TIMEOUT ({timeout})"

    aggregated = {
        "trials": n_trials,
        "n_success": n_success,
        "n_timeout": n_timeout,
    }

    numeric_fields = ["vrs", "cls", "time", "is_time"]
    for field in numeric_fields:
        values = [entry[field] for entry in successful if field in entry and isinstance(entry[field], (int, float))]
        if len(values) > 0:
            aggregated[field] = sum(values) / len(values)

    statuses = [entry["is_status"] for entry in successful if "is_status" in entry]
    if len(statuses) > 0:
        aggregated["is_status"] = Counter(statuses).most_common(1)[0][0]

    return aggregated


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Scaling experiments on random independent-set formulas.")
    parser.add_argument("--trials", type=int, default=1, help="Number of random formulas to average per (n,p).")
    args = parser.parse_args()
    if args.trials < 1:
        raise ValueError("--trials must be at least 1.")
    TRIALS = args.trials

    DATA = {}
    TIMEOUT = 300
    N_THREADS = 1
    
    n_vals = list(range(600, 3001, 600))
    p_vals = [0.5]
    JSON_DATA = { 'experimental_results': []}
    algs = {
        'bva': run_bva_on_formula,
        'factor': run_factor_on_formula,
        'cpp_bp': run_cpp_bp_on_formula,
        'bva_over_bp': run_bva_over_bp,
        'factor_over_bp': run_factor_over_bp,
        # "sbva": run_sbva_on_formula
        # 'depth_reencode': run_depth_reencode
    }

    additional_params = {
        'cpp_bp': {'n_threads': N_THREADS},
        'bva_over_bp': {'n_threads': N_THREADS},
        'factor_over_bp': {'n_threads': N_THREADS},
    }

    for n in n_vals:
        DATA[n] = {}
        for p in p_vals:
            print("-" * 30)
            is_k = max(1, 1+int(1.2 * math.log2(n)))
            print(f"Running n={n}, p={p} with trials={TRIALS}, independent-set bound k={is_k}")

            trial_entries = {"original": []}
            for alg_name in algs.keys():
                trial_entries[alg_name] = []

            for trial in range(TRIALS):
                formula = random_2cnf(n, p)
                formula_filename = f"tmp/random_2cnf_{n}_{p}_trial{trial + 1}.cnf"
                formula.serialize(formula_filename)
                print(
                    f"[trial {trial + 1}/{TRIALS}] Random CNF generated. "
                    f"vars={formula.n_vars()}, clauses={formula.n_clauses()}"
                )

                original_data = {"vrs": formula.n_vars(), "cls": formula.n_clauses()}
                try:
                    sat_time, sat_status = run_sat_with_is_constraint(
                        formula_filename,
                        base_n_vars=n,
                        k=is_k,
                        timeout=TIMEOUT,
                    )
                    original_data["is_time"] = sat_time
                    original_data["is_status"] = sat_status
                    print(f"[trial {trial + 1}] original+at_least_k runtime: {sat_time:.2f} ms ({sat_status})")
                except TimeoutError:
                    print(f"[trial {trial + 1}] original+at_least_k timed out! (timeout set to {TIMEOUT})")
                    trial_entries["original"].append(f"TIMEOUT ({TIMEOUT})")
                else:
                    trial_entries["original"].append(original_data)

                for alg_name in algs.keys():
                    additional_params_alg = dict(additional_params.get(alg_name, dict()))
                    output_formula_filename = f"tmp/random_2cnf_{n}_{p}_{alg_name}_trial{trial + 1}.cnf"
                    output_param_name = "output_formula_filename" if alg_name == "depth_reencode" else "output_file_path"
                    additional_params_alg[output_param_name] = output_formula_filename
                    try:
                        alg_vars, alg_cls, alg_time = algs[alg_name](
                            formula_filename,
                            timeout=TIMEOUT,
                            **additional_params_alg,
                        )
                        print(
                            f"[trial {trial + 1}] {alg_name} results: "
                            f"vars={alg_vars}, clauses={alg_cls}, time={alg_time}"
                        )
                        alg_data = {"vrs": alg_vars, "cls": alg_cls, "time": alg_time}

                        try:
                            sat_time, sat_status = run_sat_with_is_constraint(
                                output_formula_filename,
                                base_n_vars=n,
                                k=is_k,
                                timeout=TIMEOUT,
                            )
                            alg_data["is_time"] = sat_time
                            alg_data["is_status"] = sat_status
                            print(f"[trial {trial + 1}] {alg_name}+at_least_k runtime: {sat_time:.2f} ms ({sat_status})")
                        except TimeoutError:
                            alg_data["is_time"] = f"TIMEOUT ({TIMEOUT})"
                            alg_data["is_status"] = "TIMEOUT"
                            print(f"[trial {trial + 1}] {alg_name}+at_least_k timed out! (timeout set to {TIMEOUT})")

                        trial_entries[alg_name].append(alg_data)
                    except TimeoutError:
                        print(f"[trial {trial + 1}] {alg_name} timed out! (timeout set to {TIMEOUT})")
                        trial_entries[alg_name].append(f"TIMEOUT ({TIMEOUT})")

            data_n_p = {"is_k": is_k, "trials": TRIALS}
            data_n_p["original"] = aggregate_trial_entries(trial_entries["original"], timeout=TIMEOUT)
            for alg_name in algs.keys():
                data_n_p[alg_name] = aggregate_trial_entries(trial_entries[alg_name], timeout=TIMEOUT)


            json_result = {'n_vertices': n, 'p_edge': p, 'outcomes': data_n_p}
            DATA[n][p] = data_n_p
            JSON_DATA['experimental_results'].append(json_result)
    with open("experimental_data.json", "w") as results_file:
        json.dump(JSON_DATA, results_file, indent=4)
