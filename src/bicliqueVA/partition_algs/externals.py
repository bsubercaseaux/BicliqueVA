import time
from subprocess import TimeoutExpired, check_output, CalledProcessError, STDOUT

def system_call(command, timeout=None):
    """
    params:
        command: list of strings, ex. ["ls", "-l"]
        timeout: number of seconds to wait for the command to complete.
    returns: output, return_code
    """
    try:
        cmd_output = check_output(command, stderr=STDOUT, timeout=timeout).decode()
        return_code = 0
    except CalledProcessError as e:
        cmd_output = e.output.decode()
        return_code = e.returncode
    except TimeoutExpired:
        cmd_output = f"Command timed out after {timeout} seconds"
        return_code = (
            -4
        )  # You can use any number that is not a valid return code for a success or normal failure
    return cmd_output, return_code


def timed_run_shell(commands, timeout=None):
    """
    params:
        command: list of strings, ex. ["ls", "-l"]
        timeout: number of seconds to wait for the command to complete.
    returns: output, return_code, elapsed_time in mseconds
    """
    start_time = time.perf_counter_ns()
    output, return_code = system_call(commands, timeout=timeout)
    elapsed_time = time.perf_counter_ns() - start_time
    return output, return_code, elapsed_time / 1e6  # convert to milliseconds

def run_factor_on_formula(formula_file_path, output_file_path="bva_out.cnf", timeout=None):
    command = ["kissat", "--no-lucky", "--no-congruence", "--no-sweep", "--no-backbone", formula_file_path,
               "--no-conflicts", "--no-fastel", "-o", output_file_path, "-q"]
    output, return_code, elapsed_time = timed_run_shell(command, timeout=timeout)
    if return_code == -4:
        raise TimeoutError
    assert return_code == 0, f"factor (kissat) failed with return code {return_code}. Output:\n{output}"
    with open(output_file_path, 'r') as f:
        output = f.readline()
        tokens = output.split(' ')
        n_vars = int(tokens[2])
        n_cls = int(tokens[3])
  
    return n_vars, n_cls, elapsed_time

def run_bva_on_formula(formula_file_path, output_file_path="bva_out.cnf", timeout=None):
    command = ["sh", "bvar.sh", formula_file_path, output_file_path]
    output, return_code, elapsed_time = timed_run_shell(command, timeout=timeout)
    if return_code == -4:
        raise TimeoutError
    assert return_code == 0, f"BVA failed with return code {return_code}. Output:\n{output}"
    last_line = output.strip().split("\n")[-1]
    tokens = last_line.split(' ')
    n_vars = int(tokens[1])
    n_cls = int(tokens[3])
    return n_vars, n_cls, elapsed_time

def run_cpp_bp_on_formula(formula_file_path, output_file_path="bva_out.cnf", n_threads=1, timeout=None):
    import os
    # Get the directory of this file
    current_dir = os.path.dirname(os.path.abspath(__file__))
    # Construct path to bp binary
    bp_path = os.path.join(current_dir, "cpp", "bp")
    command = [bp_path, formula_file_path, output_file_path, str(n_threads)]
    output, return_code, elapsed_time = timed_run_shell(command, timeout=timeout)
    if return_code == -4:
        raise TimeoutError
    assert return_code == 0, f"cpp code failed with return code {return_code}. Output:\n{output}"
    print(output)
    with open(output_file_path, 'r') as f:
        output = f.readline()
        tokens = output.split(' ')
        n_vars = int(tokens[2])
        n_cls = int(tokens[3])
    return n_vars, n_cls, elapsed_time

if __name__ == "__main__":
    formula_file = "tmp/random_2cnf_1000_0.6.cnf"
    output_file = "bva_example.cnf"

    vars_bva, cls_bva, time_bva = run_bva_on_formula(formula_file, output_file, timeout=300)
    print(f"BVA Result: vars={vars_bva}, clauses={cls_bva}, time={time_bva:.2f} ms")

    n_vars, n_cls, elapsed_time = run_cpp_bp_on_formula(formula_file, output_file, timeout=300)
    n_vars2, n_cls2, elapsed_time2 = run_bva_on_formula(output_file, timeout=300)
    print(f"Bp Result: vars={n_vars}, clauses={n_cls}, time={elapsed_time:.2f} ms")
    print(f"BVA on Bp Result: vars={n_vars2}, clauses={n_cls2}, time={elapsed_time2:.2f} ms")