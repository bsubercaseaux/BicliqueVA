import time
import os
import tempfile
import subprocess
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
    # kissat uses SAT competition return codes:
    # 0 = UNKNOWN, 10 = SATISFIABLE, 20 = UNSATISFIABLE.
    # For factor preprocessing we can still use the emitted output CNF for all of them.
    valid_codes = {0, 10, 20}
    assert return_code in valid_codes, (
        f"factor (kissat) failed with return code {return_code}. Output:\n{output}"
    )
    if not os.path.exists(output_file_path):
        # Be robust if kissat solved early and skipped writing. Emit an equivalent trivial CNF.
        with open(output_file_path, "w") as f:
            if return_code == 20:
                f.write("p cnf 0 1\n0\n")
            else:
                f.write("p cnf 0 0\n")
    with open(output_file_path, 'r') as f:
        output = f.readline()
        tokens = output.split(' ')
        n_vars = int(tokens[2])
        n_cls = int(tokens[3])
  
    return n_vars, n_cls, elapsed_time

def run_bva_on_formula(formula_file_path, output_file_path="bva_out.cnf", timeout=None):
    command = ["bva", "-method=BVA2", f"-file={formula_file_path}", "-print=1", "-limit=-1"]
    
    start_time = time.perf_counter_ns()
    try:
        # Run bva directly, capturing stdout (CNF) and stderr (stats) separately
        result = subprocess.run(command, capture_output=True, text=True, timeout=timeout)
        elapsed_time = (time.perf_counter_ns() - start_time) / 1e6 # convert to milliseconds
        
        if result.returncode != 0:
             # If return code is not 0, dump all output for debugging
             assert result.returncode == 0, f"BVA failed with return code {result.returncode}. Output:\n{result.stderr}\n{result.stdout}"

    except subprocess.TimeoutExpired:
        raise TimeoutError

    # Process Stdout (CNF): mimic grep -v -e UNK -e "r "
    filtered_lines = []
    for line in result.stdout.splitlines():
        if "UNK" in line: continue
        if "r " in line: continue 
        filtered_lines.append(line)
    
    with open(output_file_path, 'w') as f:
        f.write("\n".join(filtered_lines) + "\n")

    # Process Stderr (Stats)
    # Expecting line: vars: 13 csl: 24
    stderr_output = result.stderr
    last_line = stderr_output.strip().split("\n")[-1]
    tokens = last_line.split(' ')
    n_vars = int(tokens[1])
    n_cls = int(tokens[3])
    
    return n_vars, n_cls, elapsed_time

def run_sbva_on_formula(formula_file_path, output_file_path="sbva_out.cnf", timeout=None):
    command = ["sbva", "-i", formula_file_path, "-o", output_file_path]
    if timeout is not None:
        command.extend(["-t", str(timeout)])
    
    start_time = time.perf_counter_ns()
    try:
        # Run bva directly, capturing stdout (CNF) and stderr (stats) separately
        sub_timeout = timeout + 1 if timeout is not None else None # extra second to avoid conflict with sbva's own timeout
        result = subprocess.run(command, capture_output=True, text=True, timeout=sub_timeout)
        elapsed_time = (time.perf_counter_ns() - start_time) / 1e6 # convert to milliseconds
        
        if result.returncode != 0:
             # If return code is not 0, dump all output for debugging
             assert result.returncode == 0, f"BVA failed with return code {result.returncode}. Output:\n{result.stderr}\n{result.stdout}"

    except subprocess.TimeoutExpired:
        raise TimeoutError


    with open(output_file_path, 'r') as f:
        output = f.readline()
        tokens = output.split(' ')
        n_vars = int(tokens[2])
        n_cls = int(tokens[3])
    
    return n_vars, n_cls, elapsed_time


def run_cpp_bp_on_formula(formula_file_path, output_file_path="bva_out.cnf", n_threads=1, timeout=None):
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

def _compute_remaining_timeout(timeout, elapsed_time_ms):
    if timeout is None:
        return None
    remaining = timeout - (elapsed_time_ms / 1000.0)
    if remaining <= 0:
        raise TimeoutError
    return remaining

def run_bva_over_bp(formula_file_path, output_file_path="bva_out.cnf", n_threads=1, timeout=None):
    output_dir = os.path.dirname(os.path.abspath(output_file_path))
    fd, bp_output_file_path = tempfile.mkstemp(prefix="bp_intermediate_", suffix=".cnf", dir=output_dir)
    os.close(fd)
    try:
        _, _, bp_time = run_cpp_bp_on_formula(
            formula_file_path,
            output_file_path=bp_output_file_path,
            n_threads=n_threads,
            timeout=timeout,
        )
        bva_timeout = _compute_remaining_timeout(timeout, bp_time)
        bva_vars, bva_cls, bva_time = run_bva_on_formula(
            bp_output_file_path,
            output_file_path=output_file_path,
            timeout=bva_timeout,
        )
    finally:
        if os.path.exists(bp_output_file_path):
            os.remove(bp_output_file_path)

    return bva_vars, bva_cls, bp_time + bva_time

def run_factor_over_bp(formula_file_path, output_file_path="bva_out.cnf", n_threads=1, timeout=None):
    output_dir = os.path.dirname(os.path.abspath(output_file_path))
    fd, bp_output_file_path = tempfile.mkstemp(prefix="bp_intermediate_", suffix=".cnf", dir=output_dir)
    os.close(fd)
    try:
        _, _, bp_time = run_cpp_bp_on_formula(
            formula_file_path,
            output_file_path=bp_output_file_path,
            n_threads=n_threads,
            timeout=timeout,
        )
        factor_timeout = _compute_remaining_timeout(timeout, bp_time)
        factor_vars, factor_cls, factor_time = run_factor_on_formula(
            bp_output_file_path,
            output_file_path=output_file_path,
            timeout=factor_timeout,
        )
    finally:
        if os.path.exists(bp_output_file_path):
            os.remove(bp_output_file_path)

    return factor_vars, factor_cls, bp_time + factor_time

if __name__ == "__main__":
    formula_file = "tmp/random_2cnf_1000_0.6.cnf"
    output_file = "bva_example.cnf"

    vars_bva, cls_bva, time_bva = run_bva_on_formula(formula_file, output_file, timeout=300)
    print(f"BVA Result: vars={vars_bva}, clauses={cls_bva}, time={time_bva:.2f} ms")

    n_vars, n_cls, elapsed_time = run_cpp_bp_on_formula(formula_file, output_file, timeout=300)
    n_vars2, n_cls2, elapsed_time2 = run_bva_on_formula(output_file, timeout=300)
    print(f"Bp Result: vars={n_vars}, clauses={n_cls}, time={elapsed_time:.2f} ms")
    print(f"BVA on Bp Result: vars={n_vars2}, clauses={n_cls2}, time={elapsed_time2:.2f} ms")
