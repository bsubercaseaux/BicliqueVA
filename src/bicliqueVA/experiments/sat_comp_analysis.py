import requests
import json
import zipfile
import lzma
import gzip
import bz2
import io
import os
import shutil
import statistics
import subprocess
import tempfile
from urllib.parse import urlparse
import matplotlib.pyplot as plt

try:
    from ..partition_algs.externals import run_factor_on_formula
except ImportError:
    import sys
    src_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    if src_root not in sys.path:
        sys.path.insert(0, src_root)
    from bicliqueVA.partition_algs.externals import run_factor_on_formula
# from process_2cnf import process_2cnf_part


plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Palatino Linotype', 'Palatino', 'New Century Schoolbook', 'Times New Roman', 'serif']

plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300

# plt.style.use('seaborn-talk') # Gives a nice grid background
plt.style.use("tableau-colorblind10")
plt.rcParams['axes.prop_cycle'] = plt.cycler(color=plt.cm.tab10.colors)

# MAX_VARS_CONSIDERED = 3000

# # 4. Optional: Tweak sizes for readability
# plt.rcParams['font.size'] = 14
# plt.rcParams['axes.labelsize'] = 18
# plt.rcParams['axes.titlesize'] = 20
# plt.rcParams['xtick.labelsize'] = 16
# plt.rcParams['ytick.labelsize'] = 16
# plt.rcParams['legend.fontsize'] = 18


# --- CONFIGURATION ---
FILE_WITH_LINKS = "track_main_2025.uri"
FORMULAS_DIR = "sat_comp_formulas"

# Set to an integer (e.g., 5) to only process the first N files for testing.
# Set to None to process ALL files (will take a long time).
LIMIT_FILES = None

CLAUSE_LIMIT = 1_000_000_000

CHECK_CACHE = True

KNOWN_FORMULA_EXTENSIONS = {".cnf", ".gz", ".xz", ".zip", ".bz2", ".lzma"}
ZIP_CNF_SUFFIXES = (".cnf", ".cnf.gz", ".cnf.xz", ".cnf.bz2", ".cnf.lzma")
RUN_FACTOR = True
FACTOR_TIMEOUT_SECONDS = None
FACTOR_OUTPUT_DIR = os.path.join(FORMULAS_DIR, "factorized")
STREAMING_CLAUSE_THRESHOLD = 100_000_000


def strip_formula_extensions(name):
    """Strip common SAT formula/compression extensions to get a stable formula id."""
    stripped = os.path.basename(name)
    while True:
        base, ext = os.path.splitext(stripped)
        if ext.lower() in KNOWN_FORMULA_EXTENSIONS and base:
            stripped = base
            continue
        break
    return stripped


def resolve_links_file_path():
    """Resolve the URL list from cwd first, then from this script's directory."""
    cwd_candidate = FILE_WITH_LINKS
    script_candidate = os.path.join(os.path.dirname(__file__), FILE_WITH_LINKS)
    if os.path.exists(cwd_candidate):
        return cwd_candidate
    if os.path.exists(script_candidate):
        return script_candidate
    raise FileNotFoundError(f"Could not find links file '{FILE_WITH_LINKS}'.")


def load_formula_urls(links_file_path):
    with open(links_file_path, "r") as f:
        return [line.strip() for line in f if line.strip() and not line.lstrip().startswith("#")]


def to_int_width_dict(width_dict):
    if not width_dict:
        return {}
    return {int(k): int(v) for k, v in width_dict.items()}


def compute_factor_clause_savings(original_width_counts, factor_width_counts):
    """
    Compare clause counts per width before/after factor.
    Returns raw per-width delta and positive-only "saved" counts.
    """
    original = to_int_width_dict(original_width_counts)
    factored = to_int_width_dict(factor_width_counts)
    all_widths = sorted(set(original.keys()) | set(factored.keys()))

    delta_per_width = {}
    saved_per_width = {}
    for width in all_widths:
        delta = original.get(width, 0) - factored.get(width, 0)
        delta_per_width[width] = delta
        if delta > 0:
            saved_per_width[width] = delta

    total_saved_net = sum(delta_per_width.values())
    total_saved_positive = sum(saved_per_width.values())
    saved_pct_by_width = {}
    for width, saved in saved_per_width.items():
        if total_saved_positive > 0:
            saved_pct_by_width[width] = 100.0 * saved / total_saved_positive
        else:
            saved_pct_by_width[width] = 0.0

    return (
        delta_per_width,
        saved_per_width,
        total_saved_net,
        total_saved_positive,
        saved_pct_by_width,
    )


def run_factor_and_collect_stats(file_path, file_name, original_width_counts, verbosity=0):
    os.makedirs(FACTOR_OUTPUT_DIR, exist_ok=True)
    formula_id = strip_formula_extensions(file_name)
    factor_output_path = os.path.join(FACTOR_OUTPUT_DIR, f"{formula_id}_factor.cnf")
    factor_vars, factor_total_clauses, factor_time_ms = run_factor_on_formula(
        file_path,
        output_file_path=factor_output_path,
        timeout=FACTOR_TIMEOUT_SECONDS,
    )
    with open(factor_output_path, "rb") as ff:
        factor_stat = parse_cnf_content(ff, os.path.basename(factor_output_path))
    if factor_stat is None:
        raise RuntimeError(f"Failed to parse factor output for {file_name}")

    (
        factor_clause_delta_per_width,
        factor_saved_clauses_per_width,
        factor_total_saved_clauses,
        factor_saved_by_width_total,
        factor_saved_pct_by_width,
    ) = compute_factor_clause_savings(original_width_counts, factor_stat["n_clauses_per_width"])

    if verbosity >= 1:
        print(
            f"  Factor done for {file_name}: "
            f"clauses {factor_total_clauses}, net saved {factor_total_saved_clauses}"
        )

    return {
        "factor_n_vars": factor_vars,
        "factor_total_clauses": factor_total_clauses,
        "factor_avg_width": factor_stat["avg_width"],
        "factor_n_clauses_per_width": factor_stat["n_clauses_per_width"],
        "factor_time_ms": factor_time_ms,
        "factor_clause_delta_per_width": factor_clause_delta_per_width,
        "factor_saved_clauses_per_width": factor_saved_clauses_per_width,
        "factor_total_saved_clauses": factor_total_saved_clauses,
        "factor_saved_by_width_total": factor_saved_by_width_total,
        "factor_saved_pct_by_width": factor_saved_pct_by_width,
    }


def stat_from_cached_data(file_name, cached_data):
    stat = {
        "name": file_name,
        "n_vars": int(cached_data["n_vars"]),
        "total_clauses": int(cached_data["total_clauses"]),
        "binary_clauses": int(cached_data["binary_clauses"]),
        "avg_width": float(cached_data["avg_width"]),
        "n_clauses_per_width": cached_data["n_clauses_per_width"],
    }

    optional_factor_fields = [
        "factor_n_vars",
        "factor_total_clauses",
        "factor_avg_width",
        "factor_n_clauses_per_width",
        "factor_time_ms",
        "factor_clause_delta_per_width",
        "factor_saved_clauses_per_width",
        "factor_total_saved_clauses",
        "factor_saved_by_width_total",
        "factor_saved_pct_by_width",
    ]
    for key in optional_factor_fields:
        if key in cached_data:
            stat[key] = cached_data[key]
    return stat


def cache_payload_from_stat(stat):
    payload = {
        "n_vars": stat["n_vars"],
        "total_clauses": stat["total_clauses"],
        "binary_clauses": stat["binary_clauses"],
        "avg_width": stat["avg_width"],
        "n_clauses_per_width": stat["n_clauses_per_width"],
    }
    optional_factor_fields = [
        "factor_n_vars",
        "factor_total_clauses",
        "factor_avg_width",
        "factor_n_clauses_per_width",
        "factor_time_ms",
        "factor_clause_delta_per_width",
        "factor_saved_clauses_per_width",
        "factor_total_saved_clauses",
        "factor_saved_by_width_total",
        "factor_saved_pct_by_width",
    ]
    for key in optional_factor_fields:
        if key in stat:
            payload[key] = stat[key]
    return payload


def detect_compression(file_path):
    """Detect common compression types by magic bytes."""
    with open(file_path, "rb") as f:
        header = f.read(8)

    if header.startswith(b"\x1f\x8b"):
        return "gz"
    if header.startswith(b"\xfd7zXZ\x00"):
        return "xz"
    if header.startswith(b"BZh"):
        return "bz2"
    if header.startswith(b"PK\x03\x04"):
        return "zip"
    return None


def decompress_one_layer(source_path, target_path, compression):
    if compression == "gz":
        with gzip.open(source_path, "rb") as src, open(target_path, "wb") as dst:
            shutil.copyfileobj(src, dst)
        return
    if compression == "xz":
        with lzma.open(source_path, "rb") as src, open(target_path, "wb") as dst:
            shutil.copyfileobj(src, dst)
        return
    if compression == "bz2":
        with bz2.open(source_path, "rb") as src, open(target_path, "wb") as dst:
            shutil.copyfileobj(src, dst)
        return
    if compression == "zip":
        with zipfile.ZipFile(source_path) as zf:
            file_infos = [info for info in zf.infolist() if not info.is_dir()]
            if not file_infos:
                raise ValueError(f"Zip archive has no files: {source_path}")
            preferred = next(
                (
                    info for info in file_infos
                    if info.filename.lower().endswith(ZIP_CNF_SUFFIXES)
                ),
                file_infos[0],
            )
            with zf.open(preferred, "r") as src, open(target_path, "wb") as dst:
                shutil.copyfileobj(src, dst)
        return

    raise ValueError(f"Unsupported compression type: {compression}")


def ensure_plain_cnf(file_path, verbosity=0):
    """
    Decompress in place until file is plain text CNF (handles nested containers).
    Returns the number of decompression steps applied.
    """
    steps = 0
    while True:
        compression = detect_compression(file_path)
        if compression is None:
            return steps

        tmp_fd, tmp_path = tempfile.mkstemp(
            prefix=f"{os.path.basename(file_path)}.",
            suffix=".tmp_decompressed",
            dir=os.path.dirname(file_path) or ".",
        )
        os.close(tmp_fd)
        try:
            decompress_one_layer(file_path, tmp_path, compression)
            os.replace(tmp_path, file_path)
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
        steps += 1
        if verbosity >= 1:
            print(f"Decompressed {os.path.basename(file_path)} ({compression}).")

        if steps > 8:
            raise RuntimeError(f"Too many decompression layers for {file_path}")


def normalize_formula_file(formulas_dir, file_name, verbosity=0):
    """Normalize a formula file to '<formula_id>.cnf' and ensure it is uncompressed."""
    source_path = os.path.join(formulas_dir, file_name)
    if not os.path.isfile(source_path):
        return None
    if file_name.startswith("."):
        return None
    if ".tmp_decompressed" in file_name:
        return None

    formula_id = strip_formula_extensions(file_name)
    if not formula_id:
        return None

    target_path = os.path.join(formulas_dir, f"{formula_id}.cnf")
    normalized_path = source_path

    if source_path != target_path and not os.path.exists(target_path):
        os.rename(source_path, target_path)
        normalized_path = target_path
        if verbosity >= 1:
            print(f"Renamed {file_name} -> {os.path.basename(target_path)}")
    elif source_path != target_path:
        normalized_path = target_path
        if verbosity >= 1:
            print(f"Skipping rename for {file_name}, target exists: {os.path.basename(target_path)}")

    ensure_plain_cnf(normalized_path, verbosity=verbosity)
    return os.path.basename(normalized_path)


def normalize_formula_directory(formulas_dir, verbosity=0):
    os.makedirs(formulas_dir, exist_ok=True)
    normalized = set()
    for file_name in sorted(os.listdir(formulas_dir)):
        normalized_name = normalize_formula_file(formulas_dir, file_name, verbosity=verbosity)
        if normalized_name and normalized_name.endswith(".cnf"):
            normalized.add(normalized_name)
    return sorted(normalized)


def download_missing_formulas(urls, formulas_dir, verbosity=0):
    os.makedirs(formulas_dir, exist_ok=True)

    existing_formula_ids = set()
    for file_name in os.listdir(formulas_dir):
        full_path = os.path.join(formulas_dir, file_name)
        if not os.path.isfile(full_path):
            continue
        existing_formula_ids.add(file_name)
        existing_formula_ids.add(os.path.splitext(file_name)[0])
        existing_formula_ids.add(strip_formula_extensions(file_name))

    downloaded = 0
    skipped = 0
    for url in urls:
        formula_id = strip_formula_extensions(os.path.basename(urlparse(url).path))
        if not formula_id:
            if verbosity >= 1:
                print(f"Skipping invalid formula URL: {url}")
            continue

        if formula_id in existing_formula_ids:
            skipped += 1
            continue

        target_path = os.path.join(formulas_dir, f"{formula_id}.cnf")
        was_downloaded = download_file(url, target_path)
        if not was_downloaded:
            continue
        existing_formula_ids.add(formula_id)
        existing_formula_ids.add(f"{formula_id}.cnf")
        downloaded += 1

    print(f"Formula download pass complete. Downloaded: {downloaded}, already present: {skipped}.")

def download_file(url, filename):
    """Downloads a file with a simple progress indicator."""
    if os.path.exists(filename):
        print(f"File '{filename}' already exists. Skipping download.")
        return True

    print(f"Downloading from {url}...")
    response = requests.get(url, stream=True)
    total_size = int(response.headers.get('content-length', 0))
    block_size = 1024 * 1024  # 1MB
    downloaded = 0

    if total_size > 20* 1024 * 1024 * 1024:
        print(f"Total size: {total_size // (1024*1024)} MB > 20 GB; skipping for now")
        return False

    with open(filename, 'wb') as f:
        for data in response.iter_content(block_size):
            downloaded += len(data)
            f.write(data)
            if total_size > 0:
                percent = 100 * downloaded / total_size
                print(f"\rProgress: {percent:.1f}% ({downloaded // (1024*1024)} MB)", end='')
            else:
                print(f"\rDownloaded: {downloaded // (1024*1024)} MB", end='')
    print("\nDownload complete!")
    return True

def parse_cnf_content(file_stream, file_name):
    """
    Parses a DIMACS CNF file stream.
    Returns a dict with statistics or None if invalid.
    """
    clause_widths = []
    binary_count = 0
    n_clauses_per_width = {}
    header_n_vars = None
    header_n_clauses = None
    use_streaming_width_count = False
    running_width_sum = 0
    total_clause_count = 0
    
    # Determine how to open based on extension inside the zip
    try:
        if file_name.endswith('.xz'):
            f = lzma.open(file_stream, mode='rt', encoding='utf-8')
        elif file_name.endswith('.gz'):
            f = gzip.open(file_stream, mode='rt', encoding='utf-8')
        else:
            # Assume plain text
            f = io.TextIOWrapper(file_stream, encoding='utf-8')
    except Exception as e:
        print(f"Error opening {file_name}: {e}")
        return None

    try:
        two_cnf_part = []
        n_vars = 0
        for line in f:
            line = line.strip()
            # Skip comments and header
            if line.startswith('p cnf'):
                tokens = line.split()
                header_n_vars = int(tokens[2])
                clauses = int(tokens[3])
                header_n_clauses = clauses
                if clauses > CLAUSE_LIMIT:
                    print(f"  Skipping {file_name} with {clauses} clauses (>{CLAUSE_LIMIT}).")
                    return None
                if clauses >= STREAMING_CLAUSE_THRESHOLD:
                    use_streaming_width_count = True
                    clause_widths = None
                    print(
                        f"  Large formula {file_name} with {clauses} clauses: "
                        "using streaming width count mode."
                    )
            if not line or line.startswith('c') or line.startswith('p') or line.startswith('%'):
                continue
            
            # Parse literals
            # DIMACS lines end with 0, but sometimes 0 is on a new line or missing
            # We split by whitespace
            parts = line.split()
            
            # Filter out the trailing '0' which marks end of clause
            literals = [x for x in parts if x != '0']
            
            # Handle empty clause lines ("0"), which encode UNSAT.
            if not literals:
                if '0' in parts:
                    width = 0
                    if use_streaming_width_count:
                        running_width_sum += width
                    else:
                        clause_widths.append(width)
                    total_clause_count += 1
                    if width not in n_clauses_per_width:
                        n_clauses_per_width[width] = 0
                    n_clauses_per_width[width] += 1
                continue

            n_vars = max(n_vars, max([abs(int(l)) for l in literals]))
                
            width = len(literals)
            if use_streaming_width_count:
                running_width_sum += width
            else:
                clause_widths.append(width)
            total_clause_count += 1
            
            if width == 2:
                binary_count += 1
                two_cnf_part.append(list(map(int, literals)))

            if width not in n_clauses_per_width:
                n_clauses_per_width[width] = 0
            n_clauses_per_width[width] += 1
                
    except Exception as e:
        print(f"Error parsing {file_name}: {e}")
        return None
    finally:
        f.close() # Important to close the wrapper, not the underlying zip stream directly

    if total_clause_count == 0:
        # Valid edge-case: formula with zero clauses after preprocessing.
        if header_n_clauses == 0:
            return {
                "name": file_name,
                "n_vars": header_n_vars if header_n_vars is not None else n_vars,
                "total_clauses": 0,
                "binary_clauses": 0,
                "avg_width": 0.0,
                "n_clauses_per_width": {},
            }
        return None

    # component_stats = process_2cnf_part(two_cnf_part, file_name)

    return {
        "name": file_name,
        "n_vars": n_vars,
        "total_clauses": total_clause_count,
        "binary_clauses": binary_count,
        "avg_width": (
            running_width_sum / total_clause_count
            if use_streaming_width_count
            else statistics.mean(clause_widths)
        ),
        "n_clauses_per_width": n_clauses_per_width,
        # "component_stats": component_stats
    }


def main(verbosity=0):
    links_file_path = resolve_links_file_path()
    urls = load_formula_urls(links_file_path)
    if LIMIT_FILES is not None:
        urls = urls[:LIMIT_FILES]
    download_missing_formulas(urls, FORMULAS_DIR, verbosity=verbosity)
    formula_files = normalize_formula_directory(FORMULAS_DIR, verbosity=verbosity)
    os.makedirs("results", exist_ok=True)
    os.makedirs("results/cached", exist_ok=True)
    if RUN_FACTOR:
        os.makedirs(FACTOR_OUTPUT_DIR, exist_ok=True)

    # iterate over files in sat_comp_formulas directory
    stats_list = []
    processed_files = 0
    total_files = len(formula_files)
    for file in formula_files:
        if LIMIT_FILES is not None and processed_files >= LIMIT_FILES:
            break
        file_path = os.path.join(FORMULAS_DIR, file)
        stat = None

        # check if cached_results_{file} already exists
        cache_file = f"results/cached/cached_results_{file}.txt"
        if CHECK_CACHE and os.path.exists(cache_file):
            if verbosity >= 1:
                print(f"Loading cached processing for {file_path}.")
            with open(cache_file, 'r') as cf:
                cached_data = json.load(cf)
                required_cache_fields = {
                    "n_vars",
                    "total_clauses",
                    "binary_clauses",
                    "avg_width",
                    "n_clauses_per_width",
                }
                if cached_data and required_cache_fields.issubset(cached_data.keys()):
                    stat = stat_from_cached_data(file, cached_data)
                elif verbosity >= 1:
                    print(f"  Cache file {cache_file} is missing fields, recomputing.")

        print(f"Processing (#{processed_files + 1}/{total_files}) {file_path}...")
        try:
            if stat is None:
                with open(file_path, 'rb') as f:
                    stat = parse_cnf_content(f, file)
                if stat is None:
                    print(f"Failed to parse {file_path}.")
                    continue
            if verbosity >= 1:
                print(
                    f"  Total Clauses: {stat['total_clauses']}, "
                    f"Binary Clauses: {stat['binary_clauses']}, "
                    f"Avg Width: {stat['avg_width']:.2f}"
                )

            if RUN_FACTOR:
                has_factor_stats = (
                    "factor_n_clauses_per_width" in stat
                    and "factor_saved_clauses_per_width" in stat
                    and "factor_saved_pct_by_width" in stat
                )
                if not has_factor_stats:
                    if verbosity >= 1:
                        print(f"  Running factor on {file}...")
                    factor_stat = run_factor_and_collect_stats(
                        file_path,
                        file,
                        stat["n_clauses_per_width"],
                        verbosity=verbosity,
                    )
                    stat.update(factor_stat)

            stats_list.append(stat)
            processed_files += 1
            with open(cache_file, 'w') as cf:
                json.dump(cache_payload_from_stat(stat), cf)
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            raise e

        if (LIMIT_FILES is not None) and processed_files >= LIMIT_FILES:
            break
    if stats_list:
        # Compute global statistics
        
        global_n_clauses_per_width = {}
        for s in stats_list:
            for width, count in s['n_clauses_per_width'].items():
                w = int(width)
                if w not in global_n_clauses_per_width:
                    global_n_clauses_per_width[w] = 0
                global_n_clauses_per_width[w] += count

        # compute states based on global_n_clauses_per_width
        # we have a dict width -> count
        total_count = sum(global_n_clauses_per_width.values())
        sum_w = sum(w * c for w, c in global_n_clauses_per_width.items())
        sum_w2 = sum(w * w * c for w, c in global_n_clauses_per_width.items())

        global_avg_width = sum_w / total_count
        global_stdev_width = ((sum_w2 - (sum_w**2 / total_count)) / (total_count - 1))**0.5 if total_count > 1 else 0.0

        print(f"\nProcessed {len(stats_list)} files.")
        print(f"Global Clause Width Statistics:")
        print(f"  Average Width: {global_avg_width:.2f}")
        print(f"  Std Dev Width: {global_stdev_width:.2f}")
        
        width_histogram = {}
        for w, c in global_n_clauses_per_width.items():
            width_histogram[w] = c
        # for width in all_widths:
        #     width_histogram[width] = width_histogram.get(width, 0) + 1

        print("\nClause Width Histogram:")
        total_count = sum(width_histogram.values())

        more_than_8 = 0
        for width in sorted(width_histogram.keys()):
            if width <= 8:
                print(f"  Width {width}: {width_histogram[width]}")
            else:
                more_than_8 += width_histogram[width]
        plt.bar(list(range(1, 9)), [100 *(width_histogram.get(i, 0) / total_count) for i in range(1, 9)])
        # add additional bar corresponding to widths > 8
        if more_than_8 > 0:
            print(f"  Width >8: {more_than_8}")
            plt.bar(9,100 * (more_than_8 / total_count))
            plt.xticks(list(range(1, 10)), labels=[str(i) for i in range(1, 9)] + [">8"])
        else:
            plt.xticks(list(range(1, 9)), labels=[str(i) for i in range(1, 9)])

        factor_summary = None
        if RUN_FACTOR:
            global_factor_saved_clauses_per_width = {}
            global_factor_clause_delta_per_width = {}
            formulas_with_factor_stats = 0
            total_factor_saved_clauses_net = 0
            total_factor_saved_by_width_positive = 0
            for s in stats_list:
                if "factor_saved_clauses_per_width" not in s:
                    continue
                formulas_with_factor_stats += 1
                for width, count in s["factor_saved_clauses_per_width"].items():
                    w = int(width)
                    c = int(count)
                    global_factor_saved_clauses_per_width[w] = (
                        global_factor_saved_clauses_per_width.get(w, 0) + c
                    )
                for width, count in s["factor_clause_delta_per_width"].items():
                    w = int(width)
                    c = int(count)
                    global_factor_clause_delta_per_width[w] = (
                        global_factor_clause_delta_per_width.get(w, 0) + c
                    )
                total_factor_saved_clauses_net += int(s.get("factor_total_saved_clauses", 0))
                total_factor_saved_by_width_positive += int(s.get("factor_saved_by_width_total", 0))

            global_factor_saved_pct_by_width = {}
            for width, saved in global_factor_saved_clauses_per_width.items():
                if total_factor_saved_by_width_positive > 0:
                    global_factor_saved_pct_by_width[width] = (
                        100.0 * saved / total_factor_saved_by_width_positive
                    )
                else:
                    global_factor_saved_pct_by_width[width] = 0.0

            print("\nFactor Clause Savings by Width (positive savings only):")
            for width in sorted(global_factor_saved_clauses_per_width.keys()):
                saved = global_factor_saved_clauses_per_width[width]
                pct = global_factor_saved_pct_by_width[width]
                print(f"  Width {width}: saved {saved} clauses ({pct:.2f}%)")
            print(f"  Net clauses saved (orig - factor): {total_factor_saved_clauses_net}")
            print(
                f"  Positive per-width saved clauses total: "
                f"{total_factor_saved_by_width_positive}"
            )

            factor_summary = {
                "formulas_with_factor_stats": formulas_with_factor_stats,
                "global_factor_clause_delta_per_width": global_factor_clause_delta_per_width,
                "global_factor_saved_clauses_per_width": global_factor_saved_clauses_per_width,
                "global_factor_saved_pct_by_width": global_factor_saved_pct_by_width,
                "total_factor_saved_clauses_net": total_factor_saved_clauses_net,
                "total_factor_saved_by_width_positive": total_factor_saved_by_width_positive,
            }
        
        plt.xlabel("Clause Width")
        plt.ylabel("Percentage")
        # plt.xticks(range(1, 9))
        plt.title("Clause Width Distribution (SAT'2025 Comp Formulas)")
        plt.savefig("results/clause_width_histogram.png")
        plt.close()
        print(f"Writing results to results/experiment_results.json")

        results_dict = {
            "number of formulas": len(stats_list),
            "global_avg_width": global_avg_width,
            "global_stdev_width": global_stdev_width,
            "global_n_clauses_per_width": global_n_clauses_per_width,
            "factor_summary": factor_summary,
            "stats_list": stats_list
        }
        
        with open("results/experiment_results.json", 'w') as f:
            json.dump(results_dict, f, indent=4)
    
    
    

if __name__ == "__main__":
    main()
