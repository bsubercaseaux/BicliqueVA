# BicliqueVA
Bounded Variable Addition based on Biclique partitions

## Installation

Requirements for a full execution: `kissat` in its version 4.0.4, must be executable as `kissat` in your path. `bva`, from e.g., (https://fmv.jku.at/bva/) executable as `bva` in your path, and `sbva` (https://github.com/hgarrereyn/SBVA) executable as `sbva` in your path.

### Automatic Installation

To install BicliqueVA, simply run the installation script:

```bash
bash install.sh
```

This will automatically install dependencies, compile the C++ components, and install the `biva` command globally.

### Manual Installation

If you prefer to install manually:

1. Install [uv](https://github.com/astral-sh/uv).
2. Sync dependencies:
   ```bash
   uv sync
   ```
3. Compile the C++ components:
   ```bash
   cd src/bicliqueVA/partition_algs/cpp && make && cd -
   ```
4. Install the command-line tool:
   ```bash
   uv tool install . --force
   ```

### Development (Editable) Installation

If you are developing BicliqueVA and want the `biva` command to automatically reflect changes in the source code, run:

```bash
bash install.sh --editable
```

Or manually:
```bash
uv tool install --editable . --force
```

This installs the package in "editable" mode, linking the global `biva` command directly to your repository folder.

After this, you should be able to run e.g.,

```bash
biva -f myformula.cnf
```
and see results.


## Running Experiments


For a suite of experiments on random graphs (Erd\H{os}-Renyi), run
```
uv run python3 -m bicliqueVA.experiments.scaling_random_graphs
```
