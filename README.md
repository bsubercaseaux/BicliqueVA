# BicliqueVA
Bounded Variable Addition based on Biclique decompositions

## Installation

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


## Running Experiments

For a single basic experiment, run
```
uv run python3 -m bicliqueVA.experiments.basic
```
For a suite of experiments on random graphs (Erd\H{os}-Renyi), run
```
uv run python3 -m bicliqueVA.experiments.scaling_random_graphs
```