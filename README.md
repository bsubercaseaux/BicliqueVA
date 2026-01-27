# BicliqueVA
Bounded Variable Addition based on Biclique decompositions

## Installation

```
pip3 install eznf
```

Also:

```
cd partition_algs/cpp && make
```


## Running Experiments

For a single basic experiment, run
```
python3 -m experiments.basic
```

For a suite of experiments on random graphs (Erd\H{os}-Renyi), run
```
python3 -m experiments.scaling_random_graphs
```