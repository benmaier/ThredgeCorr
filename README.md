# ThredgeCorr

Generate instances of George's thresholded edge-correlated network model.
This currently works only with the Cholesky decomposition and hence is only possible
to use up to network sizes of ≈ 150 Nodes.

## Install

Clone this repository

    git clone git@github.com:benmaier/ThredgeCorr.git

Install as development version (such that you don't have to reinstall after updating the repository)

    pip install -e ./ThredgeCorr --no-binary :all:

Alternatively, install as normal

    pip install ./ThredgeCorr

## Examples

First, import `numpy` and the `ThredgeCorrGraph` class.

```python
import numpy as np
from ThredgeCorr import ThredgeCorrGraph
```

Create an instance of the generating class.

```python
N_nodes = 100
beta = 0.49
mean_degree = 1.2

G = ThredgeCorrGraph(N_nodes, beta, mean_degree)

# Alternatively, set this up with a threshold
G = ThredgeCorrGraph(N_nodes, beta, threshold = 2.5)

# get an edge list from this generating class
edges = G.get_new_edge_list()
```

Change the threshold directly, or indirectly using the desired mean degree

```python
G.set_threshold(2.3)

# or

G.set_mean_degree(5.0)
```

Update the covariance (takes a while since the Cholesky decomposition has to be run again).

```python
G.update_covariance(0.2)
```

Compute degree sequences of 500 instances of the current configuration.

```python
ks = [ G.get_new_adjacency_matrix().sum(axis=1).flatten() for n in range(500) ]
```

You can also check out the values of the current matrices and parameters

```python
G.t # threshold
G.b # edge covariance
G.C # Covariance matrix
G.L # Cholesky decomposition
G.X # vector of multivariate Gaussian random variables
```

If you want to compute the edge index `e` of a node pair `(i, j)` (with condition `i < j`),
or the node index pair `(i, j)` from an edge index i, do

```python
e = G.edge_index(i, j)
i, j = G.node_indices(e)
```


