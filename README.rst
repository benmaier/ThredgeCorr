ThredgeCorr
===========

Generate instances of the thresholded locally-correlated edge weight
network model.

Install
-------

Use pip:

.. code:: bash

   pip install ThredgeCorr

Examples
--------

First steps
~~~~~~~~~~~

Create an instance of the network model.

.. code:: python

   from ThredgeCorr import get_networkx_graph
   N_nodes = 100
   rho = 0.49
   mean_degree = 1.2

   G = get_networkx_graph(N_nodes, covariance, mean_degree=mean_degree)

   # Alternatively, set this up with a threshold
   G = get_networkx_graph(N_nodes, covariance, threshold=0.5)

Theory
~~~~~~

To calculate network properties as explained in the paper, check out the
modules

.. code:: python

   ThredgeCorr.basic_patterns
   ThredgeCorr.degree_dist

and the scripts in
`https://github.com/benmaier/ThredgeCorr/tree/master/paper_figures`_ .

Use the more advanced class
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

   N_nodes = 100
   covariance = 0.35
   mean_degree = 2.0

   from ThredgeCorr import ThredgeCorrConstructor
   from scipy.special import binom

   TC = ThredgeCorrConstructor(N_nodes, covariance, mean_degree=mean_degree)

   # Generate edgelist
   edgelist = TC.get_edgelist()

   # generate sparse csr adjacency matrix
   A_csr = TC.get_adjacency_matrix(sparse=True) # default = False

   # generate networkx-Graph
   G = TC.get_graph()


   # Estimate the transitivity by sampling a network of 3 nodes `N` times,
   # then calculating `T = 3 * n_triangles / (3 * n_triangles + n_chains)`.
   T = TC.estimate_transitivity(int(binom(N_nodes,3))*500)

   # Change the threshold directly, or indirectly using the desired mean degree

   TC.set_threshold(2.253)
   TC.set_mean_degree(1.2)

   TC.update_covariance(0.2)

   # Compute degree sequences of 500 instances of the current configuration.
   ks = [ TC.get_adjacency_matrix().sum(axis=1).flatten() for n in range(500) ]

   # Compute degree sequence comparable to 500 instances of `N = 100` nodes.
   ks = TC.estimate_degree_sequence(N_nodes*500)

   # Compute the mean degree
   TC.get_mean_degree() # mean degree as given by the threshold

.. _`https://github.com/benmaier/ThredgeCorr/tree/master/paper_figures`: https://github.com/benmaier/ThredgeCorr/tree/master/paper_figures

