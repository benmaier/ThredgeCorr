    
import numpy as np
from scipy.special import binom
from ThredgeCorr import FastThredgeCorrGraph
from ThredgeCorr import ThredgeCorrGraph
from ThredgeCorr import get_degrees_from_edge_list
import networkx as nx
import matplotlib.pyplot as pl



from time import time
N = 150 

N_meas = 10


start = time()
A = ThredgeCorrGraph(N,0.49,.5)
[ A.get_new_edge_list() for n in range(N_meas) ]
end = time()

print("std method; N =", N, '; generating ', N_meas, 'networks took', end-start, 'seconds')

start = time()
numpymethod = FastThredgeCorrGraph(N,0.49,.5)
[ numpymethod.get_new_edge_list() for n in range(N_meas) ]
end = time()

print("numpy method; N =", N, '; generating ', N_meas, 'networks took', end-start, 'seconds')


k1 = []
C1 = []
k2 = []
C2 = []

#np_edges = T.get_n_edge_lists(500)

for meas in range(1000):
    edges = A.get_new_edge_list()
    ks = get_degrees_from_edge_list(N,edges).tolist()
    #k1.extend( ks.tolist())
    G = nx.Graph()
    G.add_nodes_from(range(N))
    G.add_edges_from(edges)

    ks = [ d[1] for d in list(G.degree())]
    k1.extend(ks)

    edges = numpymethod.get_new_edge_list()
    ks = get_degrees_from_edge_list(N,edges).tolist()
    #k1.extend( ks.tolist())
    G = nx.Graph()
    G.add_nodes_from(range(N))
    G.add_edges_from(edges)

    ks = [ d[1] for d in list(G.degree())]
    k2.extend(ks)

    #edges = np_edges[meas]
    #k2.extend(ks)

pl.hist(k1,bins=np.arange(max(k1)+1),histtype='step')
pl.hist(k2,bins=np.arange(max(k1)+1),histtype='step')
pl.xscale('log')
pl.yscale('log')

pl.show()

#print(np.array(k).mean())
#print("Transitivity from networks =", np.array(C1).mean())
#print("Transitivity from 3-cliques =", A.estimate_transitivity(int(binom(N,3))*50))
