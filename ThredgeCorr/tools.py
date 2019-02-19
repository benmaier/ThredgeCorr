from __future__ import print_function
import numpy as np
import scipy as sp
import scipy.sparse as sprs

def get_degrees_from_edge_list(N,edges):
    rowcol = np.array(edges,dtype=int)
    if len(rowcol) == 0:
        return np.zeros(N)
    A = sprs.csr_matrix((np.ones(len(edges)),(rowcol[:,0], rowcol[:,1])),shape=(N,N),dtype=int)
    A += A.T
    k = np.asarray(A.sum(axis=1)).flatten()
    return k

def get_adjacency_matrix_from_edge_list(N,edges,sparse=False):

    if sparse:
        rowcol = np.array(edges,dtype=int)
        if len(rowcol) == 0:
            return sprs.csr_matrix(([],([],[])),shape=(N,N),dtype=int)
        A = sprs.csr_matrix((np.ones(len(edges)),(rowcol[:,0], rowcol[:,1])),shape=(N,N),dtype=int)
        A += A.T
        return A
    else:
        A = np.zeros((N, N))
        for i, j in edges:
            A[i, j] = 1
        A += A.T
        return A

def ccdf(x):
    return 1 - 0.5*(1 + sp.special.erf(x/np.sqrt(2)))

