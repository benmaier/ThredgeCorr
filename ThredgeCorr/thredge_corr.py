from __future__ import print_function
import numpy as np
import scipy as sp
import scipy.sparse as sprs
from scipy.special import binom
import networkx as nx

from scipy.optimize import newton
from scipy.optimize import root

from ThredgeCorr.tools import get_degrees_from_edge_list, get_adjacency_matrix_from_edge_list, ccdf

class ThredgeCorrConstructor:

    def __init__(self,N,covariance,mean_degree=None,threshold=None):

        if (mean_degree is None) and (threshold is not None):
            self.t = threshold
        elif (mean_degree is not None) and (threshold is None):
            p = mean_degree / (N-1.0)
            self.t = sp.optimize.newton(lambda x: ccdf(x) - p, 0)

        self.N = N
        self.m = int(N*(N-1)/2)

        self.update_covariance(covariance)

        self.X = None

    def get_edgelist(self):

        z = np.random.normal(size=self.N)
        a, b = np.sqrt(1-2*self.covariance), np.sqrt(self.covariance)
        edges = [ (i,i+1+j) \
                            for i in range(self.N-1) \
                            for j in np.argwhere((   a*np.random.normal(size=(self.N-1-i)) \
                                                   + b*(z[i]+z[i+1:])\
                                                  ) > self.t).T[0] \
                ]
        return edges

    def set_threshold(self, threshold):
        self.t = threshold

    def set_mean_degree(self, mean_degree):
        p = mean_degree / (self.N-1.0)
        self.t = sp.optimize.newton(lambda x: ccdf(x) - p, 0)

    def get_mean_degree(self):
        p = ccdf(self.t)
        return p * (self.N-1)

    def update_covariance(self,covariance):
        self.covariance = covariance

    def get_adjacency_matrix(self,sparse=False):
        edges = self.get_edgelist()

        if sparse:
            rowcol = np.array(edges,dtype=int)
            if len(rowcol) == 0:
                return sprs.csr_matrix(([],([],[])),shape=(self.N,self.N),dtype=int)
            A = sprs.csr_matrix((np.ones(len(edges)),(rowcol[:,0], rowcol[:,1])),shape=(self.N,self.N),dtype=int)
            A += A.T
            return A
        else:
            A = np.zeros((self.N, self.N))
            for i, j in edges:
                A[i, j] = 1
            A += A.T
            return A

    def get_graph(self):
        G = nx.Graph()
        G.add_nodes_from(range(self.N))
        G.add_edges_from(self.get_edgelist())
        return G

    def estimate_transitivity(self,N_measurements,chunksize=None):

        if chunksize is None:
            chunksize = N_measurements

        self.C_triangle = np.ones((3,3)) * self.covariance
        np.fill_diagonal(self.C_triangle, 1)
        self.L_triangle = np.linalg.cholesky(self.C_triangle)

        n_triangles = 0
        n_chains = 0

        n_chunks = int(N_measurements // chunksize)

        for chunk in range(n_chunks+1):

            if chunk == n_chunks:
                chunksize = int(N_measurements % chunksize)
                if chunksize == 0:
                    break

            X = self.L_triangle.dot(np.random.randn(3,chunksize))

            n_edges = np.array(X>=self.t,dtype=int).sum(axis=0)
            n_triangles += np.count_nonzero(n_edges==3)
            n_chains += np.count_nonzero(n_edges==2)

        return 3.0 * float(n_triangles) / (3 * float(n_triangles) + float(n_chains))

    def estimate_triangles(self,N_measurements,chunksize=None):

        if chunksize is None:
            chunksize = N_measurements

        self.C_triangle = np.ones((3,3)) * self.covariance
        np.fill_diagonal(self.C_triangle, 1)
        self.L_triangle = np.linalg.cholesky(self.C_triangle)

        n_triangles = 0

        n_chunks = int(N_measurements // chunksize)

        for chunk in range(n_chunks+1):

            if chunk == n_chunks:
                chunksize = int(N_measurements % chunksize)
                if chunksize == 0:
                    break

            X = self.L_triangle.dot(np.random.randn(3,chunksize))

            n_edges = np.array(X>=self.t,dtype=int).sum(axis=0)
            n_triangles += np.count_nonzero(n_edges==3)

        return float(n_triangles)

    def estimate_degree_sequence(self,N_measurements,chunksize=None):
        
        if chunksize is None:
            chunksize = N_measurements

        self.C_one_node = np.ones((self.N,self.N)) * self.covariance
        np.fill_diagonal(self.C_one_node, 1)
        self.L_one_node = np.linalg.cholesky(self.C_one_node)
        
        n_chunks = int(N_measurements // chunksize)

        k = []

        for chunk in range(n_chunks+1):

            if chunk == n_chunks:
                chunksize = int(N_measurements % chunksize)
                if chunksize == 0:
                    break

            X = self.L_one_node.dot(np.random.randn(self.N,chunksize))
            k.extend( np.array(X>=self.t,dtype=int).sum(axis=0).tolist() )

        return np.array(k)


if __name__ == "__main__":

    N = 2000
    k = 10
    TCG = ThredgeCorrConstructor(N,0.5,mean_degree=k)

    G = TCG.get_graph()
    print("desired k =", k)
    print("<k> =", np.mean([ d[1] for d in G.degree() ]))
