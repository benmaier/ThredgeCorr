from __future__ import print_function
import numpy as np
import scipy as sp
import scipy.sparse as sprs
from scipy.special import binom

from scipy.optimize import newton
from scipy.optimize import root

from numpy.random import multivariate_normal

def ccdf(x):
    return 1 - 0.5*(1 + sp.special.erf(x/np.sqrt(2)))

def _F(X,n,beta):
    a,b,c = X
    ans = np.zeros(3)
    ans[0] = -1 + a**2 + 2*b**2*(-2 + n) + c**2*(-1 - 2*(-2 + n) + ((-1 + n)*n)/2.)
    ans[1] = 2*a*b + 2*b*c*(-3 + n) + b**2*(-2 + n) + (c**2*(12 - 7*n + n**2))/2. - beta
    ans[2] = 4*b**2 + 4*b*c*(-4 + n) + (c*(4*a + c*(20 - 9*n + n**2)))/2.
    return ans

def _jacobian(X,n,beta):
    a,b,c = X
    J = np.zeros((3,3))
    J[0,0] = 2*a
    J[0,1] = 4*b*(-2 + n)
    J[0,2] = 2*c*(-1 - 2*(-2 + n) + ((-1 + n)*n)/2.)
    J[1,0] = 2*b
    J[1,1] = 2*a + 2*c*(-3 + n) + 2*b*(-2 + n)
    J[1,2] = 2*b*(-3 + n) + 2*c*(-2*(-3 + n) - n + ((-1 + n)*n)/2.)
    J[2,0] = 2*c
    J[2,1] = 8*b + 4*c*(-4 + n)
    J[2,2] = 2*a + 4*b*(-4 + n) + 2*c*(-6 - 4*(-4 + n) + ((-1 + n)*n)/2.)
    return J

class ThredgeCorrGraph:

    def __init__(self,N,covariance,mean_degree=None,threshold=None,build_cholesky_matrix=True):


        if (mean_degree is None) and (threshold is not None):
            self.t = threshold
        elif (mean_degree is not None) and (threshold is None):
            p = mean_degree / (N-1.0)
            self.t = sp.optimize.newton(lambda x: ccdf(x) - p, 0)

        self.N = N
        self.m = int(N*(N-1)/2)

        self.update_covariance(covariance,build_cholesky_matrix)

        self.X = None

    def set_threshold(self, threshold):
        self.t = threshold

    def set_mean_degree(self, mean_degree):
        p = mean_degree / (self.N-1.0)
        self.t = sp.optimize.newton(lambda x: ccdf(x) - p, 0)

    def mean_degree(self):
        p = ccdf(self.t)
        return p * (self.N-1)

    def edge_index(self,i,j):
        return int(i*self.N + j - (i+2)*(i+1) * 0.5)

    def node_indices_slow(self,e):
        marker = e
        for i in range(self.N):
            delta = self.N - 1 - i
            if marker - delta < 0:
                break
            else:
                marker -= delta
        j = int(e - i * self.N + 0.5*(i+1)*(i+2))
        return i, j

    def node_indices(self,e):
        N = self.N

        lower_bound = (N-1.5) - np.sqrt( (N-1.5)**2 -2*e-4+2*N )
        i = int(np.ceil(lower_bound))
        j = int(e - i * N + 0.5*(i+1)*(i+2))

        return i,j

    def update_covariance_slow(self,covariance):
        """Create the covariance matrix and the cholesky matrix subsequently"""

        self.b = b = covariance

        C = np.eye(self.m)

        for e1 in range(self.m):
            i1, j1 = self.node_indices(e1)
            for e2 in range(e1+1, self.m):
                i2, j2 = self.node_indices(e2)

                if (i1 == i2) or (i1 == j2) or (j1 == i2) or (j1 == j2):
                    C[e1, e2] = b
                    C[e2, e1] = b


        self.C = C
        self.L = np.linalg.cholesky(C)

    def update_covariance(self,covariance,build_cholesky=True):
        """Create the covariance matrix and the cholesky matrix subsequently"""

        self.b = b = covariance

        C = np.eye(self.m)

        for node in range(self.N):
            for neighbor1 in range(node+1,self.N-1):
                edge1 = self.edge_index(node,neighbor1)
                for neighbor2 in range(neighbor1+1,self.N):
                    edge2 = self.edge_index(node,neighbor2)
                    edge3 = self.edge_index(neighbor1,neighbor2)

                    C[edge1, edge2] = b
                    C[edge2, edge1] = b
                    C[edge1, edge3] = b
                    C[edge3, edge1] = b
                    C[edge2, edge3] = b
                    C[edge3, edge2] = b
                    
        self.C = C
        if build_cholesky:
            self.L = np.linalg.cholesky(C)


    def generate_weight_vector(self):
        #self.X = self.L.dot( np.random.normal(0,1,self.m) )
        self.X = self.L.dot( np.random.randn(self.m) )

    def threshold_and_get_edgelist(self,threshold):
        if self.X is None:
            self.generate_weight_vector()

        ndx = np.where(self.X>=threshold)[0]

        edges = []
        for e in ndx:
            edges.append( self.node_indices(e) )

        return edges

    def get_new_edge_list(self):
        self.X = None
        return self.threshold_and_get_edgelist(self.t)

    def get_new_adjacency_matrix(self,sparse=False):
        self.X = None
        edges = self.threshold_and_get_edgelist(self.t)

        if sparse:
            rowcol = np.array(edges,dtype=int)
            if len(rowcol) == 0:
                return np.zeros(N)
            A = sprs.csr_matrix((np.ones(len(edges)),(rowcol[:,0], rowcol[:,1])),shape=(N,N),dtype=int)
            A += A.T
            return A
        else:
            A = np.zeros((self.N, self.N))
            for i, j in edges:
                A[i, j] = 1
            A += A.T
            return A

    def estimate_transitivity(self,N_measurements,chunksize=None):

        if chunksize is None:
            chunksize = N_measurements

        self.C_triangle = np.ones((3,3)) * self.b
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

        self.C_triangle = np.ones((3,3)) * self.b
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

        self.C_one_node = np.ones((self.N,self.N)) * self.b
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

class FastThredgeCorrGraph(ThredgeCorrGraph):


    def __init__(self,N,covariance,mean_degree=None,threshold=None):

        ThredgeCorrGraph.__init__(self,N,covariance,
                                  mean_degree = mean_degree,
                                  threshold = threshold,
                                  build_cholesky_matrix = False)

    def update_covariance(self,covariance,build_cholesky_matrix=False):

        self.b = covariance
        self.parameters = root(_F,[0.7, 0.02, -4e-5],jac=_jacobian,args=(self.N,self.b),tol=1e-16).x

    def generate_weight_vector(self):
        a,b,c = self.parameters
        A = np.zeros((self.N,self.N))
        y = np.random.normal(size=self.m)
        s = np.sum(y)
        w = np.zeros(self.N)
        for i in range(self.N):
            for j in range(i):
                w[i] += y[self.edge_index(j,i)]
            for j in range(i+1,self.N):
                w[i] += y[self.edge_index(i,j)]
        self.X = (a-2*b+c)*y + c*s
        for i in range(self.N):
            for j in range(i+1,self.N):
                self.X[self.edge_index(i,j)] += (b-c)*(w[i]+w[j])
	

class NumpyThredgeCorrGraph(ThredgeCorrGraph):

    def __init__(self,N,covariance,mean_degree=None,threshold=None):

        ThredgeCorrGraph.__init__(self,N,covariance,
                                  mean_degree = mean_degree,
                                  threshold = threshold,
                                  build_cholesky_matrix = False)

    def get_n_edge_lists(self,n):

        X = multivariate_normal(np.zeros(self.C.shape[:1]),self.C,size=n)

        edges = []
        for meas in range(n):
            ndx = np.where(X[meas,:]>self.t)[0]
            these_edges = [ self.node_indices(e) for e in ndx ]
            edges.append(these_edges)

        return edges


def get_degrees_from_edge_list(N,edges):
    rowcol = np.array(edges,dtype=int)
    if len(rowcol) == 0:
        return np.zeros(N)
    A = sprs.csr_matrix((np.ones(len(edges)),(rowcol[:,0], rowcol[:,1])),shape=(N,N),dtype=int)
    A += A.T
    k = np.asarray(A.sum(axis=1)).flatten()
    return k
        

if __name__ == "__main__":
    
    N = 100 

    N_meas = 1

    from time import time

    start = time()
    A = ThredgeCorrGraph(N,0.49,.5)
    [ A.get_new_edge_list() for n in range(N_meas) ]
    end = time()

    print("std method; N =", N, '; generating ', N_meas, 'networks took', end-start, 'seconds')

    #start = time()
    #numpymethod = NumpyThredgeCorrGraph(N,0.49,.5)
    #numpymethod.get_n_edge_lists(N_meas)
    #end = time()

    #print("numpy method; N =", N, '; generating ', N_meas, 'networks took', end-start, 'seconds')


    import networkx as nx
    import matplotlib.pyplot as pl

    k1 = []
    C1 = []
    k2 = []
    C2 = []

    #np_edges = T.get_n_edge_lists(500)
    
    for meas in range(50):
        edges = A.get_new_edge_list()
        ks = get_degrees_from_edge_list(N,edges).tolist()
        #k1.extend( ks.tolist())
        G = nx.Graph()
        G.add_nodes_from(range(N))
        G.add_edges_from(edges)

        ks = [ d[1] for d in list(G.degree())]
        C1.append(nx.transitivity(G))
        k1.extend(ks)


        #edges = np_edges[meas]
        #k2.extend(ks)

    #print(np.array(k).mean())
    print("Transitivity from networks =", np.array(C1).mean())
    print("Transitivity from 3-cliques =", A.estimate_transitivity(int(binom(N,3))*50))

    from rocsNWL.drawing import draw

    pl.figure()
    pl.hist(k1,histtype='step',bins=max(k1))
    k2 = A.estimate_degree_sequence(N*50)
    pl.hist(k2,histtype='step',bins=max(k2))
    pl.xscale('log')
    pl.yscale('log')
    print(np.mean(k1))

    #pl.figure()
    #pl.hist(C1,histtype='step')
    #pl.hist(C2,histtype='step')

    #pl.figure()

    #G = nx.Graph()
    #G.add_nodes_from(range(N))
    #G.add_edges_from(edges)
    ##draw(G,labels=list(range(N)))
    #draw(G)

    fig,ax = pl.subplots(1,2,figsize=(10,5))
    ax[0].set_title("the covariance matrix $C$")
    ax[0].imshow(A.C)
    ax[1].set_title("the covariance matrix $C^2$")
    C_squared = A.C.dot(A.C)
    ax[1].imshow(C_squared)

    #ax[1].set_title("sophisticated")
    #ax[1].spy(A.C2)

    print(A.L[-3:,-3:])

    pl.figure()
    pl.spy(A.L)

    pl.show()

    """
    for node in range(N-1):
        for neigh in range(node+1,N):
            print("===========")
            print(node, neigh)
            print(A.edge_index(node,neigh))
            print(A.node_indices(A.edge_index(node,neigh)))
            """
