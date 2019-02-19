import numpy as np
import scipy as sp
import scipy.sparse as sprs

from scipy.optimize import newton

def edge_index(N,i,j):
    return int(i*N + j - (i+2)*(i+1) * 0.5)

def node_indices_slow(N,e):
    marker = e
    for i in range(N):
        delta = N - 1 - i
        if marker - delta < 0:
            break
        else:
            marker -= delta
    j = int(e - i * N + 0.5*(i+1)*(i+2))
    return (i,j)

def node_indices_bounds(N,e):

    upper_bound = (N-.5) - np.sqrt( (N-.5)**2 -2*e)
    lower_bound = (N-1.5) - np.sqrt( (N-1.5)**2 -2*e-4+2*N )

    print("i_lower =", lower_bound, np.ceil(lower_bound))
    print("i_upper =", upper_bound)


def node_indices_fast(N,e):

    lower_bound = (N-1.5) - np.sqrt( (N-1.5)**2 -2*e-4+2*N )
    i = int(np.ceil(lower_bound))
    j = int(e - i * N + 0.5*(i+1)*(i+2))

    return (i,j)


    



if __name__ == "__main__":
    N = 10

    all_edges_slow = []
    all_edges_fast = []

    for node in range(N-1):
        for neigh in range(node+1, N):
            e = edge_index(N,node,neigh)
            i,j = node_indices_slow(N,e)
            all_edges_slow.append( i == node and j == neigh )
            i,j = node_indices_fast(N,e)
            all_edges_fast.append( i == node and j == neigh )

    print("slow algorithm found the right edges:", all(all_edges_slow))
    print("fast algorithm found the right edges:", all(all_edges_fast))


    from time import time
    from progressbar import ProgressBar as PB

    N_meas = 100

    N = 200
    e = []
    print("======== delta method (slow) ========")
    print("N_meas =", N_meas)

    start = time()

    bar = PB()
    for meas in bar(range(N_meas)):
        for node in range(N-1):
            for neigh in range(node+1, N):
                e = edge_index(N,node,neigh)
                i,j = node_indices_slow(N,e)

    end = time()
    print("took",end-start,"seconds")


    print("======== bounds method (fast) ========")
    print("N_meas =", N_meas)

    start = time()

    bar = PB()

    for meas in bar(range(N_meas)):
        for node in range(N-1):
            for neigh in range(node+1, N):
                e = edge_index(N,node,neigh)
                i,j = node_indices_fast(N,e)

    end = time()
    print("took",end-start,"seconds")


