import numpy as np
import scipy as sp
import networkx as nx
from scipy.stats import expon

from ThredgeCorr.tools import ccdf

def get_networkx_graph(N, covariance, threshold=None,mean_degree=None):
    """Create an instance of a networkx-Graph from the locally correlated
    Threshold model.

    Parameters
    ----------
    N : int
        Number of nodes
    covariance : float
        The local covariance of the generated edge weights
    threshold : float, default : None
        Edges are constructed if edge weights are greater than
        this threshold. If `None`, then ``mean_degree`` is
        expected to be given.
    mean_degree : float, default : None
        The desired mean degree. If given, the threshold will
        be compute from this value.
        If `None`, then ``Threshold`` is
        expected to be given.

    Returns
    -------
    G : networkx.Graph
        The constructed network.
    """

    if (mean_degree is None) and (threshold is not None):
        t = threshold
    elif (mean_degree is not None) and (threshold is None):
        p = mean_degree / (N-1.0)
        t = sp.optimize.newton(lambda x: ccdf(x) - p, 0)

    rho = covariance

    assert(0<=rho)
    assert(rho<=0.5)

    G, z = nx.empty_graph(N), np.random.normal(size=N)
    a, b = np.sqrt(1-2*rho), np.sqrt(rho)

    for i in range(N-1):
        for j in np.argwhere((a*np.random.normal(size=(N-1-i)) + b*(z[i]+z[i+1:])) > t).T[0]:
            G.add_edge(i,i+1+j)

    return G

from scipy.stats import norm

#def get_cdf_serrano(N, covariance, significance_threshold):
#
#    rho = covariance
#    G, z, y = nx.empty_graph(N), norm(size=N), np.random.normal(size=(N*(N-1))//2)
#    a, b = np.sqrt(1-2*rho), np.sqrt(rho)
#
#    edge_weights = np.zeros_like(y)
#    node_strengths = np.zeros_like(z)
#
#    for i in range(N-1):
#        edge_weights += np.argwhere((a*np.random.normal(size=(N-1-i)) + b*(z[i]+z[i+1:])
#            neigh = 
#            G.add_edge(i,i+1+j)
#
#    assert(0<=rho)
#    assert(rho<=0.5)
#
#def get_serrano(N, covariance, significance_threshold):
#    pass
#
#
#def get_serrano_thresholded_networkx_graph(N, covariance, significance_threshold, scipy_stats_distribution='expon'):
#
#    if scipy_stats_distribution == 'expon':
#        dist = expon()
#    else:
#        dist = scipy_stats_distribution
#
#    assert(0<=rho)
#    assert(rho<=0.5)
#
#    alpha_tilde = -np.log(significance_threshold)
#
#    G, z = nx.empty_graph(N), dist.rvs(size=N)
#    a, b = np.sqrt(1-2*rho), np.sqrt(rho)
#
#    for i in range(N-1):
#        for j in np.argwhere((a*dist.rvs(size=(N-1-i)) + b*(z[i]+z[i+1:])) > t).T[0]:
    

if __name__ == "__main__":

    k = 10
    N = 200
    N_meas = 100

    print("N =", N)
    print("desired k =", k)
    print("doing",N_meas,"measurements... ")

    all_k = 0
    for meas in range(N_meas):

        G = get_networkx_graph(N,0.5,mean_degree=k)
        all_k += np.mean([ d[1] for d in G.degree() ])

    print("<k> =", all_k/N_meas)

