#change the following according to your needs
import numpy as np

from ThredgeCorr.degree_dist import pk, pk_asymptotic
from ThredgeCorr.basic_patterns import solve_t

from scipy.stats import binom


def simulation_code(kw):

    n = int(kw['n'])
    rho = float(kw['rho'])
    k = float(kw['k'])

    kmax = n-1

    if k < n-1:
        t = solve_t(k, n)
        if rho > 0.0:
            p_exact = pk(n, t, rho, kmax)
            p_asymptotic = pk_asymptotic(np.arange(1,kmax,dtype=float), n, t, rho)
        else:
            P = k / (n-1.0)
            rv = binom(n-1, P)
            p_exact = rv.pmf(np.arange(kmax+1))
            p_asymptotic = np.array([])
    else:
        p_exact = np.array([])
        p_asymptotic = np.array([])

    return [p_exact, p_asymptotic]
