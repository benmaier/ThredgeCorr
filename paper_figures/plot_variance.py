import numpy as np
#import matplotlib.pyplot as pl
from bfmplot import pl
import bfmplot as bp
from ThredgeCorr.basic_patterns import *
import os

fn = 'data/variance_data.npz'
data_exists = os.path.exists(fn)

if data_exists:
    data = np.load(fn)
    variances = data['variances']


### variance plot
n = 10**5
NMAX = 20 
rho = np.linspace(0.0,0.5,200)
k = 2**np.arange(0,10,2)
k = np.array(k, dtype=float)

if not data_exists:
    variances = np.zeros((len(rho),len(k)))
    for i in range(len(rho)):
        for j in range(len(k)):
            t = solve_t(k[j],n)
            variances[i,j] = variance(t,rho[i],n,NMAX)

fig, ax = pl.subplots(1,1,figsize=(4,2.5))
pl.xlabel(r'local edge weight correlation $\rho$')
pl.ylabel(r'node degree variance $\mathrm{Var}(k)$')

for i in range(len(k)):
    pl.plot(rho, variances[:,i], c=bp.brewer_qualitative[1])

bp.strip_axis(ax)
ax.set_yscale('log')
#ax.set_xscale('log')
ax.set_xlim([0,0.5])
fig.tight_layout()

for i in range(len(k)):
    bp.add_curve_label(ax, rho, variances[:,i], r'$\left\langle k\right\rangle={:d}$'.format(int(k[i])), label_pos_rel=0.4)

if not data_exists:
    np.savez(fn,variances=variances)

#pl.legend(loc=4)
pl.savefig('figures/variance.pdf')


pl.show()
