from ThredgeCorr.basic_patterns import *
from ThredgeCorr.degree_dist import *
import numpy as np
import networkx as nx
import bfmplot as bp
from bfmplot import pl

import sys

fig = pl.figure(1,figsize=(8,2.5))

### high school network
k = np.loadtxt('degree_sequences/comm50.degree_sequence')
k = np.array(k,dtype=int)
kmax = np.max(k)+1
n = len(k)
t = solve_t( np.mean(k) ,n )
rho = solve_rho( np.mean(k), np.mean(k**2), n )

ax = fig.add_subplot(1, 3, 1)
pl.hist(k,bins=np.arange(0,kmax,2),density=True,alpha=1,align='left',color=bp.brewer_qualitative[1],label='data',histtype='step')
p = pk(n,t,rho,kmax,X,W)
pl.plot(np.arange(kmax),p,lw=1,c=bp.brewer_qualitative[2],label='model fit')
pl.xlim(0,30)
ax.text(0.04, 0.95, '(a)', va='top', ha='left', transform=ax.transAxes)
bp.strip_axis(ax)
leg = pl.legend()
bp.align_legend_right(leg)

ax.set_ylabel('probability density $p_k$')
ax.set_xlabel('node degree $k$')

#pl.show()
#sys.exit(0)
### collab network
k = np.loadtxt('degree_sequences/cond-mat.degree_sequence')
k = np.array(k,dtype=int)
kmax = np.max(k)+1
n = len(k)
t = solve_t( np.mean(k) ,n )
rho = solve_rho( np.mean(k), np.mean(k**2), n )

ax = fig.add_subplot(1, 3, 2)
b = np.logspace(np.log10(1),np.log10(kmax+10), 12)
pl.hist(k,bins=np.arange(0,kmax,5),density=True,alpha=1,align='mid',color=bp.brewer_qualitative[1],histtype='step')
#pl.hist(k,bins=b,normed=True,alpha=0.5)
p = pk(n,t,rho,kmax,X,W)
pl.plot(np.arange(kmax),p,lw=1,c=bp.brewer_qualitative[2])

ax.set_yscale('log')
#ax.set_xscale('log')
#pl.xlim(-2,kmax+10)
pl.ylim(10**-6,0.5)

ax.text(0.8, 0.95, '(b)', va='top', ha='left', transform=ax.transAxes)
ax.set_xlabel('node degree $k$')
ax.set_ylabel('probability density $p_k$')
bp.strip_axis(ax)




### proteins
k = np.loadtxt('degree_sequences/out.reactome.degree_sequence')
k = np.array(k,dtype=int)

kmax = np.max(k)+1
n = len(k)
t = solve_t( np.mean(k) ,n )
rho = solve_rho( np.mean(k), np.mean(k**2), n )

p = pk(n,t,rho,kmax,X,W)


tail = 1
ax = fig.add_subplot(1,3,3)
b = np.logspace(np.log10(tail),np.log10(kmax+10), 12)
pl.hist(k[k>=tail],bins=b,density=True,alpha=1,color=bp.brewer_qualitative[1],histtype='step')
pl.plot(np.arange(tail,kmax),p[tail:]/(1-np.sum(p[:tail])),lw=1,c=bp.brewer_qualitative[2])
ax.set_yscale('log')
ax.set_xscale('log')
pl.xlim(tail,10**3)
ax.text(0.8, 0.95, '(c)', va='top', ha='left', transform=ax.transAxes)
bp.strip_axis(ax)

fig.tight_layout()
ax.set_xlabel('node degree $k$')
ax.set_ylabel('probability density $p_k$')
pl.savefig('figures/degree_dist_fits.pdf')

pl.show()

