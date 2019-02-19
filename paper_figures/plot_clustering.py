import numpy as np
#import matplotlib.pyplot as pl
from bfmplot import pl as pl
from ThredgeCorr.basic_patterns import *
import bfmplot as bp
import matplotlib
import os

data_exists = os.path.exists('data/clustering_data.npz')

if data_exists:
    data = np.load('data/clustering_data.npz')
    clustering_k = data['clustering_k']
    clustering_n = data['clustering_n']
    clustering_r = data['clustering_r']
    clustering_T = data['clustering_T']

    data_arrays = [ clustering_k, clustering_n, clustering_r, clustering_T ]

matplotlib.rcParams['font.size']= 12

### triangle plots
fig, axes = pl.subplots(1,4,figsize=(12,3),
        #sharey=True
        )

##
axes[0].set_xlabel(r'mean degree $\left\langle k\right\rangle$')
axes[0].set_ylabel(r'clustering coefficient $C$')


for ax in axes:
    ax.set_yscale('log')

for ax in axes[:2]:
    ax.set_xscale('log')

axes[3].set_xscale('log')

axes[1].set_xlabel(r'number of nodes $n$')
#axes[1].set_ylabel(r'mean global clustering $\left\langle C\right\rangle$')

axes[2].set_xlabel(r'local edge weight correlation $\rho$')
#axes[2].set_ylabel(r'mean global clustering $\left\langle C\right\rangle$')

axes[3].set_xlabel(r'number of nodes $n$')
axes[3].set_ylabel(r'triangles per node $T$')
##

ax = axes[0]
col = bp.brewer_qualitative[1]

n = 10**5
NMAX = 20 # increase this if needed

k_0 = 10**np.linspace(0,4.99,50)
rho_0 = np.linspace(0,0.5,6)

ax.set_xlim([k_0.min(), k_0.max()] )

if not data_exists:
    clustering_k = np.zeros((len(k_0),len(rho_0)))

    for i in range(len(k_0)):
        for j in range(len(rho_0)):
            t = solve_t(k_0[i],n)
            clustering_k[i,j] =  triangles(t,rho_0[j],NMAX) / two_stars(t,rho_0[j],NMAX)

for i in range(len(rho_0)):
    ax.plot(k_0, clustering_k[:,i], label=r'$\rho={0:3.1f}$'.format(rho_0[i]),c=col,lw=1.5)

ax.set_ylim([1e-5,1])

#################################3

ax = axes[1]

NMAX = 20
k_1 = 4 
rho_1 = np.linspace(0,0.5,6)
n_1 = 10**np.linspace(2,5,50)
ax.set_xlim([n_1.min(), n_1.max()] )

if not data_exists:
    clustering_n = np.zeros((len(n_1),len(rho_1)))

    for i in range(len(n_1)):
        for j in range(len(rho_1)):
            t = solve_t(k_1,n_1[i])
            clustering_n[i,j] = triangles(t,rho_1[j],NMAX) / two_stars(t,rho_1[j],NMAX)


for i in range(6):
    ax.plot(n_1, clustering_n[:,i], label=r'$\rho={0:3.1f}$'.format(rho_1[i]),c=col,lw=1.5)

ax.set_ylim([1e-5,1])


################################
ax = axes[2]

NMAX = 20
n_2 = 10**5
rho_2 = np.linspace(0,0.5,50)
k_2 = np.logspace(0,8,5,base=4)
ax.set_xlim([rho_2.min(), rho_2.max()] )
if not data_exists:
    clustering_r = np.zeros((len(rho_2),len(k_2)))

    for i in range(len(rho_2)):
        for j in range(len(k_2)):
            t = solve_t(k_2[j],n_2)
            clustering_r[i,j] = triangles(t,rho_2[i],NMAX) / two_stars(t,rho_2[i],NMAX)

for i in range(len(k_2)):
    ax.plot(rho_2, clustering_r[:,i], label=r'$\left\langle k=\right\rangle$' + str(k_2[i]),c=col,lw=1.5)
    
ax.set_ylim([1e-5,1])
##################################
ax = axes[3]

ax.set_xlim([n_1.min(), n_1.max()] )

if not data_exists:
    clustering_T = np.zeros((len(n_1),len(rho_1)))

    for i in range(len(n_1)):
        for j in range(len(rho_1)):
            t = solve_t(k_1,n_1[i])
            clustering_T[i,j] = (n_1[i]-1)*(n_1[i]-2)*triangles(t,rho_1[j],NMAX) * 0.5


col = bp.brewer_qualitative[2]
for i in range(6):
    ax.plot(n_1, clustering_T[:,i], label=r'$\rho={0:3.1f}$'.format(rho_1[i]),c=col)


for ax in axes:
    bp.strip_axis(ax)


fig.tight_layout()
#pl.subplots_adjust(wspace=0.15)

for i in range(len(rho_0)):
    label_pos_rel = 0.2 - i/60.0
    va = 'center'
    if i == len(rho_1)-1:
        #label_pos_rel = 0.8
        va = 'baseline'
    bp.add_curve_label(axes[0], 
                       k_0, clustering_k[:,i], 
                       r'$\rho={0:3.1f}$'.format(rho_0[i]),
                       label_pos_rel = label_pos_rel,
                       bbox_pad = 0.1,
                       va = va,
                       fontsize=9,
                    )

for i in range(len(rho_1)):
    label_pos_rel = 0.8 + i/60.0
    va = 'center'
    if i == len(rho_1)-1:
        #label_pos_rel = 0.8
        va = 'baseline'
    bp.add_curve_label(axes[1], 
                       n_1, clustering_n[:,i],                        
                       r'$\rho={0:3.1f}$'.format(rho_1[i]),
                       label_pos_rel = label_pos_rel,
                       bbox_pad=0.5,
                       va = va,
                       fontsize=9,
                    )

for i in range(len(k_2)):
    if i == 0:
        s = r'$\left\langle k \right\rangle=1$'
    else:
        s = r'$\left\langle k \right\rangle=4^{:d}$'.format(2*i)
    bp.add_curve_label(axes[2], 
                       rho_2, clustering_r[:,i],                     
                       s,
                       label_pos_rel = 0.2 + i /50.0,
                       fontsize=9,
                    )

for i in range(len(rho_1)):
    label_pos_rel = 0.7 + i/60.0
    va = 'center'
    #if i == len(rho_1)-1:
    #    #label_pos_rel = 0.8
    #va = 'baseline'
    bp.add_curve_label(axes[3], 
                       n_1, clustering_T[:,i],                        
                       r'$\rho={0:3.1f}$'.format(rho_1[i]),
                       label_pos_rel = label_pos_rel,
                       bbox_pad=1,
                       va = va,
                       fontsize=9,
                    )


if not data_exists:
    np.savez('data/clustering_data.npz',
            clustering_k=clustering_k, 
            clustering_n=clustering_n, 
            clustering_r=clustering_r,
            clustering_T=clustering_T
            )

axes[0].text(0.95, 0.05, '(a)', va='bottom', ha='right', transform=axes[0].transAxes)
axes[1].text(0.05, 0.05, '(b)', va='bottom', ha='left', transform=axes[1].transAxes)
axes[2].text(0.95, 0.05, '(c)', va='bottom', ha='right', transform=axes[2].transAxes)
axes[3].text(0.05, 0.05, '(d)', va='bottom', ha='left', transform=axes[3].transAxes)
pl.savefig('figures/clustering.pdf')
pl.show()
