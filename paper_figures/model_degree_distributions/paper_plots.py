from load_data import data
import qsuite_config as cf
import numpy as np

import bfmplot as bp
from bfmplot import pl

import matplotlib as mpl
from matplotlib.ticker import MaxNLocator

mpl.rcParams['font.size'] = 9
ks = cf.ks
ns = cf.ns
rhos = cf.rhos
fig, axes = pl.subplots(1,2,figsize=(8,2.5))

LEGN = 0
LEGK = 0
RHON = 0
RHOK = 2

labels = ['0', '0.05', '0.15', '0.5']
pos = [ 67,3.5,2,2]

ax = axes[0]

for iax, ax in enumerate(axes):

    for i_n, n in enumerate(ns):
        for i_mk, mk in enumerate(ks):
            for i_rho, rho in enumerate(rhos):


                if mk < n-1:

                    if rho == 0:
                        k_ex = np.arange(n)
                    else:
                        k_ex = np.arange(n-1)
                    
                    k_as = np.arange(1,n-1)

                    p_ex = data[i_mk][i_n][i_rho][0]
                    p_as = data[i_mk][i_n][i_rho][1]

                    if len(p_ex) > 0:
                        if rho == 0:
                            col = c=bp.brewer_qualitative[2]
                            ls = '-'
                            label = r'random graph'
                        else:
                            ls = '-'
                            col = bp.brewer_qualitative[0]
                            label = r'exact solution'
                        if iax != 1 or i_rho > 1:
                            label = None
                        ax.plot(k_ex, p_ex, c=col, ls=ls, label = label)
                    if len(p_as) > 0:
                        if iax != 1 or i_rho > 1:
                            label = None
                        else:
                            label = r'asymptotic solution'
                        ax.plot(k_as, p_as, ls='--', c=bp.brewer_qualitative[1], label = label)


            
            ax.set_xlabel('node degree $k$')
            ax.set_ylabel('probability $p_k$')

            if iax == 0:
                ax.set_xscale('log')
                ax.set_yscale('log')
                ax.set_ylim(1e-5,1)
                ax.set_xlim(1,1.8e3)
            else:
                ax.set_ylim(0,0.05)
                ax.set_xlim(0,140)


        if iax == 1:
            leg = ax.legend()
            bp.align_legend_right(leg)

        #axes[0,i_n].set_title(r'$n={:d}$'.format(int(n)))
    bp.strip_axis(ax)

fig.tight_layout()
#pl.subplots_adjust(wspace=0.1,hspace=0.1)

#axes[2,0].axis('off')
#axes[3,0].axis('off')

ax = axes[0]
for i_n, n in enumerate(ns):
    for i_mk, mk in enumerate(ks):
        #print(list(xticklabels))
        #xticklabels[-1] = ''
        #ax.set_xticklabels(xticklabels)
        """
        if i_mk == len(ks) -1:
            pl.sca(ax)
            _locs, _labels = pl.xticks()
            if i_n < len(ns) - 1:
                _locs = _locs[1:-2]
                _labels = _labels[1:-2]
            else:
                _locs = _locs[2:-1]
                _labels = _labels[2:-1]
            #print(_locs, _labels[-1])
            #_labels[-1] = ''
            pl.xticks(_locs, _labels)

        if i_n == 0 and i_mk > 0:
            pl.sca(ax)
            _locs, _labels = pl.yticks()
            print(_locs, _labels)
            _locs = _locs[1:-2]
            _labels = _labels[1:-2]
            print(_locs, _labels)
            pl.yticks(_locs, _labels)
        """

        for i_rho, rho in enumerate(rhos):
            if rho == 0:
                k_ex = np.arange(n)
            else:
                k_ex = np.arange(n-1)
            p_ex = np.array([ float(p) for p in data[i_mk][i_n][i_rho][1]])

ax.text(34, 1e-4, r'$\rho=0$')
ax.text(3.5, 1e-4, r'$\rho=0.05$')
ax.text(1.2, 6e-3 , r'$\rho=0.15$')
ax.text(1.8, 9e-2, r'$\rho=0.5$')
ax.text(0.8, 1.0 , '(a)', ha='right',va='top',transform = ax.transAxes)
ax = axes[1]
ax.text(107, 0.037, r'$\rho=0$')
ax.text(43, 0.011, r'$\rho=0.05$')
ax.text(14.2, 0.016 , r'$\rho=0.15$')
ax.text(5.6, 0.024, r'$\rho=0.5$')
ax.text(1.0, 1.0 , '(b)', ha='right',va='top',transform = ax.transAxes)

fig.savefig('example_degree_distributions.pdf')

pl.show()
                
                
