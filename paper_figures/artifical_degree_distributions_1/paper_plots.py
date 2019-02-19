from load_data import data
import qsuite_config as cf
import numpy as np

import bfmplot as bp
from bfmplot import pl

import matplotlib as mpl
from matplotlib.ticker import MaxNLocator

mpl.rcParams['font.size'] = 9
ks = cf.ks[:-1]
ns = cf.ns[1:]
rhos = cf.rhos
fig, axes = pl.subplots(3,3,figsize=(9,6))

LEGN = 2
LEGK = 0
RHON = 0
RHOK = 2

labels = ['0', '0.025', '0.05', '0.1', '0.25', '0.5']
pos = [ 202,100,42,17,4,3]
for i_n, n in enumerate(ns):
    for i_mk, mk in enumerate(ks):
        ax = axes[i_mk, i_n]
        for i_rho, rho in enumerate(rhos):


            if mk < n-1:

                if rho == 0:
                    k_ex = np.arange(n)
                else:
                    k_ex = np.arange(n-1)
                
                k_as = np.arange(1,n-1)

                p_ex = data[i_mk][i_n+1][i_rho][0]
                p_as = data[i_mk][i_n+1][i_rho][1]

                if len(p_ex) > 0:
                    if rho == 0:
                        col = c=bp.brewer_qualitative[2]
                        ls = '-'
                        label = r'ER ($\rho=0$)'
                    else:
                        ls = '-'
                        col = bp.brewer_qualitative[0]
                        label = r'exact solution'
                    if i_n != LEGN or i_mk != LEGK or i_rho > 1:
                        label = None
                    ax.plot(k_ex, p_ex, c=col, ls=ls, label = label)
                if len(p_as) > 0:
                    if i_n != LEGN or i_mk != LEGK or i_rho > 1:
                        label = None
                    else:
                        label = r'asymptotic solution'
                    ax.plot(k_as, p_as, ls='--', c=bp.brewer_qualitative[1], label = label)


        
        if i_mk == len(ks) -1:
            ax.set_xlabel('node degree $k$')
        if i_n == 0:
            ax.set_ylabel('probability $p_k$')

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim(1e-6,1)

        if i_n == LEGN and i_mk == LEGK:
            ax.text(0.99,0.14,r'$\left\langle k\right\rangle={:d}$'.format(int(mk)),
                            va = 'bottom',
                            ha = 'right',
                            transform = ax.transAxes,
                            )
            ax.text(0.99,0.03,r'$n={:d}$'.format(int(n)),
                            va = 'bottom',
                            ha = 'right',
                            transform = ax.transAxes,
                            )
        else:
            ax.text(0.99,0.99,r'$\left\langle k\right\rangle={:d}$'.format(int(mk)),
                            va = 'top',
                            ha = 'right',
                            transform = ax.transAxes,
                            )
            ax.text(0.99,0.88,r'$n={:d}$'.format(int(n)),
                            va = 'top',
                            ha = 'right',
                            transform = ax.transAxes,
                            )
        bp.strip_axis(ax)

        if i_mk < len(ks) -1:

            _labels = [item.get_text() for item in ax.get_xticklabels()]

            empty_string_labels = ['']*len(_labels)
            ax.set_xticklabels(empty_string_labels)

        if i_n > 0:
            pl.sca(ax)
            pl.yticks([],[])

        if i_mk == LEGK and i_n == LEGN:
            leg = ax.legend()
            bp.align_legend_right(leg)

        ax.set_xlim(1,n-1)

    #axes[0,i_n].set_title(r'$n={:d}$'.format(int(n)))

fig.tight_layout()
pl.subplots_adjust(wspace=0.1,hspace=0.1)

#axes[2,0].axis('off')
#axes[3,0].axis('off')

for i_n, n in enumerate(ns):
    for i_mk, mk in enumerate(ks):
        ax = axes[i_mk, i_n]
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
            p_ex = np.array([ float(p) for p in data[i_mk][i_n+1][i_rho][0]])
            if i_n == RHON and i_mk == RHOK:
                bp.add_curve_label(ax, k_ex, p_ex, r'$\rho={}$'.format(labels[i_rho]), label_pos_abs = pos[i_rho])

fig.savefig('example_degree_distributions.pdf'.format(int(n)))

pl.show()
                
                
