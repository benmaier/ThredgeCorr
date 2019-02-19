import numpy as np
#import matplotlib.pyplot as plt
from bfmplot import pl as plt
from ThredgeCorr.basic_patterns import *

### mean degree plots
n = 10**5
t = np.linspace(-6,6,800)
density = np.array([ edges(tt) for tt in t ]).astype(float)
k = (n-1)*density

fig = plt.figure(1,figsize=(8,4))
ax = fig.add_subplot(1, 2, 1)
ax.plot(t,density)
plt.xlim(-3,3)
plt.xlabel(r'$t$')
plt.ylabel(r'$1 - \Phi\left( t \right)$')

ax = fig.add_subplot(1, 2, 2)
ax.plot(t,k)
plt.xlim(2.5,4.5)
plt.ylim(float((n-1)*edges(4.5)*0.9), float(1.1*(n-1)*edges(2.5)) )
ax.set_yscale('log')
plt.xlabel(r'$t$')
plt.ylabel(r'$E \left[ k \right]$')

plt.savefig('figures/edges.pdf')


### variance plot
n = 10**5
NMAX = 20 
rho = np.linspace(0,0.5,200)
k = 2**np.arange(0,10,2)
k = np.array(k, dtype=float)
variances = np.zeros((200,5))
for i in range(200):
    for j in range(5):
        t = solve_t(k[j],n)
        variances[i,j] = variance(t,rho[i],n,NMAX)

fig = plt.figure(2,figsize=(6,6))
ax = fig.add_subplot(1, 1, 1)
plt.xlabel(r'$\rho$',size=16)
plt.ylabel(r'Var[k]',size=14)

for i in range(5):
    plt.plot(rho, variances[:,i], label=r'$k=$' + str(k[i]) )

ax.set_yscale('log')
#ax.set_xscale('log')
plt.legend(loc=4)
plt.savefig('figures/variance.pdf')


### triangle plots

n = 10**5
NMAX = 20 # increase this if needed

k = 10**np.linspace(0,4.99,20)
rho = np.linspace(0,0.5,6)
triangle_density = np.zeros((20,6))

for i in range(20):
    for j in range(6):
        t = solve_t(k[i],n)
        triangle_density[i,j] =  triangles(t,rho[j],NMAX)

fig = plt.figure(3,figsize=(12,6))
ax = fig.add_subplot(1, 2, 1)
plt.xlabel(r'$k$',size=16)
plt.ylabel(r'density of triangles',size=14)

for i in range(6):
    plt.plot(k, triangle_density[:,i], label=r'$\rho=$' + str(rho[i]) )

ax.set_xscale('log')
ax.set_yscale('log')
plt.legend(loc=4)


NMAX = 20
k = 2 
rho = np.linspace(0,0.5,6)
n = 10**np.linspace(2,5,20)
triangle_per_node = np.zeros((20,6))

print(rho)

for i in range(20):
    for j in range(6):
        t = solve_t(k,n[i])
        triangle_per_node[i,j] = (n[i]-1)*(n[i]-2)*triangles(t,rho[j],NMAX)

ax = fig.add_subplot(1, 2, 2)
plt.xlabel(r'$n$',size=16)
plt.ylabel(r'triangles per node',size=14)

for i in range(6):
    plt.plot(n, triangle_per_node[:,i], label=r'$\rho=$' + str(rho[i]) )

ax.set_xscale('log')
ax.set_yscale('log')
#plt.legend(loc=3)

plt.savefig('figures/triangles.pdf')

#plt.show()
plt.show()
