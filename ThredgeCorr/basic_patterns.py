from mpmath import mp
import numpy as np

### precision: increase this if higher precision necessary
mp.dps = 30

### standard functions, as defined in the text
def phi(x):
    return( mp.npdf(x) )

def Phi(x):
    return( mp.ncdf(x) )

def Phi_inv(x):
    return( -mp.sqrt(2)*mp.erfinv(mp.mpf(1)-mp.mpf(2)*x) )

def hermite(n,x):
    '''
    probabilists' Hermite polynomials, H_n(x)
    '''
    return( mp.power(2,-mp.mpf(0.5)*n)*mp.hermite(n,x/mp.sqrt(2)) )


### density fuctions
def edges(t):
    return( mp.mpf(1)-Phi(t) )

def two_stars(t,rho,NMAX):
    # compute Hermite polynomials 
    Hm1 = (mp.mpf(1)-Phi(t))/phi(t)
    H = [hermite(N,t) for N in range(NMAX+1)] + [Hm1]
    # compute partial sum (up to NMAX)
    ans = mp.mpf(0)
    for N in range(NMAX+1):
        ans = ans + mp.power(rho,N) * mp.power(H[N-1],2) / mp.factorial(N)
    ans = ans*mp.power(phi(t),2)
    return( ans )

def triangles(t,rho,NMAX):
    # compute Hermite polynomials and factorials
    Hm1 = (mp.mpf(1)-Phi(t))/phi(t)
    H = [hermite(N,t) for N in range(NMAX+1)] + [Hm1]
    factorial = [mp.factorial(i) for i in range(NMAX+1)]
    # compute partial sum (up to NMAX)
    ans = mp.mpf(0)
    for N in range(NMAX+1):
        rhoN = mp.power(rho,N)
        for i in range(N+1):
            for j in range(N+1-i):
                k = N-i-j
                ans = ans + rhoN * ( (H[N-i-1]/factorial[i]) *
                        (H[N-j-1]/factorial[j]) * (H[N-k-1]/factorial[k]))
    ans = ans*mp.power(phi(t),3)
    return(ans)


### variance of degree distribution plots
def variance(t,rho,n,NMAX):
    return( (n-1)*(n-2)*two_stars(t,rho,NMAX) + (n-1)*edges(t) - ((n-1)*edges(t))**2 )


### find parameter from desired mean degrees
def dts_drho(t,rho,NMAX):
    H = [mp.power(hermite(N,t),2) for N in range(NMAX)] 
    ans = mp.mpf(0)
    for N in range(NMAX):
        ans = ans +  mp.power(rho,N) * H[N] / mp.factorial(N)
    return( ans*phi(t)*phi(t) )

def solve_t(k,n):
    return( - Phi_inv( mp.mpf(k)/mp.mpf(n-1) ) )

def solve_rho(k,k2,n,NMAX=40):
    t = solve_t(k,n)
    mu = mp.mpf(k2-k)/mp.mpf((n-1)*(n-2))
    rho = 0.25
    drho = 1.
    for i in range(50):
        drho = -(two_stars(t,rho,NMAX) - mu)/dts_drho(t,rho,NMAX)
        rho = rho + drho
    return(rho)



### functions for fitting:
def p0(n,t,rho):
    return( pk(n,t,rho,1,X,W)[0] )

def mean_degree(n,t):
    return((n-1)*edges(t))

def mean_degree_sq(n,t,rho,NMAX=40):
    return((n-1)*(n-2)*two_stars(t,rho,NMAX) + (n-1)*edges(t)) 

def F(n,t,rho):
    p = p0(n,t,rho)
    return( np.array( [mean_degree(n,t) , mean_degree_sq(n,t,rho) ] )/(1-p) )

def J(n,t,rho,eps):
    f = F(n,t,rho)
    return( np.array([(F(n,t + eps,rho)-f)/eps, (F(n,t,rho+eps)-f)/eps]) )

def find_params_fixed_n(k,k2,n):
    mu = np.array([k,k2])
    x = np.array([ solve_t(k,n), solve_rho(k,k2,n) ] )
    for i in range(8):
        f = F(n,x[0],x[1])
        j = J(n,x[0],x[1],10**-8)
        x = x - np.dot(np.linalg.inv(j.astype(float).T),f-mu)
        print(f-mu)
    return(x)

def num_nonzero(k,k2,n):
    x = find_params_fixed_n(k,k2,n)
    return( (n-1)*(1-p0(n,x[0],x[1])) )

def find_n(k,k2,n_obs,n_low,n_high):
    low = num_nonzero(k,k2,n_low) - n_obs
    high = num_nonzero(k,k2,n_high) - n_obs
    while (n_high-n_low)>1:
        print((n_low,n_high))
        n_mid = int( (n_high+n_low)/2 )
        mid = num_nonzero(k,k2,n_mid) - n_obs
        if mid == 0:
            return n_mid,n_mid
        if mid > 0:
            high = mid
            n_high = n_mid
        else:
            low = mid
            n_low = n_mid
    return(n_low,n_high)
        


