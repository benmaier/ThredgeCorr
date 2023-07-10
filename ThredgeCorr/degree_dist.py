from mpmath import mp
import numpy as np
import scipy.special

### precision: increase this if higher precision necessary
mp.dps = 500

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

def Gauss_Hermite_points(N,num_its=10):
    '''
    find the points and weights for Gauss Hermite quadrature
    scipy already has an algorithm to do this but it's not accurate
    so we use Netwon's method to improve
    '''
    X = list(scipy.special.he_roots(N)[0])
    for i in range(num_its):
        for root in range(N):
            X[root] = X[root] - hermite(N,X[root]) / (N * hermite(N-1,X[root]) )
    W = [(mp.factorial(N)*mp.sqrt(2*mp.pi))/mp.power(N*hermite(N-1,x),2) for x in X]
    return(X,W)

# ================================
X, W = Gauss_Hermite_points(40)
# ================================

def f(y,k,n,t,rho):
    P = Phi(y)
    ans = k*mp.log(1-P) + (n-1-k)*mp.log(P)
    ans = ans - ((t-mp.sqrt(1-rho)*y)**2)/(2*mp.mpf(rho))
    return ans

def dfdy(y,k,n,t,rho):
    return phi(y)* ( -k/(1-Phi(y)) + (n-1-k)/Phi(y) ) + mp.sqrt(1-rho)*(t-mp.sqrt(1-rho)*y)/rho

def d2fdy2(y,k,n,t,rho):
    ans = -y*phi(y)*( -k/(1-Phi(y)) + (n-1-k)/Phi(y) )
    ans = ans + phi(y)**2 * ( -k/((1-Phi(y))**2) - (n-1-k)/(Phi(y)**2) )
    ans = ans + (-(mp.mpf(1)-rho)/rho)
    return ans

def find_y0(k,n,t,rho):
    kk = min(max(1,k),n-2)
    #initial guess -- should be very close
    y0 = Phi_inv( (n-1.-kk)/(n-1.) )
    # a few Newton-method steps to get closer
    for i in range(10):
        y0 = y0 - dfdy(y0,k,n,t,rho)/d2fdy2(y0,k,n,t,rho)
    return y0


def ln_I(k,n,t,rho,X=X,W=W):
    N = len(X)
    y0 = find_y0(k,n,t,rho)
    fy0 = f(y0,k,n,t,rho)
    fppy0 = d2fdy2(y0,k,n,t,rho)
    r = lambda y: f(y,k,n,t,rho) - fy0 + ((y-y0)**2) *(abs(fppy0)/2.)
    ln_ans = mp.log(sum([W[i]*mp.exp(r(mp.sqrt(1./abs(fppy0))*X[i] +y0)) for i in range(N)]))
    ln_ans = ln_ans + fy0 + 0.5*mp.log(1-rho) - 0.5*mp.log(2*mp.pi*rho*abs(fppy0))
    return ln_ans

def log_binomial(n,k):
    return mp.loggamma(n+1) - mp.loggamma(k+1) - mp.loggamma(n-k+1)

def lnpk(n,t,rho,k_max,X=X,W=W):
    return np.array([log_binomial(n-1,k)+ln_I(k,n,t,rho,X,W) for k in range(k_max)])

def pk(n,t,rho,k_max,X=X,W=W):
    return np.array([mp.exp(log_binomial(n-1,k)+ln_I(k,n,t,rho,X,W)) for k in range(k_max)])

def log_likelihoods(n,t,rho,data,X=X,W=W):
    return np.array([log_binomial(n-1,k)+ln_I(k,n,t,rho,X,W) for k in data])

def pk_asymptotic(k, n_, t_, rho_):
    n = int(n_)
    t = float(t_)
    rho = float(rho_)

    a0 = 2.30753
    a1 = 0.27061
    b1 = 0.99229
    b2 = 0.04481

    x = 1-k/(n-1)
    low_x = x[x<=0.5]
    high_x = 1-x[x>0.5]
    ck = []

    for ix, x_ in enumerate([low_x,high_x]):
        s = np.sqrt(np.log(1/x_**2))
        this_phi = (a0+a1*s) / (1+b1*s + b2*s**2) - s
        prefac = 1 - 2*ix
        ck.append(prefac*this_phi)

    ck = np.concatenate(ck[::-1])

    fac = 1/(n-1)*np.sqrt((1-rho)/rho) * np.exp(-t**2/2/rho)
    return fac * np.exp((2*rho-1)/(2*rho) * ck**2 + t*np.sqrt(1-rho)/rho * ck)



if __name__=="__main__":

    n = 10
    k = np.arange(n)
    t = 1
    rho = 0.4
    logp = log_likelihoods(n, t, rho, k)

    mp.dps = 6

    import matplotlib.pyplot as pl

    pl.plot(logp)
    pl.show()

