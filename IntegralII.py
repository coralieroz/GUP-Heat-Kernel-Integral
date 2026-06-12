"""
 --- Investigate integral with values for s, alpha, beta. ---
"""

import numpy as np
import matplotlib.pyplot as plt
import warnings
from scipy.integrate import quad
from scipy.special import j1

def integrand(p,s,alpha,beta):
    q=(alpha*p)**2
    T1=( (p**2)*(1+q) )/beta
    A1=-s*(p**2)*( 1+(2*q)/3 )
    A2=( beta*p )/( 1+q )
    return T1*np.exp(A1)*j1(A2)

def Z_grid(s,a,beta):
    Z=np.empty((s.size,beta.size))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")     # Ignores convergence warnings when doing integrals.
        for i,e1 in enumerate(s):
            for j,e2 in enumerate(beta):
                Z[i,j]=quad(integrand,0,np.inf,(e1,a,e2),limit=200)[0]
        return Z

"""
 --- s>=0, alpha>=0, beta=2rsin(omega/2)   ---
 --- r>=0, omega=tau-tau'... -inf<beta<inf ---
"""

n=100                                           # Number of points
s=np.linspace(1e-6,100,n)
alpha=[1e-15,1e-9,1e-3]
beta=np.linspace(-100,100,n)                    # Extend to negative numbers after.

S,B=np.meshgrid(s,beta)
fig = plt.figure(figsize=(15, 5))
for k, a in enumerate(alpha, start=1):
    Z  = Z_grid(s,a,beta)
    ax = fig.add_subplot(1, len(alpha), k, projection='3d')
    surf = ax.plot_surface(S, B, Z, cmap='viridis', linewidth=0, antialiased=True)
    ax.set_xlabel('s')
    ax.set_ylabel(r'$\beta$')
    ax.set_zlabel('Z')
    ax.set_title(fr'$\alpha$ = {a}')
    fig.colorbar(surf, ax=ax, shrink=0.6, pad=0.1)

fig.tight_layout()
plt.show()