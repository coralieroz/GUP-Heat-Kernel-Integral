"""
 --- Investigate integral with values for s, alpha, beta. ---
"""

import numpy as np
import matplotlib.pyplot as plt
import warnings
from scipy.integrate import quad
from scipy.special import j1

def integrand(p,s,alpha,beta):
    # T1*j1(A2) = p^2*(1+q)/beta * j1(A2) = p^3*j1(A2)/A2, which is finite at
    # beta=0 (limit p^3/2), unlike the naive T1=p^2*(1+q)/beta form.
    q=(alpha*p)**2
    A1=-s*(p**2)*( 1+(2*q)/3 )
    A2=( beta*p )/( 1+q )
    ratio = 0.5 if abs(A2) < 1e-8 else j1(A2)/A2
    return (p**3)*ratio*np.exp(A1)

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
s=np.linspace(0.1,100,n)                        # s_min=0.1: below this, exp(-s p^2) barely
                                                 # damps the oscillatory j1 term and quad(limit=200)
                                                 # no longer converges (see PeakInvestigation.py).
alpha=[1e-15,1e-9,1e-3]
beta=np.linspace(-100,100,n)                    # Extend to negative numbers after.

S,B=np.meshgrid(s,beta,indexing='ij')   # 'ij' keeps S,B aligned with Z[s_idx, beta_idx]
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