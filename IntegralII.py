"""
 --- Investigate integral with values for s, alpha, beta. ---
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.special import j1

def integrand(p,s,alpha,beta):
    q=(alpha*p)**2
    T1=( (p**2)*(1+q) )/beta
    A1=-s*(p**2)*( 1+(2*q)/3 )
    A2=( beta*p )/( 1+q )
    return T1*np.exp(A1)*j1(A2)

"""
 --- s>=0, alpha>=0, beta=2rsin(omega/2)   ---
 --- r>=0, omega=tau-tau'... -inf<beta<inf ---
"""

s=np.arange(1e-6,2,1e-4)
alpha=[1e-15,1e-9,1e-3]
beta=np.arange(1e-6,2,1e-4)                  # Extend to negative numbers after.

ax=plt.figure().add_subplot(projection='3d')
ax.set_xlabel('s')
ax.set_ylabel(r'$\beta$')
ax.set_zlabel('Z')

for a in alpha:
    Z,err=zip(*[quad(integrand,0,np.inf,(s[i],a,beta[i]),limit=200) for i in range(len(s))])
    ax.plot(s,beta,Z,label=fr'$\alpha$ = {a}')

ax.legend()
plt.show()