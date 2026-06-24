"""
 --- Investigate the sharp positive peak seen in IntegralNumeric.py. ---

 The peak sits at the smallest s value of the grid (s_min = 1e-6), where the
 integrand decays very slowly with p (the exponential damping exp(-s p^2) is
 almost switched off), so the p-integral is highly oscillatory and only
 conditionally convergent. quad() with a finite 'limit' may not converge here.

 Two checks distinguish a boundary-cutoff artifact from an inherent feature:

  1. Lower-bound test: recompute Z on a grid whose s lower bound is moved
     away from 0. If the peak height collapses and its location moves off
     the s-boundary, the peak height was set by the cutoff s_min, not by Z
     itself.

  2. Resolution test: recompute Z on the *same* (s_min=1e-6) domain at a
     different grid resolution. If the peak value/location at s=s_min
     changes wildly (even in sign) with resolution, the values there are
     unconverged numerical noise rather than a smooth, well-defined feature.
"""

import numpy as np
import matplotlib.pyplot as plt
import warnings
from scipy.integrate import quad
from scipy.special import j1

# --- integrand: beta=0-safe rewrite (see IntegralNumeric.py / plot_surface_stable_integrand.py) ---
def integrand(p,s,alpha,beta):
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


def peak_info(s,beta,Z):
    """Largest-magnitude (s,beta,Z) on the grid, and whether it sits on the lower s-boundary."""
    idx = np.unravel_index(np.argmax(np.abs(Z)), Z.shape)
    return s[idx[0]], beta[idx[1]], Z[idx], idx[0]==0


def plot_Z(ax, s, beta, Z, title):
    S,B = np.meshgrid(s,beta)
    surf = ax.plot_surface(S, B, Z.T, cmap='viridis', linewidth=0, antialiased=True)
    sp,bp,zp,_ = peak_info(s,beta,Z)
    ax.scatter([sp],[bp],[zp],color='r',s=40,depthshade=False)
    ax.set_xlabel('s')
    ax.set_ylabel(r'$\beta$')
    ax.set_zlabel('Z')
    ax.set_title(title)
    return surf


alpha = 1e-15     # alpha value for which the peak is most prominent

"""
 --- Part 1: two grids with different lower bounds on s -----------------
 --- (same beta range and same resolution n1)                          ---
"""
n1 = 50
beta1 = np.linspace(-100,100,n1)
s_min_values = [1e-6, 1e-1]

Z_lb = []
print("--- Part 1: dependence on the lower bound of the s-grid ---")
for s_min in s_min_values:
    s = np.linspace(s_min,100,n1)
    Z = Z_grid(s,alpha,beta1)
    Z_lb.append((s,Z))
    sp,bp,zp,on_bdy = peak_info(s,beta1,Z)
    print(f"  s_min={s_min:g}: peak Z={zp:.4g} at s={sp:.4g}, beta={bp:.4g}  (on s-boundary: {on_bdy})")

"""
 --- Part 2: two grids with different resolutions -----------------------
 --- (same bounds as IntegralNumeric.py: s_min=1e-6, beta in [-100,100]) ---
"""
n_values = [30, 80]

Z_res = []
print("\n--- Part 2: dependence on grid resolution (s_min fixed at 1e-6) ---")
for n in n_values:
    s = np.linspace(1e-6,100,n)
    beta = np.linspace(-100,100,n)
    Z = Z_grid(s,alpha,beta)
    Z_res.append((s,beta,Z))
    sp,bp,zp,on_bdy = peak_info(s,beta,Z)
    print(f"  n={n}: peak Z={zp:.4g} at s={sp:.4g}, beta={bp:.4g}  (on s-boundary: {on_bdy})")

# --- Plots ---
fig = plt.figure(figsize=(13,10))

ax1 = fig.add_subplot(2,2,1,projection='3d')
surf1 = plot_Z(ax1, Z_lb[0][0], beta1, Z_lb[0][1], fr'$\alpha$={alpha:g}, $s_{{min}}$={s_min_values[0]:g}, n={n1}')
fig.colorbar(surf1, ax=ax1, shrink=0.6, pad=0.1)

ax2 = fig.add_subplot(2,2,2,projection='3d')
surf2 = plot_Z(ax2, Z_lb[1][0], beta1, Z_lb[1][1], fr'$\alpha$={alpha:g}, $s_{{min}}$={s_min_values[1]:g}, n={n1}')
fig.colorbar(surf2, ax=ax2, shrink=0.6, pad=0.1)

ax3 = fig.add_subplot(2,2,3,projection='3d')
surf3 = plot_Z(ax3, Z_res[0][0], Z_res[0][1], Z_res[0][2], fr'$\alpha$={alpha:g}, $s_{{min}}$=1e-6, n={n_values[0]}')
fig.colorbar(surf3, ax=ax3, shrink=0.6, pad=0.1)

ax4 = fig.add_subplot(2,2,4,projection='3d')
surf4 = plot_Z(ax4, Z_res[1][0], Z_res[1][1], Z_res[1][2], fr'$\alpha$={alpha:g}, $s_{{min}}$=1e-6, n={n_values[1]}')
fig.colorbar(surf4, ax=ax4, shrink=0.6, pad=0.1)

fig.suptitle('Peak investigation: lower-bound dependence (top) vs resolution dependence (bottom)')
fig.tight_layout()
plt.show()
