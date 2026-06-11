import numpy as np
from scipy.integrate import quad
from scipy.special import j0, j1
import matplotlib.pyplot as plt

x=np.arange(0,10,0.001)
y1=[j1(i,out=None) for i in x]
y0=[j0(i,out=None) for i in x]

plt.plot(x,y1)
plt.plot(x,y0)
plt.show()