import matplotlib.pyplot as plt
import numpy as np
from numpy import log2, zeros, mean, var, sum, arange, array, cumsum
from numpy.linalg import inv
from math import sqrt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from scipy import stats

x,M0=np.loadtxt('project1h/Results/VMCHarmonic.dat', delimiter=' ', usecols =(0,1), unpack=True)
x1,M1=np.loadtxt('project1h/Results/VMCHarmonic1.dat', delimiter=' ', usecols =(0,1), unpack=True)


plt.plot(x,M0,label='With the Jaystrow factor')
plt.plot(x1,M1,label='Without the Jaystrow factor')
plt.legend(fontsize=18)
plt.tick_params(labelsize=18)
plt.xlabel('r',fontsize=18)
plt.ylabel('F(r)',fontsize=18)
plt.show()
