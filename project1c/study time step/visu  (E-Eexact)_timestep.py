import matplotlib.pyplot as plt
import numpy as np
import random as rd

tm1,ecart1=np.loadtxt('project1c/study time step/Results/VMCHarmonic1.dat', delimiter=' ', usecols =(0,1), unpack=True)








plt.plot(tm1,ecart1,'-o')
plt.title('$N=2000$ and $Dim=1$',fontsize='22')
plt.tick_params(labelsize=18)
plt.ylabel('$f(\delta_t)$',fontsize='22')
plt.xlabel('$\delta_t$',fontsize='22')
plt.show()