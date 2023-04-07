import matplotlib.pyplot as plt
import numpy as np
import random as rd

alpha,E,Eexact,N,Dim=np.loadtxt('project1c/VMC_E/VMCHarmonic.dat', delimiter=' ', usecols =(0,1,2,3,4), unpack=True)



# plt.subplot(2,1,1)
plt.plot(alpha,E,'-o',label='VMC using analytical expression of the local energy')
plt.plot(alpha,Eexact, label='Analytical solution of $E_L(ɑ)$')
plt.ylabel('Dimensionless energy',fontsize='22')
plt.legend(fontsize=18)
plt.tick_params(labelsize=18)
plt.title('Number of particles :'+str(int(N[0]))+', Dimensions :'+str(int(Dim[0])),fontsize='22')
plt.xlabel('ɑ' ,fontsize='22')

plt.show()