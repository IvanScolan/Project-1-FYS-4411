import matplotlib.pyplot as plt
import numpy as np
import random as rd
'''This code allows us to visualize the results'''


alpha_a,E,Eexact_a,T_a=np.loadtxt('project1b/CPU time analysis/Results/analytical_exp.dat', delimiter=' ', usecols =(0,1,2,3), unpack=True)
alpha_n,E_num,Eexact_n,T_n=np.loadtxt('project1b/CPU time analysis/Results/numerical_exp.dat', delimiter=' ', usecols =(0,1,2,3), unpack=True)


plt.plot(alpha_a,T_a,'-o',label='VMC using analytical expression of the local energy')
plt.plot(alpha_n,T_n,'-r+',label='VMC using a numerical calculation of the local energy ')
plt.ylabel('CPU time',fontsize='22')
plt.legend(fontsize=18)
plt.tick_params(labelsize=18)
plt.xlabel('N' ,fontsize='22')
plt.show()