import matplotlib.pyplot as plt
import numpy as np
import random as rd

alpha,E,Eexact,tm=np.loadtxt('project1c/study time step/Results/VMCHarmonic.dat', delimiter=' ', usecols =(0,1,2,3), unpack=True)

var_alpha=1
num_tm=1
for j in range( len(tm)-1):
    if tm[j]==tm[j+1]:
        var_alpha+=1
    elif j==len(tm)-2:
        break
    else:
        var_alpha=1
        num_tm+=1

alphav=np.zeros((num_tm,var_alpha))
Ev=np.zeros((num_tm,var_alpha))
Eexactv=np.zeros((num_tm,var_alpha))
var=np.zeros((num_tm,var_alpha))

for i in range(num_tm):
    for j in range(var_alpha):
        alphav[i][j]=alpha[i*var_alpha+j]
        Ev[i][j]=E[i*var_alpha+j]
        Eexactv[i][j]=Eexact[i*var_alpha+j]
        var[i][j]=abs(Eexact[i*var_alpha+j]-E[i*var_alpha+j])





plt.subplot(2,1,1)
for t in range(num_tm):
    plt.plot(alphav[t],Ev[t],'-o',label='$\delta_t$='+str(tm[t*var_alpha]))
plt.plot(alphav[0],Eexactv[0],label='Exact solution',linewidth=2)
plt.title('$<E_L>(α)$ for different values of time step',fontsize='22')
plt.ylabel('Dimensionless energy',fontsize='22')
plt.legend(fontsize=18)
plt.tick_params(labelsize=18)


plt.subplot(2,1,2)
for t in range(num_tm):
    plt.plot(alphav[t],var[t],'-o',label='$\delta_t$='+str(tm[t*var_alpha]))
plt.xlabel('$α$',fontsize='22')
plt.legend(fontsize=18)
plt.title('Absolute deviation of the computed energy from the exact value',fontsize='22')
plt.tick_params(labelsize=18)

plt.show()