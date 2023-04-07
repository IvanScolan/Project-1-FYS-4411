import matplotlib.pyplot as plt
import numpy as np
import random as rd

alpha,E,Eexact,tm=np.loadtxt('project1c/VMC_E/VMCHarmonic.dat', delimiter=' ', usecols =(0,1,2,3), unpack=True)

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

for i in range(num_tm):
    for j in range(var_alpha):
        alphav[i][j]=alpha[i*num_tm+j]
        Ev[i][j]=E[i*num_tm+j]
        Eexactv[i][j]=Eexact[i*num_tm+j]



for t in range(num_tm):
    plt.plot(alphav[t],Ev[t],'-o',label='Time step='+str(tm[t*num_tm]))


plt.plot(alphav[0],Eexactv[0],label='analytical',linewidth=6)
    
  







plt.xlabel('$α$',fontsize='22')
plt.ylabel('Dimensionless energy',fontsize='22')
plt.legend(fontsize=18)
plt.tick_params(labelsize=18)


plt.show()