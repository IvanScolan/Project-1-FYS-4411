import matplotlib.pyplot as plt
import numpy as np
import random as rd

alpha,var,N,D=np.loadtxt('project1c/study new method/Results/VMCHarmonic.dat', delimiter=' ', usecols =(0,1,2,3), unpack=True)


def turbu(L,n):
    if n<=1:
        H=np.zeros(len(L))
        for i in range(1,len(L)-1):
            H[i]=abs(L[i-1]-L[i])+abs(L[i]-L[i+1])
        H[0]=H[1]
        H[len(L)-1]=H[len(L)-2]
        return turbu(H,n+1)
    else:
        return L




varf=turbu(var,0)
ground_state,indice=varf.min(),alpha[varf.argmin()]
print('The minimum is obtained for  ɑ='+str(indice))

# plt.subplot(2,1,1)
plt.plot(alpha,var,'-o',label='Var(ɑ) ')
plt.title(' N='+str(N[0])+' and Dim='+str(D[0]),fontsize='22')
plt.legend(fontsize=18)
plt.tick_params(labelsize=18)


# plt.subplot(2,1,2)
# plt.plot(alpha,varf,'-o',label='$T(Var)(ɑ)$ ')
# plt.legend(fontsize=18)
# plt.tick_params(labelsize=18)

plt.xlabel('ɑ' ,fontsize='22')
plt.show()


