import matplotlib.pyplot as plt
import numpy as np
import random as rd


alpha,E,N,D=np.loadtxt('project1b/Fast code/Results/VMCHarmonic.dat', delimiter=' ', usecols =(0,1,2,3), unpack=True)

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

Ef=turbu(E,0)
ground_state,indice=Ef.min(),alpha[Ef.argmin()]
print('The ground state energy is : '+str(ground_state))
print('It is obtained for  ɑ='+str(indice))

Eexact=N[0]*D[0]*(1/2*0.5+1/(8*0.5))
print('The exact value for the energy is '+str(Eexact))

plt.plot(alpha,Ef,'-o',label='$T(<E_L>)$ ')
plt.ylabel('Dimensionless energy',fontsize='22')
plt.title('N='+str(N[0])+' and Dim='+str(D[0]),fontsize='22')
plt.legend(fontsize=18)
plt.tick_params(labelsize=18)
plt.xlabel('ɑ' ,fontsize='22')
plt.show()