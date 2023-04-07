import matplotlib.pyplot as plt
import numpy as np
import random as rd

tm1,ecart1=np.loadtxt('project1c/study time step/Results/VMCHarmonic1.dat', delimiter=' ', usecols =(0,1), unpack=True)

tm_zoom,ecart_zoom=np.loadtxt('project1c/study time step/Results/VMCHarmonic2.dat', delimiter=' ', usecols =(0,1), unpack=True)
# def find_value(L,l):
#     for i in range(len(L)):
#         if l<=L[i]:
#             return i
#     return

# print(find_value(tm,1))
# tm_zoom=tm[:find_value(tm,1)]
# ecart_zoom=ecart[:find_value(tm,1)]







plt.subplot(2,1,1)
plt.plot(tm1,ecart1,'-o')

plt.title('$N=2000$ and $Dim=1$',fontsize='22')
plt.tick_params(labelsize=18)
plt.ylabel('$f(\delta_t)$',fontsize='22')
plt.subplot(2,1,2)
plt.plot(tm_zoom,ecart_zoom,'-o')

plt.ylabel('$f(\delta_t)$',fontsize='22')
plt.tick_params(labelsize=18)
plt.xlabel('$\delta_t$',fontsize='22')
plt.show()