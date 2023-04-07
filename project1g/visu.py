import matplotlib.pyplot as plt
import numpy as np
from numpy import log2, zeros, mean, var, sum, arange, array, cumsum
from numpy.linalg import inv
from math import sqrt


E0=np.loadtxt('project1g/Results/VMCHarmonic0.dat', delimiter=' ', usecols =(0), unpack=True)
E1=np.loadtxt('project1g/Results/VMCHarmonic1.dat', delimiter=' ', usecols =(0), unpack=True)
E2=np.loadtxt('project1g/Results/VMCHarmonic2.dat', delimiter=' ', usecols =(0), unpack=True)
E3=np.loadtxt('project1g/Results/VMCHarmonic3.dat', delimiter=' ', usecols =(0), unpack=True)


E=np.zeros(len(E0)+len(E1)+len(E2)+len(E3))

# E[:len(E0)]=E0
# E[len(E0):len(E1)+len(E0)]=E1
# E[len(E0)+len(E1):len(E0)+len(E1)+len(E2)]=E2
# E[len(E0)+len(E1)+len(E2):len(E0)+len(E1)+len(E2)+len(E3)]=E3


def block(x):
    # preliminaries
    n = len(x)
    d = int(log2(n))
    s, gamma = zeros(d), zeros(d)
    mu = mean(x)
    # estimate the auto-covariance and variances
    # for each blocking transformation
    for i in arange(0,d):
        n = len(x)
        # estimate autocovariance of x
        gamma[i] = (n)**(-1)*sum( (x[0:(n-1)]-mu)*(x[1:n]-mu) )
        # estimate variance of x
        s[i] = var(x)
        # perform blocking transformation
        x = 0.5*(x[0::2] + x[1::2])
    # generate the test observator M_k from the theorem
    M = (cumsum( ((gamma/s)**2*2**arange(1,d+1)[::-1])[::-1] ) )[::-1]
    # we need a list of magic numbers
    q =array([6.634897,  9.210340,  11.344867, 13.276704, 15.086272, 
              16.811894, 18.475307, 20.090235, 21.665994, 23.209251,
              24.724970, 26.216967, 27.688250, 29.141238, 30.577914, 
              31.999927, 33.408664, 34.805306, 36.190869, 37.566235,
              38.932173, 40.289360, 41.638398, 42.979820, 44.314105, 
              45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

    # use magic to determine when we should have stopped blocking
    for k in arange(0,d):
        if(M[k] < q[k]):
            break
    if (k >= d-1):
        print("Warning: Use more data")
    return mu, s[k]/2**(d-k)


(mean, var) = block(E)
std = sqrt(var)

import pandas as pd
from pandas import DataFrame
data ={'Mean':[mean], 'STDev':[std]}
frame = pd.DataFrame(data,index=['Values'])
print(frame)


