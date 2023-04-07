# Common imports
import os
# Where to save the figures and data files
DATA_ID = "project1h/Results"
if not os.path.exists(DATA_ID):
    os.makedirs(DATA_ID)



def data_path(dat_id):
    return os.path.join(DATA_ID, dat_id)

from math import exp, sqrt,pi
import numpy as np
import matplotlib.pyplot as plt
from decimal import *
rng = np.random.default_rng()
from time import process_time
from mpi4py import MPI
from scipy import integrate


def g(alpha,R):
    S=np.linalg.norm(R)
    return exp(-alpha*(S**2))

def f(a,ri,rj):
    n=np.linalg.norm(ri-rj)
    if a>= n :
        return 0
    else:
        return 1-a/n

# Trial wave function for the Harmonic oscillator with interraction adapted to a two body system
def WaveFunction(R2i,R2j):
    R=np.zeros((2,2))
    R[0]=R1
    R[1]=(R2i,R2j)
    return g(alpha,R)*f(a,R[0],R[1])


##Alows us to compute  a 6 dimensional integral, however it is clearly not efficient enough
def integN3(psi,r,x_max,x_min,N,d,alpha,a):
    nn=10
    h=(x_max-x_min)/nn
    R=np.ones((N,Dim))
    R[0]=r
    R[1:N]=R[1:N]*x_min
    S=0
    for l in range(nn):
        for l in range(nn):
            for l in range(nn):
                for l in range(nn):
                    for l in range(nn):
                        for l in range(nn):
                            R[1][0]+=h
                            S+=abs(psi(R,alpha,N,a))**2
                        R[1][1]+=h
                        R[1][0]=x_min
                    R[1][1]=x_min
                    R[1][2]+=h
                R[1][2]=x_min
                R[2][0]+=h
            R[2][0]=x_min
            R[2][1]+=h
        R[2][1]=x_min
        R[2][2]+=h
    return S*h**((N-1)*d)

##Alows us to compute  a 2 dimensional integral, however it is clearly not efficient enough
def integN2(psi,r,x_max,x_min,alpha,a):
    nn=100
    h=(x_max-x_min)/nn
    R=np.zeros((2,2))
    R[0]=r
    R[1]=x_min
    S=0
    for l in range(nn):
        R[1][1]+=h
        for l in range(nn):
            R[1][0]+=h
            S+=abs(psi(R,alpha,a))**2
        R[1][0]=x_min
    return S*h**2




'''Declaration of variable'''
##About the system
N=2
##Here we need to set up Dim=3
Dim=2
alpha=0.22
beta=2.82843
gamma=2.82843

##Vector used for vectrorize 
BETA=np.ones((N,Dim))
GAMMA=np.ones((N,Dim))
for i in range(N):
    BETA[i][2]=beta
    GAMMA[i][2]=gamma



"""Onebody density with Jastrow factor"""

outfile = open(data_path("VMCHarmonic.dat"),'w')
n=50
x_max=5
x_min=-x_max
x=np.linspace(0,x_max,n)


a=0.0043
R1=np.zeros(2)
W=[]
S=0##will alows us to normalize by approximating the integral of F(r)
for x_i in range(n):
    R1[0]=x[x_i]
    D=integrate.nquad(WaveFunction, [[x_min, x_max],[x_min, x_max]])[0]
    W.append(D)
    S+=D

for w in range(len(x)):
    outfile.write('%f %f\n' %(x[w],W[w]/S))

"""Onebody density without Jastrow factor"""


outfile1 = open(data_path("VMCHarmonic1.dat"),'w')
n=50

a=0
W=[]
S=0##will alows us to normalize by approximating the integral of F(r)
for x_i in range(n):
    R1[0]=x[x_i]
    D=integrate.nquad(WaveFunction, [[x_min, x_max],[x_min, x_max]])[0]
    W.append(D)
    S+=D

for w in range(len(x)):
    outfile1.write('%f %f\n' %(x[w],W[w]/S))