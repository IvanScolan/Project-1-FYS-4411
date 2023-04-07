# Common imports
import os
from math import exp, sqrt, pi
from random import random, seed
import numpy as np
import matplotlib.pyplot as plt
from decimal import *

# Where to save the figures and data files
DATA_ID = "project1b/Complete code for N particles and D dimension/Results"
if not os.path.exists(DATA_ID):
    os.makedirs(DATA_ID)

def data_path(dat_id):
    return os.path.join(DATA_ID, dat_id)

outfile = open(data_path("VMCHarmonic.dat"),'w')


'''Here we set up the differents functions'''
# Trial wave function
def WaveFunction(r,alpha,N,d):
    S=0
    for i in range(N):
        for j in range(d):
            S += r[i][j]*r[i][j]
    return (sqrt((8*alpha)/pi))**(N*d)*exp(-alpha*S)


# Analytic expression of the local energy
def LocalEnergy(r,alpha,N,d):
    R=0
    for i in range(N):
        for j in range(d):
            R+=((0.5-2*alpha*alpha)*r[i][j]*r[i][j])
    return R+N*alpha*d


##Numeric evaluation of the local energy
def Der_kin(psi,r,alpha,N,d):
    delta=0.001
    S=0
    rnew1=r
    rnew2=r
    psip=0.0
    psim=0.0

    for i in range(N):
        for j in range(d):
            g=r[i][j]
            rnew1[i][j]=delta+g
            psip=WaveFunction(rnew1,alpha,N,d)
            rnew2[i][j]=-delta+g
            psim=WaveFunction(rnew2,alpha,N,d)
            S+=1/(delta*delta)*(psip-2*psi+psim)
            rnew1[i][j]=g
            rnew2[i][j]=g
    return -0.5/WaveFunction(r,alpha,N,d)*S

##Norm**2 of a vector
def norm(r,N,d):
    S=0
    for i in range(N):
        for j in range(d):
            S += r[i][j]*r[i][j]
    return S

##Variational Monte-Carlos programm
def MCS():
    #Mont-Carlos cycles
    MCC=100000
    #Stepsize
    step=0.6
    #Initial value of alpha
    alpha=0.45
    #Initial position
    R_old=np.zeros((N,Dim),np.double)
    R_new=np.zeros((N,Dim),np.double)
    
    ##Start of the MC sampling
    for t in range(var_of_alpha):
        dalpha=+0.007
        alpha+=dalpha
        E=0.0
        Ee=0.0
        Enum_tri=0.0
        for i in range(N):
            for j in range(Dim):
                R_old[i][j]=step*(random()-0.5)
        wold = WaveFunction(R_old,alpha,N,Dim)
        ##start of of the MC cycles
        for j in range(MCC):
            
            for i in range(N):
                for d in range(Dim):
                    R_new[i][d]=R_old[i][d]+step*(random()-0.5)
            wnew= WaveFunction(R_new,alpha,N,Dim) 
            if random() <= wnew**2 / wold**2:
                for i in range(N):
                    for d in range(Dim):
                        R_old[i][d]=R_new[i][d]
                wold=wnew
            DeltaE = LocalEnergy(R_old,alpha,N,Dim)#Here we evaluate the local energy analyticaly
            DeltaEnum=Der_kin(wold,R_old,alpha,N,Dim)#Here we evaluate the local energy numericaly
            nor=norm(R_old,N,Dim)
            E+=DeltaE
            Enum_tri+=DeltaEnum +0.5*nor
        E/=MCC
        Eexact=N*Dim*(1/2*alpha+1/(8*alpha))#Computation of the exact value of the local energy for the given alpha
        Enum_tri/=MCC
        outfile.write('%f %f %f %f %f %f\n' %(alpha,E,Eexact,Enum_tri,N,Dim))
    return 


'''Here the programm really begins'''

var_of_alpha=20
N=1
Dim=1

MCS()



