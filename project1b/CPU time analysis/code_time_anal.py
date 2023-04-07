# Common imports
import os
from math import exp, sqrt, pi
from random import random, seed
import numpy as np
import matplotlib.pyplot as plt
from decimal import *
from time import process_time


DATA_ID = "project1b/CPU time analysis/Results"

if not os.path.exists(DATA_ID):
    os.makedirs(DATA_ID)


def data_path(dat_id):
    return os.path.join(DATA_ID, dat_id)

outfile1 = open(data_path("analytical_exp.dat"),'w')
outfile2 = open(data_path("numerical_exp.dat"),'w')



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


##Variational Monte-Carlos programm using the analytical expression of the local energy
def MCS_analytical():
    #Mont-Carlos cycles
    MCC=1000
    #Stepsize
    step=0.6
    #Initial value of alpha
    alpha=0.45
    #Initial position
    
    #initial energy
    seed()
    
    ##Start of the MC sampling
    for N in NL:
        R_old=np.zeros((N,Dim),np.double)
        R_new=np.zeros((N,Dim),np.double)
        dalpha=+0.02
        alpha+=dalpha
        E=0.0
        T=0
        for i in range(N):
            for j in range(Dim):
                R_old[i][j]=step*(random()-0.5)
        wold = WaveFunction(R_old,alpha,N,Dim)
        ##time analysis
        t1_start = process_time()
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
            DeltaE = LocalEnergy(R_old,alpha,N,Dim)
           
            E+=DeltaE
        E/=MCC
        Eexact=N*Dim*(1/2*alpha+1/(8*alpha))
        t1_stop = process_time()
        T=t1_stop-t1_start
        outfile1.write('%f %f %f %f \n' %(N,E,Eexact,T))
    return 

##Variational Monte-Carlos programm using the numeric calculation of the local energy
def MCS_num():
    #Mont-Carlos cycles
    MCC=1000
    #Stepsize
    step=0.6
    #Initial value of alpha
    alpha=0.45
    
    ##Start of the MC sampling
    for N in NL:
        R_old=np.zeros((N,Dim),np.double)
        R_new=np.zeros((N,Dim),np.double)
        Enum_tri=0.0
        T=0
        for i in range(N):
            for j in range(Dim):
                R_old[i][j]=step*(random()-0.5)
        wold = WaveFunction(R_old,alpha,N,Dim)
        ##time analysis
        t1_start = process_time()
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
            DeltaEnum=Der_kin(wold,R_old,alpha,N,Dim)
            nor=norm(R_old,N,Dim)
            Enum_tri+=DeltaEnum +0.5*nor
        Eexact=N*Dim*(1/2*alpha+1/(8*alpha))
        Enum_tri/=MCC
        t1_stop = process_time()
        T=t1_stop-t1_start
        outfile2.write('%f %f %f %f \n' %(N,Enum_tri,Eexact,T))
    return 

'''Here the programm really begins'''

#We compute the CPU time for different N
NL=[i for i in range(10)]
Dim=2

MCS_num()
MCS_analytical()







