# Common imports
import os
from math import exp, sqrt, pi
from random import random, seed
import numpy as np
import matplotlib.pyplot as plt
from decimal import *
from time import process_time

"""This version of the code is an optimized one by using the vectorial operation of the numpy library"""

DATA_ID = "project1b/Fast code/Results"

if not os.path.exists(DATA_ID):
    os.makedirs(DATA_ID)
def data_path(dat_id):
    return os.path.join(DATA_ID, dat_id)

outfile = open(data_path("VMCHarmonic.dat"),'w')


'''Here we set up the differents functions'''


# Trial wave function for the Harmonic oscillator in one dimension
def WaveFunction(r,alpha,N,d):
    S=np.linalg.norm(r)
    return (sqrt((8*alpha)/pi))**(N*d/2)*exp(-alpha*S**2)
# Local energy for the Harmonic oscillator in one dimension
def LocalEnergy(r,alpha,N,d):
    S=np.linalg.norm(r)
    return (0.5-2*alpha*alpha)*(S**2)+N*alpha*d


def MCS():
    #Mont-Carlos cycles
    MCC=100000
    #Stepsize
    step=0.2
    #Initial value of alpha
    alpha=0.45
    #Initial position
    R_old=np.zeros((N,Dim),np.double)
    R_new=np.zeros((N,Dim),np.double)
    
    ##Start of the MC sampling
    for t in range(var_of_alpha):
        dalpha=+0.001
        alpha+=dalpha
        E=0.0
        B=(np.random.random((N,Dim))-0.5)*step
        R_old=B
        wold = WaveFunction(R_old,alpha,N,Dim)
        ##start of of the MC cycles
        for j in range(MCC):
            A=(np.random.random((N,Dim))-0.5)*step
            R_new=R_old+A
            wnew= WaveFunction(R_new,alpha,N,Dim) 
            if random() <= wnew**2 / wold**2:
                R_old=R_new
                wold=wnew
            DeltaE = LocalEnergy(R_old,alpha,N,Dim)
            E+=DeltaE
        E/=MCC
        outfile.write('%f %f %f %f\n' %(alpha,E,N,Dim))
    return 


'''Here the programm really begins'''

var_of_alpha=100
N=2
Dim=2


MCS()



