# Common imports
import os
from math import exp, sqrt,pi
from random import random, seed,  normalvariate
import numpy as np
import matplotlib.pyplot as plt
from decimal import *
rng = np.random.default_rng()

DATA_ID = "project1c/study time step/Results"
if not os.path.exists(DATA_ID):
    os.makedirs(DATA_ID)

def data_path(dat_id):
    return os.path.join(DATA_ID, dat_id)

outfile = open(data_path("VMCHarmonic.dat"),'w')




# Trial wave function for the Harmonic oscillator in one dimension
def WaveFunction(r,alpha,N,d):
    S=np.linalg.norm(r)
    return (sqrt((8*alpha)/pi))**(N*d/2)*exp(-alpha*S**2)
# Local energy for the Harmonic oscillator in one dimension
def LocalEnergy(r,alpha,N,d):
    S=np.linalg.norm(r)
    return (0.5-2*alpha*alpha)*(S**2)+N*alpha*d

#Quantum force
def QuantumForce(r,alpha):
    return -4*alpha*r

#Green function ratio for 1D

def GreenFunction(x,y,Dt,D,alpha):
    X=np.linalg.norm(x)**2
    Y=np.linalg.norm(y)**2
    P=(-Y+X)*(4*alpha*D*Dt*alpha-2*alpha)
    return P



def MCS(TimeStep):
    #Mont-Carlos cycles
    MCC=100000
    #Initial value of alpha
    alpha=0.45
    #Value of D for importance sampling
    D=0.5
    #Initial position
    R_old=np.zeros((N,Dim),np.double)
    R_new=np.zeros((N,Dim),np.double)
    #Quantum force
    QFold=np.zeros((N,Dim),np.double)
    Qfnew=np.zeros((N,Dim),np.double)

    GreensFunction = 0.0
    ##Start of the MC sampling
    for t in range(var_of_alpha):
        alpha+=0.005
        E=0.0
        B=rng.normal(0.0, 1.0, size=(N, Dim))*sqrt(TimeStep)
        R_old=B
        wold = WaveFunction(R_old,alpha,N,Dim)
        QFold=QuantumForce(R_old, alpha)
        ##start of of the MC cycles
        for j in range(MCC):
            GreensFunction = 0.0
            A=rng.normal(0.0, 1.0, size=(N, Dim))*sqrt(TimeStep)
            R_new=R_old+A+ QFold*TimeStep*D
            GreensFunction= GreenFunction(R_old,R_new,TimeStep,D,alpha)
            wnew= WaveFunction(R_new,alpha,N,Dim)
            Qfnew=QuantumForce(R_new, alpha)
            GreensFunction=exp(GreensFunction)
            if np.random.random() <= (wnew**2 / wold**2)*GreensFunction:
                R_old=R_new
                QFold=Qfnew
                wold=wnew
            DeltaE = LocalEnergy(R_old,alpha,N,Dim)
            E+=DeltaE
        E/=MCC
        Eexact=N*Dim*(1/2*alpha+1/(8*alpha))
        outfile.write('%f %f %f %f \n' %(alpha,E,Eexact,TimeStep))
    return 


var_of_alpha=25
N=10
Dim=2

Timestep=[0.5,4]
for j in Timestep:
    MCS(j)

