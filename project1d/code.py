# Common imports
import os

DATA_ID = "project1d/Results"

if not os.path.exists(DATA_ID):
    os.makedirs(DATA_ID)

def data_path(dat_id):
    return os.path.join(DATA_ID, dat_id)


outfile = open(data_path("VMCHarmonic.dat"),'w')

from math import exp, sqrt,pi
from random import random, seed,  normalvariate
import numpy as np
import matplotlib.pyplot as plt
from decimal import *
# Trial wave function for the Harmonic oscillator
def WaveFunction(r,alpha,N,d):
    S=0
    for i in range(N):
        for j in range(d):
            S += r[i][j]*r[i][j]
    return (sqrt((8*alpha)/pi))**(N*d/2)*exp(-alpha*S)

# Local energy for the Harmonic oscillator
def LocalEnergy(r,alpha,N,d):
    R=0
    for i in range(N):
        for j in range(d):
            R+=((0.5-2*alpha*alpha)*r[i][j]*r[i][j])
    return R+N*alpha*d

#Quantum force
def QuantumForce(r,alpha,N,d):
    F=np.zeros((N,d),np.double)
    for i in range(N):
        for j in range(d):
            F[i][j]=-4*alpha*r[i][j]
    return F

#Green function ratio for 1D

def GreenFunction(x,y,Dt,D,alpha):
    P=(-y**2+x**2)*(4*alpha*D*Dt*alpha-2*alpha)
    # P=2*alpha*(y-x)*(2*alpha*D*Dt*(y-x)-y+x)
    return P

#Derivative of the wave-function with respect to alpha

def derWave_func(r,alpha,N,d):
    S=0
    for i in range(N):
        for j in range(d):
            S+=r[i][j]**2
    return 1/alpha-S
    


def MCS(alpha):
    #Mont-Carlos cycles
    MCC=10000
    #Time step
    TimeStep=0.05
    #Value of D for importance sampling
    D=0.5
    #Initial position
    R_old=np.zeros((N,Dim),np.double)
    R_new=np.zeros((N,Dim),np.double)
    #Quantum force
    QFold=np.zeros((N,Dim),np.double)
    Qfnew=np.zeros((N,Dim),np.double)

    seed()
    GreensFunction = 0.0
    #Energy for the local energy
    E=0.0
    #Derivatives for the derivatives part in the expression of the derivative of the average local energy
    Dw=0.0
    Dww=0.0
    for i in range(N):
        for k in range(Dim):
            R_old[i][k]=normalvariate(0.0,1.0)*sqrt(TimeStep)
    wold = WaveFunction(R_old,alpha,N,Dim)
    QFold=QuantumForce(R_old, alpha, N,Dim)
    ##start of of the MC cycles
    for j in range(MCC):
        for i in range(N):
            GreensFunction = 0.0
            for k in range(Dim):
                R_new[i][k]=R_old[i][k]+normalvariate(0.0,1.0)*sqrt(TimeStep)+ QFold[i][k]*TimeStep*D
                GreensFunction += GreenFunction(R_old[i][k],R_new[i][k],TimeStep,D,alpha)
            wnew= WaveFunction(R_new,alpha,N,Dim)
            Qfnew=QuantumForce(R_new, alpha, N,Dim)
            GreensFunction=exp(GreensFunction)
            if random() <= (wnew**2 / wold**2)*GreensFunction:
                for k in range(Dim):
                    R_old[i][k]=R_new[i][k]
                    QFold[i][k]=Qfnew[i][k]
                wold=wnew
        DeltaE = LocalEnergy(R_old,alpha,N,Dim)
        DeltaDer= derWave_func(R_old,alpha,N,Dim)
        #Local energy
        E+=DeltaE
        #Derivative
        Dw+=DeltaDer
        #Product derivative and local energy
        Dww+=DeltaE*DeltaDer
    E/=MCC
    Dw/=MCC
    Dww/=MCC
    grad=2*(Dww-Dw*E)
    return E,grad


var_of_alpha=20
N=500
Dim=2
alpha=0.4
Energy = 0
EDerivative = 0
eta = 0.01
Niterations = 50
delta=0.1

epsilon=0.001
memory=0
eff=0
#
for iter in range(Niterations):
    eff+=1
    Energy, EDerivative = MCS(alpha)
    alphagradient = EDerivative 
    alpha -= eta*alphagradient+delta*(alphagradient-memory)
    if abs(alpha-memory)<=epsilon:
        break
    else:
        memory=alpha

print('ɑ='+str(alpha))
print('efficiency='+str(eff))

print(Energy)