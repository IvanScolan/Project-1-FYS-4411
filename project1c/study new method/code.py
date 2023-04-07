# Common imports
import os
from math import exp, sqrt,pi
from random import random, seed,  normalvariate
import numpy as np
import matplotlib.pyplot as plt
from decimal import *
rng = np.random.default_rng()

DATA_ID = "project1c/study new method/Results"
if not os.path.exists(DATA_ID):
    os.makedirs(DATA_ID)

def data_path(dat_id):
    return os.path.join(DATA_ID, dat_id)

outfile = open(data_path("VMCHarmonic.dat"),'w')




# Trial wave function for the Harmonic oscillator
def WaveFunction(r,alpha,N,d):
    S=np.linalg.norm(r)
    return (sqrt((8*alpha)/pi))**(N*d/2)*exp(-alpha*S**2)
# Local energy for the Harmonic oscillator
def LocalEnergy(r,alpha,N,d):
    S=np.linalg.norm(r)
    return (0.5-2*alpha*alpha)*(S**2)+N*alpha*d

#Quantum force
def QuantumForce(r,alpha):
    return -4*alpha*r

#Green function ratio

def GreenFunction(x,y,Dt,D,alpha):
    X=np.linalg.norm(x)**2
    Y=np.linalg.norm(y)**2
    P=(-Y+X)*(4*alpha*D*Dt*alpha-2*alpha)
    return P

def find_value(E):
    m=2**19
    i=0
    for e in range(len(E)):
        if E[e]<m:
            i=e 
            m=E[e]
    return i




def MCS_find_alpha_opt(TimeStep,alpha):
    #Mont-Carlos cycles
    MCC=10
    #Initial value of alpha
    Ee=np.zeros(len(alpha))
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
    for t in range(len(alpha)):
        E=0.0
        B=rng.normal(0.0, 1.0, size=(N, Dim))*sqrt(TimeStep)
        R_old=B
        wold = WaveFunction(R_old,alpha[t],N,Dim)
        QFold=QuantumForce(R_old, alpha[t])
        ##start of of the MC cycles
        for j in range(MCC):
            GreensFunction = 0.0
            A=rng.normal(0.0, 1.0, size=(N, Dim))*sqrt(TimeStep)
            R_new=R_old+A+ QFold*TimeStep*D
            GreensFunction= GreenFunction(R_old,R_new,TimeStep,D,alpha[t])
            wnew= WaveFunction(R_new,alpha[t],N,Dim)
            Qfnew=QuantumForce(R_new,alpha[t])
            GreensFunction=exp(GreensFunction)
            if np.random.random() <= (wnew**2 / wold**2)*GreensFunction:
                R_old=R_new
                QFold=Qfnew
                wold=wnew
            DeltaE = LocalEnergy(R_old,alpha[t],N,Dim)
            E+=DeltaE
        E/=MCC
        Ee[t]=E
    return Ee



def MCS_last(alpha):
    #Mont-Carlos cycles
    MCC=2**19
    #Value of D for importance sampling
    D=0.5
    TimeStep=0.43
    #Initial position
    R_old=np.zeros((N,Dim),np.double)
    R_new=np.zeros((N,Dim),np.double)
    #Quantum force
    QFold=np.zeros((N,Dim),np.double)
    Qfnew=np.zeros((N,Dim),np.double)

    GreensFunction = 0.0
    ##Start of the MC sampling
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
        Qfnew=QuantumForce(R_new,alpha)
        GreensFunction=exp(GreensFunction)
        if np.random.random() <= (wnew**2 / wold**2)*GreensFunction:
            R_old=R_new
            QFold=Qfnew
            wold=wnew
        DeltaE = LocalEnergy(R_old,alpha,N,Dim)
        E+=DeltaE
    E/=MCC
    return E


Alpha=np.linspace(0.45,0.55,100)
N=500
Dim=1

tm=np.linspace(0.42,0.45,100)

Moy=np.zeros(len(Alpha))
Var1=np.zeros(len(Alpha))

##Here start the code to compute the average value and the variance
for t in range(len(tm)):
    E=MCS_find_alpha_opt(tm[t],Alpha)
    Moy+=E
    Var1+=E**2
Moy/=len(tm)
Var1/=len(tm)
Var=Var1-Moy**2


for i in range(len(Alpha)):
    outfile.write('%f %f %f %f  \n' %(Alpha[i],Var[i],N,Dim))


alpha_opt=Alpha[Var.argmin()]

print('Optimal value for α found is '+str(alpha_opt))
E=MCS_last(alpha_opt)
print(E)

