# Common imports
import os
# Where to save the figures and data files
DATA_ID = "project1e/Results"
if not os.path.exists(DATA_ID):
    os.makedirs(DATA_ID)



def data_path(dat_id):
    return os.path.join(DATA_ID, dat_id)

outfile = open(data_path("VMCHarmonic.dat"),'w')

from math import exp, sqrt,pi
import numpy as np
import matplotlib.pyplot as plt
from decimal import *
rng = np.random.default_rng(12345)

# Trial wave function for the Harmonic oscillator
def WaveFunction(r,alpha,N,d):
    S=np.linalg.norm(r)
    return exp(-alpha*S**2)

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


#Derivative of the wave-function with respect to alpha
def derWave_func(r,alpha,N,Dim):
    return -np.linalg.norm(r)**2


#VMC code used for the gradient descent
def MCS_grad(alpha):
    #Mont-Carlos cycles
    MCC=5000
    #Time step
    TimeStep=0.43
    #Value of D for importance sampling
    D=0.5
    #Initial position
    R_old=np.zeros((N,Dim),np.double)
    R_new=np.zeros((N,Dim),np.double)
    #Quantum force
    QFold=np.zeros((N,Dim),np.double)
    Qfnew=np.zeros((N,Dim),np.double)

    GreensFunction = 0.0
    #Energy for the local energy
    En=0.0
    #Derivatives for the derivatives part in the expression of the derivative of the average local energy
    Dw=0.0
    Dww=0.0
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
        DeltaDer= derWave_func(R_old,alpha,N,Dim)
        #Local energy
        En+=DeltaE
        #Derivative
        Dw+=DeltaDer
        #Product derivative and local energy
        Dww+=DeltaE*DeltaDer
    En/=MCC
    Dw/=MCC
    Dww/=MCC
    grad=2*(Dww-Dw*En)
    return En,grad

#VMC code used after finding the optimal alpha
def MCS_last(alpha,n):
    #Mont-Carlos cycles
    MCC=n
    #Time step
    TimeStep=0.01
    #Value of D for importance sampling
    D=0.5
    #Initial position
    R_old=np.zeros((N,Dim),np.double)
    R_new=np.zeros((N,Dim),np.double)
    #Quantum force
    QFold=np.zeros((N,Dim),np.double)
    QFnew=np.zeros((N,Dim),np.double)
    GreensFunction = 0.0
    E=0.0
    Bb=rng.normal(0.0, 1.0, size=(N, Dim))*sqrt(TimeStep)
    R_old=Bb
    wold = WaveFunction(R_old,alpha,N,Dim)
    QFold=QuantumForce(R_old, alpha)
    ##start of of the MC cycles
    for j in range(MCC):
        Aa=rng.normal(0.0, 1.0, size=(N, Dim))*sqrt(TimeStep)
        R_new=R_old+Aa+ QFold*TimeStep*D
        wnew= WaveFunction(R_new,alpha,N,Dim)
        QFnew=QuantumForce(R_new, alpha)
        GreensFunction=GreenFunction(R_old,R_new,TimeStep,D,alpha)
        GreensFunction=exp(GreensFunction)
        if rng.random() <= (wnew**2 / wold**2)*GreensFunction:
            R_old=R_new
            QFold=QFnew
            wold=wnew
        DeltaE = LocalEnergy(R_old,alpha,N,Dim)
        outfile.write(' {:.12f}\n' .format(DeltaE) )
    return

'''Declaration of variable'''
##About the system
N=500
Dim=1

alpha=0.45
Energy = 0
EDerivative =0
eta = 0.001
Niterations = 50
delta=0.00001

epsilon=0.000001
memory=0
eff=0

for iter in range(Niterations):
    eff+=1
    Energy, EDerivative = MCS_grad(alpha)
    alphagradient = EDerivative
    alpha = alpha- eta*alphagradient+delta*(alpha-memory)
    memory=alpha




print('ɑ='+str(alpha))
print('efficiency='+str(eff))
print(Energy)
n=2**19
MCS_last(0.4999999999942274,n)




