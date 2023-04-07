# Common imports
import os
from math import exp, sqrt,pi
import numpy as np
import matplotlib.pyplot as plt
from decimal import *
from mpi4py import MPI
rng = np.random.default_rng()

# Where to save the figures and data files
DATA_ID = "project1f/Results"
if not os.path.exists(DATA_ID):
    os.makedirs(DATA_ID)
def data_path(dat_id):
    return os.path.join(DATA_ID, dat_id)




'''Here we set up the differents functions'''
# Trial wave function
def WaveFunction(r,alpha,N,d):
    S=np.linalg.norm(r)
    return (sqrt((8*alpha)/pi))**(N*d/2)*exp(-alpha*S**2)


# Analytic expression of the local energy
def LocalEnergy(r,alpha,N,d):
    S=np.linalg.norm(r)
    return (0.5-2*alpha*alpha)*(S**2)+N*alpha*d

#Quantum force
def QuantumForce(r,alpha):
    return -4*alpha*r

#Green function
def GreenFunction(x,y,Dt,D,alpha):
    X=np.linalg.norm(x)**2
    Y=np.linalg.norm(y)**2
    P=(-Y+X)*(4*alpha*D*Dt*alpha-2*alpha)
    return P

#Derivative of the wave function
def derWave_func(r,alpha):
    return 1/alpha-np.linalg.norm(r)**2
    

##Variational Monte-Carlos programm for the gradient descent
def MCS_grad(alpha,MCC1):
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

    GreensFunction = 0.0
    #Energy for the local energy
    E=0.0
    #Derivatives for the derivatives part in the expression of the derivative of the average local energy
    Dw=0.0
    Dww=0.0
    B=rng.normal(0.0, 1.0, size=(N, Dim))*sqrt(TimeStep)
    R_old=B
    wold = WaveFunction(R_old,alpha,N,Dim)
    QFold=QuantumForce(R_old, alpha)
    ##start of of the MC cycles
    for j in range(MCC1):
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
        DeltaDer= derWave_func(R_old,alpha)
        #Local energy
        E+=DeltaE
        #Derivative
        Dw+=DeltaDer
        #Product derivative and local energy
        Dww+=DeltaE*DeltaDer
    E/=MCC1
    Dw/=MCC1
    Dww/=MCC1
    grad=2*(Dww-Dw*E)
    return E,grad

##Variational Monte-Carlos programm for the final round
def MCS_last(alpha,MCC2):
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
    GreensFunction = 0.0
    B=rng.normal(0.0, 1.0, size=(N, Dim))*sqrt(TimeStep)
    R_old=B
    wold = WaveFunction(R_old,alpha,N,Dim)
    QFold=QuantumForce(R_old, alpha)
    ##start of of the MC cycles
    for j in range(MCC2):
        GreensFunction = 0.0
        A=rng.normal(0.0, 1.0, size=(N, Dim))*sqrt(TimeStep)
        R_new=R_old+A+ QFold*TimeStep*D
        GreensFunction = GreenFunction(R_old,R_new,TimeStep,D,alpha)
        wnew= WaveFunction(R_new,alpha,N,Dim)
        Qfnew=QuantumForce(R_new, alpha)
        GreensFunction=exp(GreensFunction)
        if np.random.random() <= (wnew**2 / wold**2)*GreensFunction:
            R_old=R_new
            QFold=Qfnew
            wold=wnew
        DeltaE = LocalEnergy(R_old,alpha,N,Dim)
        outfile.write('%f \n' %(DeltaE))        
    return


'''Here the programm really begins'''
##About the system
N=3
Dim=2

##For the code

MCC1=1000###Number of MCC for the gradient descent



alpha=0.4
Energy = 0
EDerivative = np.zeros((2), np.double)
eta = 0.01
Niterations = 50
delta=0.1

epsilon=0.001
memory=0
eff=0
#


for iter in range(Niterations):
    eff+=1
    Energy, EDerivative = MCS_grad(alpha,MCC1)
    alphagradient = EDerivative 
    alpha -= eta*alphagradient+delta*(alphagradient-memory)
    if abs(alpha-memory)<=epsilon:
        break
    else:
        memory=alpha

# print('ɑ='+str(alpha))
# print('efficiency='+str(eff))
# print(Energy)
'''Last VMC'''

##Parallelization
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()
MCC2=2**20##Number of MCC for the final analysis
outfile = open(data_path(f"VMCHarmonic{rank}.dat"),'w')


if rank == 0:
    ave=MCC2//nprocs
    data = [ave for p in range(nprocs)]
else:
    data = None

data = comm.scatter(data, root=0)


outfile = open(data_path(f"VMCHarmonic{rank}.dat"),'w')
MCS_last(alpha,data)

##Command to execute
#
##mpiexec -n 4 python project1f\code.py






