# Common imports
import os
# Where to save the figures and data files
DATA_ID = "project1g/Results"
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

def g(alpha,R):
    S=np.linalg.norm(R*BETA)
    return exp(-alpha*(S**2))

def f(a,ri,rj):
    n=np.linalg.norm(ri-rj)
    if a>= n :
        return 0
    else:
        return 1-a/n

# Trial wave function for the Harmonic oscillator with interraction
def WaveFunction(R,alpha,N):
    
    B=1
    for i in range(N):
        for j in range(i):
            n=np.linalg.norm(R[i]-R[j])
            if a>= n :
                return 0
            else:
                B*=(1-a/n)
    
    return g(alpha,R)*B

##Numeric evaluation of the local energy
def Kinetic_energy(psi,r,alpha,N,d):
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
            psip=WaveFunction(rnew1,alpha,N)
            rnew2[i][j]=-delta+g
            psim=WaveFunction(rnew2,alpha,N)
            S+=1/(delta*delta)*(psip-2*psi+psim)
            rnew1[i][j]=g
            rnew2[i][j]=g
    return -0.5/WaveFunction(r,alpha,N)*S

"""This is an attempt to implement the analytic expression for the local energy"""
# def Kinetic_energy(r,alpha,beta,N,a,d):
#     A=-4*alpha*N-2*alpha*beta*N+4*alpha**2*np.linalg.norm(r*BETA)**2
#     B=np.zeros((N,d))
#     for i in range(N):
#         for l in range(i):
#             r_il=np.linalg.norm(r[i]-r[l])
#             B[i]+=(r[i]-r[l])*a/(r_il**3-a*r_il**2)
#         for l in range(i+1,N):
#             r_il=np.linalg.norm(r[i]-r[l])
#             B[i]+=(r[i]-r[l])*a/(r_il**3-a*r_il**2)
#     C=-4*alpha*np.linalg.norm(r*BETA*B)
#     D=0

#     for k in range(N):
#         for i in range(N):
#             if k!=i:
#                 r_ki=np.linalg.norm(r[k]-r[i])
#                 for j in range(N):
#                     if k!=j:
#                         r_kj=np.linalg.norm(r[k]-r[j])
#                         D+=np.linalg.norm((r[k]-r[i])*(r[k]-r[j]))*a**2/((r_ki**3-a*r_ki**2)*(r_kj**3-a*r_kj**2))
#     E=0
#     for i in range(N):
#         for l in range(i):
#             r_il=np.linalg.norm(r[i]-r[l])
#             E+=-2*a*(2*r_il-a)/(r_il**2-r_il*a)**2+2*a/(r_il**3-a*r_il**2)
#         for l in range(i+1,N):
#             r_il=np.linalg.norm(r[i]-r[l])
#             E+=-2*a*(2*r_il-a)/(r_il**2-r_il*a)**2+2*a/(r_il**3-a*r_il**2)
#     return -0.5*(A+C+D+E)
"""#######"""


def pot_non_int(R):
    S=np.linalg.norm(R*GAMMA)
    return 0.5*S**2




#Green function ratio

def GreenFunction(R_old,R_new,TimeStep,D,QFold,QFnew):
    S=0.5*(QFold+QFnew)*(D*TimeStep*0.5*(QFold-QFnew)-R_new+R_old)
    Sp=np.linalg.norm(S)
    return Sp

#Derivative of the wave function with respect to alpha

def derWave_func(r,alpha):
    S=np.linalg.norm(r*BETA)
    return -alpha*S**2
    

#Quantum force
def QuantumForce(r,alpha,N,d):
    A=-4*alpha*r*BETA
    B=np.zeros((N,d))
    for i in range(N):
        for l in range(i):
            r_il=np.linalg.norm(r[i]-r[l])
            B[i]+=(r[i]-r[l])*a/(r_il**3-a*r_il**2)
        for l in range(i+1,N):
            r_il=np.linalg.norm(r[i]-r[l])
            B[i]+=(r[i]-r[l])*a/(r_il**3-a*r_il**2)
    return A+2*B


##Variational Monte-Carlos programm for the gradient descent
def MCS_grad(alpha):
    #Mont-Carlos cycles
    MCC=1000
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
    E=0.0
    #Derivatives for the derivatives part in the expression of the derivative of the average local energy
    Dw=0.0
    Dww=0.0
    B=rng.normal(0.0, 1.0, size=(N, Dim))*sqrt(TimeStep)
    R_old=B
    wold = WaveFunction(R_old,alpha,N)
    QFold=QuantumForce(R_old,alpha,N,Dim)
    ##start of of the MC cycles
    for j in range(MCC):
        GreensFunction = 0.0
        A=rng.normal(0.0, 1.0, size=(N, Dim))*sqrt(TimeStep)
        R_new=R_old+A+ QFold*TimeStep*D
        wnew= WaveFunction(R_new,alpha,N)
        Qfnew=QuantumForce(R_new,alpha,N,Dim)
        GreensFunction= GreenFunction(R_old,R_new,TimeStep,D,QFold,Qfnew)
        GreensFunction=exp(GreensFunction)
        if np.random.random() <= (wnew**2 / wold**2)*GreensFunction:
            R_old=R_new
            QFold=Qfnew
            wold=wnew
        DeltaE = Kinetic_energy(wold,R_old,alpha,N,Dim)+pot_non_int(R_old)
        DeltaDer= derWave_func(R_old,alpha)
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
    return grad

##Variational Monte-Carlos programm for the final round
def MCS_last(alpha,n):
    #Mont-Carlos cycles
    MCC=n
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
    B=rng.normal(0.0, 1.0, size=(N, Dim))*sqrt(TimeStep)
    R_old=B
    wold = WaveFunction(R_old,alpha,N)
    QFold=QuantumForce(R_old,alpha,N,Dim)
    ##start of of the MC cycles
    for j in range(MCC):
        A=rng.normal(0.0, 1.0, size=(N, Dim))*sqrt(TimeStep)
        R_new=R_old+A+ QFold*TimeStep*D
        wnew= WaveFunction(R_new,alpha,N)
        Qfnew=QuantumForce(R_new,alpha,N,Dim)
        GreensFunction = GreenFunction(R_old,R_new,TimeStep,D,QFold,Qfnew)
        GreensFunction=exp(GreensFunction)
        if np.random.random() <= (wnew**2 / wold**2)*GreensFunction:
            R_old=R_new
            QFold=Qfnew
            wold=wnew
        DeltaE = pot_non_int(R_old)+Kinetic_energy(wold,R_old,alpha,N,Dim)
        outfile.write('%f \n' %(DeltaE))
    return

'''Declaration of variable'''
##About the system
N=5
##Here we need to set up Dim=3
Dim=3
alpha=0.22
beta=2.82843
gamma=2.82843
a=0.0043

##Vector used for vectrorize 
BETA=np.ones((N,Dim))
GAMMA=np.ones((N,Dim))
for i in range(N):
    BETA[i][2]=beta
    GAMMA[i][2]=gamma


# EDerivative =0
# eta = 0.01
# delta=0.001
# Niterations = 200


# epsilon=1
# memory=0
# eff=0
# #

# for iter in range(Niterations):
#     alphagradient   = MCS_grad(alpha)
#     print(alpha, alphagradient)
#     alpha -= eta*alphagradient-delta*(alpha-memory)
#     memory=alpha



# print('ɑ='+str(alpha))
# print('efficiency='+str(eff))
# # print(Energy)

##Parallelization
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()
MCC2=2**19##Number of MCC for the final analysis
outfile = open(data_path(f"VMCHarmonic{rank}.dat"),'w')
BETA=np.ones((N,Dim))
GAMMA=np.ones((N,Dim))
for i in range(N):
    BETA[i][2]=beta
    GAMMA[i][2]=gamma



if rank == 0:
    ave=MCC2//nprocs
    data = [ave for p in range(nprocs)]
else:
    data = None

data = comm.scatter(data, root=0)


outfile = open(data_path(f"VMCHarmonic{rank}.dat"),'w')
MCS_last(alpha,data)

##Command to execute

##mpiexec -n 4 python project1g\code.py




