
#FDM1
#import numpy
#from matplotlib import pyplot
#numpy.set_printoptions(precision=3)
#L=1
#M=100
#dr=float(L)/float(J-1)
#x_grid=numpy.array([M*dr for M in range(M)])

import scipy.special
import numpy as np
import numpy
from scipy import linalg

#PARAMETERS
h=20
k=12
th=0.01
r1=0.10
r2=0.47
r2c=(r2+(th/2))
L=(r2-r1)
Lc=(L+(th/2))
T1=100
Tinf=20

#EQUATIONS
m=(sqrt((2*h)/(12*th)))
Ab=(pi*(r1)^2)
Af=(2*pi*(((r2c)^2)-(r1)^2))
Ap=(Lc*th)
C2=(((2*r1)/(m))/(((r2c)^2)-(r1)^2))

#BESSEL FUNCTIONS
k0=scipy.special.k0(m*r1)
k1=scipy.special.k1(m*r1)
k1c=scipy.special.k1(m*r2c)
I0=scipy.special.i0(m*r1)
I1=scipy.special.i1(m*r1)
I1c=scipy.special.i1(m*r2c)

#FIN EFFICIENCY AND EFFECTIVENESS
eta=((C2*(((k1*I1c)-(I1*k1c))/((I0*k1c)+(k0*I1c))))) #EFFICIENCY
show("Fin efficiency is: {0:%}".format(eta,digits=8))

#show(r2c/r1) #LINE TO FOLLOW ON EFFICIENCY V.S. PSI PLOT

#xi=(((Lc)^(1.5))*(sqrt((h)/(k*Ap)))) #PSI
#show(xi)

EFFCT=((Af/Ab)*FE) #EFFECTIVENESS
epsilon=n(EFFCT,digits=8)
show("Fin effectiveness is: {0:1}".format(epsilon))

EFFCTO=((h*(Ab+(FE*Af)*(T1-Tinf)/(h*Ab*(T1-Tinf))))) #OVERALL EFFECTIVENESS
#show(EFFCTO)





#FDM2
N=input ('Enter the number of unknowns (N) for problem:');
dr=L/N
dr2=dr*dr
rn=dr , dr, L

np.transpose(rn)

A=numpy.zeros(shape=(N,N))
b=numpy.zeros(shape=(N,0))

#First BC
#A(1,1)=1; b(1)=T1

#CORE OF FIN
x=[]
for i in xrange(1,N-1):
    for j in xrange(1, N-1):
        [A(i,j-1) = 1-((dr)/(2*rn)),
        A(i,j)=-2,
        A(i,j+1)=1+((dr)/(2*rn)),
        b(i)=dr2*m^2(1-Tinf)
        x.append(i,j)



# Last BC
#A(N,N-1)=1; A(N,N-1)=1 b(N-1)=Tinf

#solution=linalg.solve(A,b)
#print solution

#xn=[0;xn]; Tn=[T1;Tn];







