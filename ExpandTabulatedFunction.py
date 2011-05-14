#!/usr/bin/env python

#

from matplotlib import pyplot as mpl
from numpy import*
from dataFunction import *
from dataFuncVsTime import *
from scipy import special as sp
cheb = sp.eval_chebyt
P= sp.legendre
pi = arccos(0.)*2.
from JeffSpec import *



def trapezoidIntegrateFromPoints(xs, ys):

    N = len (xs)
    
    sum = 0
    for i in range(N-1):
        sum = sum +  trapezoidal_piece(ys[i],ys[i+1],xs[i],xs[i+1])
                                       

    return sum


def xtointerval(x,xmin,xmax):
    return 2.0*(x - (xmin+xmax)/2.0 )/(xmax-xmin)



def trapezoidal_piece(fa,fb, a, b):
    return (b-a) * ( fa/2 + fb/2  ) 
 

def simpson(f, a ,b, n):
    """f=name of function, a=initial value, b=end value, n=number of double intervals of size 2h"""
 
    n *= 2
    h = (b - a) / n;
    S = f(a)
 
    for i in range(1, n, 2):
        x = a + h * i
        S += 4 * f(x)
 
    for i in range(2, n-1, 2):
        x = a + h * i
        S += 2 * f(x)
 
    S += f(b)
    F = h * S / 3
 
    return F



eos=dataFunction()

eosFile=open("LSepsPcs.dat",'r')

print eosFile.readline()
print eosFile.readline()
print eosFile.readline()
print eosFile.readline()
print eosFile.readline()

eos.readFuncDataFromFile(eosFile,[0],[1])

print eos.data

exmin=eos.points[0]
exmax=eos.points[len(eos.points) -1]

eps = log(eos.points[1]) - log(eos.points[0])

print "eps: ", eps
xs =xtointerval(log(eos.points),log(exmin)-eps,log(exmax)+eps)


 
print trapezoidIntegrateFromPoints(xs,eos.data + 1.e-2)

N = 30
fks=[]
for i in range(N):
    chebi = []
    for j in range(len(xs)):
        chebi.append(
            log(eos.data[j] + 1.e-2)  *cheb(i,xs[j])
            )
    
    

    fks.append(
        trapezoidIntegrateFromPoints(xs,
                                     log(eos.data + 1.e-2) 
                                     * cheb(i,xs) * (2.0 - i % (N-1))/pi
                                     )
        )
    
  
        
print fks
      
ys=[]  
for x in xs:
    ys.append(NthPartialSum(x,N-1, fks))


def interpFunc(x):
    


mpl.plot(xs, log(eos.data + 1.e-2) ,xs, ys)

#mpl.show()
#mpl.semilogx(eos.points,eos.data)

mpl.show()
