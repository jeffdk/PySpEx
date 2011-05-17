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

    length = len (xs)
    
    sum = 0
    for i in range(length-1):
        sum = sum +  trapezoidal_piece(ys[i],ys[i+1],xs[i],xs[i+1])
                                       

    return sum


def xtointerval(x,xmin,xmax):
    return 2.0*(x - (xmin+xmax)/2.0 )/(xmax-xmin)



def trapezoidal_piece(fa,fb, a, b):
    return (b-a) * ( fa/2.0 + fb/2.0  ) 
 

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
realxs=xtointerval(log(eos.points),log(exmin),log(exmax))

 
print trapezoidIntegrateFromPoints(xs,eos.data + 1.e-2)

N = 7
fks=[]
for i in range(N):
    ci = 1.0
    if  (i == N-1 ):
        ci = 2.0
    
    chebi = []
    for j in range(len(xs)):
        chebi.append(
            log(eos.data[j] + 1.e-2)  *cheb(i,xs[j])
            )
    
    

    fks.append(
        trapezoidIntegrateFromPoints(xs,
                                     log(eos.data + 1.e-2) / (1.0-xs*xs)**(1./2.)
                                     * cheb(i,xs) * (ci)/pi
                                     )
        )
   

################### TEST INTEGRATOR #######################
#
#  tests integration of function f(x)=x
#



testxs = arange(0.0, 1.0, 1e-2)

testintegrand = []

for i in range(len(testxs)):
    portion = testxs[0:i]
    testintegrand.append (
        trapezoidIntegrateFromPoints(portion, portion)
        )

#mpl.plot(testxs,testintegrand-1./2.*testxs*testxs)
#mpl.show()
# done with test integration
#######################


##########################
#Test chebbies

#for i in range(N):
#    mpl.plot(xs,cheb(i,xs))

#mpl.show()

#
#############################


def interpFunc(x):
    func = log(eos.data + 1.e-2) 
    index =index_from_value(x,realxs )

    a = realxs[index -1]
    b = realxs[index]
    fa = func[index-1]
    fb = func[index]
    value = fa + (fb - fa)/(b-a) * (x -a)

    #print x, xs[index]
    return value

print interpFunc(.5)

pseudofks= calc_fks(interpFunc,N-1)
 

print pseudofks
print fks

pseudopoints = matrixTransform(linalg.inv(dft_matrix(N)),pseudofks)
cpoints = -cos(arange(0,N)*pi/(N-1))
dpseudopoints =  matrixTransform(deriv_matrix(N),pseudopoints)
dpseudofks = matrixTransform(dft_matrix(N), dpseudopoints)

ys=[]  
ypseudo=[]
dypseudo=[]
for x in xs:
    ys.append(NthPartialSum(x,N, fks))
    ypseudo.append(NthPartialSum(x,N, pseudofks))
    dypseudo.append(NthPartialSum(x,N, dpseudofks))

mpl.plot(xs, log(eos.data + 1.e-2) ,xs,ypseudo)#, xs,ys,cpoints,pseudopoints)
mpl.legend(["data","psedo-interp"])#,"galerk-interp","pseudo-points"])
mpl.grid(True)
mpl.show()

mpl.plot(xs, eos.data  ,xs, exp(ypseudo)  - 1.e-2)#, xs,ys,cpoints,pseudopoints)
mpl.legend(["data","psedo-interp"])#,"galerk-interp","pseudo-points"])
#mpl.show()
#mpl.semilogx(eos.points,eos.data)
mpl.grid(True)
mpl.show()

mpl.plot(cpoints,dpseudopoints,xs,dypseudo)
mpl.legend(["dpseudo-pts","dpseudo-interp"])
mpl.grid(True)
mpl.show()


mpl.semilogy(range(N),absolute(pseudofks),range(N), absolute(fks))
mpl.show()
