#!/usr/bin/env python
#
#  Jeff Kaplan 5/2011
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
from IntegrateTabulatedData import *


eos=dataFunction()
#press=dataFunction()

eosFile=open("LSepsPcs.dat",'r')

print eosFile.readline()
print eosFile.readline()
print eosFile.readline()
print eosFile.readline()
print eosFile.readline()

eos.readFuncDataFromFile(eosFile,[0],[1])
#press.readFuncDataFromFile(eosFile,[0],[2])



print eos.data


expandInLogXspace=1
funcExpanding = log(eos.data + 1.0e-2)


exmin=eos.points[0]
exmax=eos.points[len(eos.points) -1]


eps = eos.points[1] - eos.points[0]

print "eps: ", eps
xs =xtointerval(eos.points,exmin-eps,exmax+eps)
realxs=xtointerval(eos.points,exmin,exmax)

if expandInLogXspace:
    eps = log(eos.points[1]) - log(eos.points[0])
    xs =xtointerval(log(eos.points),log(exmin)-eps,log(exmax)+eps)
    realxs=xtointerval(log(eos.points),log(exmin),log(exmax))



fdeos = []

for i in range(len(realxs)-1):
    fdeos.append((funcExpanding[i+1]-funcExpanding[i])/(realxs[i+1]-realxs[i]))

fdeos.append(fdeos[len(eos.points)-2])
 
print trapezoidIntegrateFromPoints(xs,funcExpanding)

#### N is number of terms in expansion
N = 7

fks=[]
for i in range(N):
    ci = 1.0
    if  (i == N-1 ):
        ci = 2.0
    
    chebi = []
    for j in range(len(xs)):
        chebi.append(
            funcExpanding  *cheb(i,xs[j])
            )
    
    

    fks.append(
        trapezoidIntegrateFromPoints(xs,
                                     funcExpanding / (1.0-xs*xs)**(1./2.)
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

# for i in range(N):
#     mpl.plot(xs,cheb(i,xs))

# mpl.show()

#
#############################


def interpFunc(x):
    func =funcExpanding
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
print "COEFFICIENTS:", fks

pseudopoints = matrixTransform(linalg.inv(dft_matrix(N)),pseudofks)
cpoints = -cos(arange(0,N)*pi/(N-1))
dpseudopoints =  matrixTransform(deriv_matrix(N),pseudopoints)
dpseudofks = matrixTransform(dft_matrix(N), dpseudopoints)
d2pseudopoints =  matrixTransform((deriv_matrix(N)), dpseudopoints)
d2pseudofks = matrixTransform(dft_matrix(N),d2pseudopoints)


## appears det deriv matrix is nearly 0, so inversion is bad
## so we can't integrate by inverting the deriv matrix  =(

print "dpseudofks"
print dpseudofks

print 
print "d2pseudopoints:"
print d2pseudopoints

ys=[]  
ypseudo=[]
dypseudo=[]
d2pseudo = []
fdpseudoys=[0.0]
count = 0
for x in xs:
    ys.append(NthPartialSum(x,N, fks))
    ypseudo.append(NthPartialSum(x,N, pseudofks))
    dypseudo.append(NthPartialSum(x,N, dpseudofks))
    d2pseudo.append(NthPartialSum(x,N, d2pseudofks))
    if count >= 1.0:
        fdpseudoys.append(( ys[count] - ys[count-1] )/ (xs[count] - xs[count-1] ))
    count +=1;
        


mpl.plot(xs, funcExpanding ,xs,ypseudo)#, xs,ys,cpoints,pseudopoints)
mpl.legend(["data","psedo-interp"])#,"galerk-interp","pseudo-points"])
mpl.grid(True)
mpl.show()



# mpl.plot(xs, d2pseudo, cpoints, d2pseudopoints )#, xs,ys,cpoints,pseudopoints)
# mpl.legend(["data","d2psedo-interp"])#,"galerk-interp","pseudo-points"])
# #mpl.show()
# #mpl.semilogx(eos.points,eos.data)
# mpl.grid(True)
# mpl.show()

mpl.plot(cpoints,dpseudopoints,xs,dypseudo,realxs,fdeos)
mpl.legend(["dpseudo-pts","dpseudo-interp","finite-difference"])
mpl.grid(True)
mpl.show()


mpl.semilogy(range(N),absolute(pseudofks),range(N), absolute(fks))
mpl.show()
