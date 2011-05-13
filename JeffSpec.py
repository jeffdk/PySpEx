#!/usr/bin/env python
# 
# Jeff Kaplan 
#
# 1/2010
#
#  Classes & functions for pseudo-spectral methods
#

##############!!!!!!!!! FUNCTIONS !!!!!!!! #############

from scipy import special as sp
import math
from numpy import *


##Get ChebyshevTs and Legendre Ps and rename
cheb = sp.eval_chebyt
P= sp.legendre
pi = arccos(0.)*2.

ERROR_TOL=1.0e-10


#Heaviside lambda 'triangle' function
def HeavisideLambda(x):
    x=x 
    value = 0.0
    if x < 0.0 and x > - 1.0:
        value = x + 1.0
    elif x == 0.0:
        value = 1.0
    elif(x > 0.0 and x < 1.0):
        value = 1.0 - x
    return value

#Frompyfunc allows a function to operate array wise!
H=frompyfunc(HeavisideLambda,1,1)

#HeavisideFunction offset by off
def Hoffset(x,off):
    return H(x-off)

#2.5 is to make it fit in the interval better
VectorSin=frompyfunc(sin,1,1)
def SinOffset (x, off):
    return VectorSin(2.5*(x-off))


#Extended Heaviside lambda 'triangle' function
# Continuous instead of zero outside of -1,1
def ExtendedHeavisideLambda(x):
    x=x 
    value = 0.0
    if x < 0.0:
        value = x + 1.0
    elif x == 0.0:
        value = 1.0
    elif(x > 0.0):
        value = 1.0 - x
    return value

#Frompyfunc allows a function to operate array wise!
EH=frompyfunc(ExtendedHeavisideLambda,1,1)

def jump(x):
    value=0.0
    if x < 0.0:
        value =0.0
    elif x > 0.0:
        value = 1.0
    elif x == 1.0:
        value = 1.0
    return value
jmp=frompyfunc(jump,1,1)
def jmpoffset(x,off):
    return jmp(x-off)

#HeavisideFunction offset by off
def EHoffset(x,off):
    return EH(x-off)

##A smooth function on -1,1
def smoothfuncA(x):
    return tanh(x**2)*sin((pi*2.*x/3.)**( 13/ 9 ) )#sin(pi*x*x*2.)
smoothfunc=frompyfunc(smoothfuncA,1,1)

def SmoothFuncOffset(x,xoff):
    return smoothfunc(x-xoff)



##A Constant function on -1,1
def constfuncA(x):
    return .5
constfunc=frompyfunc(constfuncA,1,1)

##Arraywise abs
aabs=frompyfunc(abs,1,1)


def forwardEulerStep(rhs,u,deltaT):
    



    dts = rhs(u)
    u_new = u +  deltaT * dts

    return u_new

def rungeKutta4Step(rhs,u,deltaT):
    
    

    k1 = deltaT * rhs(u)
    k2 = deltaT * rhs( u + k1/2.0 )
    k3 = deltaT * rhs( u + k2/2.0 )
    k4 = deltaT * rhs( u + k3 )

    u_new = u + k1/6.0 + k2/3.0 + k3 /3.0 + k4/6.0

    return u_new

def filterTopk(psi, k):
    N=len(psi)
    if  k > N:
        print "CANT FILTER MORE THAN SIZE!!"
        exit()

    matr= dft_matrix(N)

    psi_powers = matrixTransform(matr,psi)
    for i in range(N-k-1,N):
        psi_powers[i] = psi_powers[i]*0.0

    inver = linalg.inv(matr)

    psi = matrixTransform(inver,psi_powers)

    return psi

def waveEqRhs1D(ui):


    dts =  waveEqRhs1DNoBrdy(ui)

    if (len(ui[0]) == len(ui[1])) and  (len(ui[2])==  len(ui[1])):
        N = len(ui[0])
    else:
        print "THHPPP MISTMATCH OF DIMS IN waveEqRhs1D"
        exit()
    
    cheb_xis =[]
    for i in range(N):
        cheb_xis.append(-cos(pi*i/N))
    cheb_xis = array(cheb_xis)

    dtu_plus  = dts[1]+dts[2]
    dtu_minus = dts[1]-dts[2]
    
    from matplotlib import pyplot as mpl
    #mpl.plot(cheb_xis, dtu_plus, cheb_xis, dtu_minus)
    #mpl.legend(["u+", "u-"])
    #mpl.show()
    dtu_plus[0] = 0.0 
    dtu_minus[N-1] = 0.0
    #mpl.plot(cheb_xis, dtu_plus, cheb_xis, dtu_minus)
    #mpl.legend(["u+", "u-"])
    #mpl.show()
    dts[1] = (dtu_plus + dtu_minus)/2.0
    dts[2] = (dtu_plus - dtu_minus)/2.0

    return array(dts)

def waveEqRhs1DNoBrdy(ui):
    
    Psi = ui[0]
    Pi  = ui[1]
    Phi = ui[2]

    N = -1
    dtPsi = []
    dtPi  = []
    dtPhi = []
    if (len(Psi) == len(Pi)) and  (len(Pi)==  len(Phi)):
        N = len(Psi)
    else:
        print "THHPPP MISTMATCH OF DIMS IN waveEqRhs1DNoBrdy"
        exit()

   # print N, len(Psi), len(Pi), len(Phi)

    dtPsi =  Pi*-1.0
    dtPi  =  matrixTransform(deriv_matrix(N),Phi)
    dtPhi =  matrixTransform(deriv_matrix(N),Pi)
    
   # print len(
   # print N, len(dtPsi), len(dtPi), len(dtPhi), "THIS LINE"

 
    return array([dtPsi, dtPi, dtPhi])

# returns N+1xN+1 dft matrix, transforms from collocation to 
# spectral values
def dft_matrix(N):
    N=N-1 
    result = []
    for k  in range(N+1):
        #weird coeffs see eq 3.9
        ck = 1.0
        if k == 0 or k == N:
            ck = 2.0
        
        result.append([])
        #strictly fk is  fk* N ck/2 in eq 3.10
        for n in range(N+1):
            cn= 1.0
            if n == 0 or n == N:
                cn = 2.0
            #nth colocation point given by
            xn=-cos(pi*n/N)
            Tk=cheb(k,xn)
            result[k].append ( 1.0/cn *Tk * 2.0/(N*ck))
            #if k == 35:
            #print "k: ", k,  "\tn: ", n,  "\txn: ", xn,          "\t\t\tTk(xn): ", Tk(xn), "\tfk: ", fk,"\tfunc(xn): ", func(xn)
        #result[k] = result[k] 
        
    return array (result)




def matrixTransform(dftMatrix,points):
    result=[]
    #dimensions good continue
#    print "matrix xform"
#    print len(dftMatrix) 
#    print len(dftMatrix[0])
#    print len(points)
#    print "-----"
    if (len(dftMatrix) == len(points)) and  (len(dftMatrix)==  len(dftMatrix[0])):
        N=len(points)
        result = arange(N)
        result = result * 0.0
        
        #print 'yay'
        sum = 0.0
        for k in range(N):
            sum =0.0
            for n in range(len(points)):
                sum = sum + dftMatrix[k][n]*points[n]

            result[k] = (sum)


    else:
        print "THHPPP MISTMATCH OF DIMS IN matrixTransform"
        exit()
    return result

def deriv_matrix(N):
    N=N-1 #hack for off index
    result = []
    for i in range(N+1):
        result.append([])
        xi = -cos(pi*i/N)
        ci = 1
        if ( i == 0) or ( i == N ):
            ci =2
        for j in range(N+1):
            cj = 1
            if ( j == 0) or ( j == N ):
                cj =2
            xj = -cos(pi*j/N)
            guy = 0.0
            if (i == j) and (j==0):
                #print "case i = j = 0:: i is ", i, " j is ", j
                guy = 1./6. * (1. + 2. * ( N)**2)            
            elif (i == j) and (j == N ):
                guy = -1./6. * (1. + 2. * ( N)**2)
                #print "case i = j = N:: i is ", i, " j is ", j
            elif i == j:
                #print "case i = j !=0,N:: i is ", i, " j is ", j
                guy = 0.5*xj /(1.-xj**2)
            else:
                #print "case i !=      :: i is ", i, " j is ", j
                guy = (-1)**(1+i+j)*ci / (cj * ( xi-xj))
            if (i>N):
                guy = 0
            result[i].append(guy)

    return array(result)
########!!!!!!!!! METHODS !!!!!!!!!!!#########

#AS INSTRUCTED BY PAPER Kidder et al. 2000 sec III
#fks are spectral coeffcients
def calc_fks(func, order):

    N = order

    result = []
    for k  in range(N+1):
        #weird coeffs see eq 3.9
        ck = 1.0
        if k == 0 or k == N:
            ck = 2.0
        
        #fk is current coeff to be appended to result
        fk = 0.0
        #strictly fk is  fk* N ck/2 in eq 3.10
        for n in range(N+1):
            cn= 1.0
            if n == 0 or n == N:
                cn = 2.0
            #nth colocation point given by
            xn=-cos(pi*n/N)
            Tk=cheb(k,xn)
            fk = fk + 1.0/cn *func(xn)*Tk
            #if k == 35:
            #print "k: ", k,  "\tn: ", n,  "\txn: ", xn,          "\t\t\tTk(xn): ", Tk(xn), "\tfk: ", fk,"\tfunc(xn): ", func(xn)
        fk = fk * 2.0/(N*ck)
        result.append(fk)
    return result



#note n should match N
def NthPartialSum(x,N,coefs):
    val=0.0
    for i in range(N+1):
        Ti=cheb(i,x)
        val = val + coefs[i]*Ti
    return val

## Does a log2(N) search for the index value in an monotonic array near r
def index_from_value(r,array ):
    min =0
    max = len(array)-1
    start =array[0]
    end = array[max]
    
    if r < start or r > end:
        print "ERROR in index_from_value, value: '",r, "' is out of range: [", start, ",",end, "]"

    steps = 0
    while min !=max and (min +1) !=max:
        steps = steps+1            
        middle = int( (max-min)/2.) + min
        if  r > array[middle]:
            min = middle
        elif r < array[middle]:
            max = middle
        elif r == array[middle]:
            return_val = middle
            return return_val
        else:
            print >> sys.stderr,"Impossible error finding index from r!"
            sys.exit()
        #print "Returning r index: ", max
    return max


########## NORM FUNCTIONS


from scipy import  integrate as integrate
def L2Norm(func, coeffs):
    N = len(coeffs) -1 
    approx= lambda x: NthPartialSum(x,N,coeffs)
    

    value= integrate.quad( lambda x: (approx(x)-func(x))**2,
                           -1.0,
                           1.0,
                           epsabs=ERROR_TOL,
                           epsrel=ERROR_TOL,
                           limit=50)

    
    return sqrt(value)

from scipy import optimize as opt
def LInfinityNorm(func, coeffs):
    N = len(coeffs) -1 
    approx= lambda x: NthPartialSum(x,N,coeffs)
    funcToMinimize = lambda x: -aabs(approx(x)-func(x))
    
    ftm=frompyfunc(funcToMinimize,1,1)

    ##Random irrational guess for the max 
    #xmin = opt.fminbound(ftm,-1.0,1.0,xtol=ERROR_TOL)

    #xmin = opt.brute(func,[(-1.0+ERROR_TOL,1.0-ERROR_TOL)])

    xmin,fmin = bruteForceBracketMinimize(funcToMinimize,-1.0+ERROR_TOL,1.0-ERROR_TOL,14,ERROR_TOL)
     
    #print xmin

    return -fmin,ERROR_TOL


#Samples function at many different points, then uses minimum point
# as starting location for a downhill minimization function
#returns [xmin, fmin]
def bruteForceBracketMinimize(func,x1,x2,initialSamples,tol):
    fmin = 10e20
    xmin = 0.0
    delta=(x2-x1)
    for i in range(initialSamples + 1):
        xi = x1 + i* delta/initialSamples
        fi = func(xi)
        if fi < fmin:
            fmin= fi
            xmin = xi
    
    returned = opt.fmin(func,xmin,ftol=tol,full_output=1)

    xmin = ((returned[0]))[0]
    fmin = returned[1]

    if returned[4]:
        print "WARNING flag IN bruteForceBracketMinimize: #", returned[4]
        print "1 : Maximum number of function evaluations made. 2 : Maximum number of iterations reached."

    return [xmin,fmin]
    
