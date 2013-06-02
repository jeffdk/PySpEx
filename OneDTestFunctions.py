#!/usr/bin/env python
# 
# Jeff Kaplan 
#
#
#  These functions are useful in testing spectral methods on
#  the interval [ -1, 1 ]
#

from numpy import *

pi = arccos(0.)*2.


#Heaviside lambda 'triangle' function
def HeavisideLambdaScalar(x):
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
HeavisideLambda=frompyfunc(HeavisideLambdaScalar,1,1)

#HeavisideFunction offset by off
def HeavisideLambdaOffset(x,off):
    return HeavisideLambda(x-off)


VectorSin=frompyfunc(sin,1,1)
def SinOffset (x, off):
    return VectorSin(2.5*(x-off))


#Extended Heaviside lambda 'triangle' function
# Continuous instead of zero outside of -1,1
def ExtendedHeavisideLambdaScalar(x):
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
ExtendedHeavisideLambda=frompyfunc(ExtendedHeavisideLambdaScalar,1,1)

#HeavisideFunction offset by off
def ExtendedHeavisideLambdaOffset(x,off):
    return ExtendedHeavisideLambda(x-off)

def JumpScalar(x):
    value=0.0
    if x < 0.0:
        value =0.0
    elif x > 0.0:
        value = 1.0
    elif x == 1.0:
        value = 1.0
    return value

Jump=frompyfunc(JumpScalar,1,1)
def JumpOffset(x,off):
    return Jump(x-off)


##A smooth function on -1,1
def SmoothFuncScalar(x):
    return tanh(x**2)*sin((pi*2.*x/3.)**( 13/ 9 ) )#sin(pi*x*x*2.)
SmoothFunc=frompyfunc(SmoothFuncScalar,1,1)

def SmoothFuncOffset(x,xoff):
    return SmoothFunc(x-xoff)


## A function which is C^1 smooth (second deriv discontinuous)
def C1FunctionScalar(x):
    x = x + 1.0  # translate to center
    value =  x - 1.0
    if x < 1.0:
        value = 0.5* x ** 2 - 0.5
    return value
C1Function = frompyfunc(C1FunctionScalar,1,1)
def C1FunctionOffset(x, off):
    return C1Function(x - off)

## A function which is C^1 smooth (second deriv discontinuous)
def C1FunctionOddScalar(x):
    value =  -x ** 2
    if x < 0.0:
        value =  x ** 2 
    return value
C1FunctionOdd = frompyfunc(C1FunctionOddScalar,1,1)
def C1FunctionOddOffset(x, off):
    return C1FunctionOdd(x - off)


##A Constant function on -1,1
def ConstFuncScalar(x):
    return 0.5
ConstFunc=frompyfunc(ConstFuncScalar,1,1)

##Arraywise abs
aabs=absolute
