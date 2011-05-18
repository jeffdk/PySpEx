#!/usr/bin/env python
# Jeff Kaplan 5/2011
#

from numpy import*


#Should really eliminate trapezoidal_piece
def trapezoidIntegrateFromPoints(xs, ys):

    length = len (xs)
    
    sum = 0
    for i in range(length-1):
        sum = sum +  trapezoidal_piece(ys[i],ys[i+1],xs[i],xs[i+1])
                                       

    return sum


def trapezoidal_piece(fa,fb, a, b):
    return (b-a) * ( fa/2.0 + fb/2.0  ) 
 

#Copied from wikipedia; have not tested yet
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
