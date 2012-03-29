#!/usr/bin/env python
# Jeff Kaplan 3/2012
#

from numpy import*
from scipy.integrate import odeint

C_CGS = 2.9979245800e10
MSUN_CGS = 1.98892e33
G_CGS = 6.67260e-8
PI = math.pi

#choose physical constants

c    = C_CGS
Msun = MSUN_CGS
G    = G_CGS

# EoS so hard coded right now
Kappa=100.
Gamma=2.
def densityFromPressure(P):
    return pow((P/Kappa),(1./Gamma))

def pressureFromDensity(rho):
    return Kappa* pow(rho,Gamma)

#Initial conditions

#must start from non-zero radius so eqns not singular
r0 = 1e-3 

rhoc=.001
Pc = pressureFromDensity(rhoc)
Mc = 4./3. * PI * r0**3


dx = 0.01

# vars = M(r), P(r)
varsc=[Mc,Pc]



#vars[0] is M, vars[1] is r
def derivTOV(vars,r):
    M = vars[0]
    P = vars[1]
    rho = densityFromPressure(P)
    
    dPdr = -G / r**2 * ( rho + P/c**2) * ( M + 4.*PI * r**3 * P /c**2) 
    dPdr = dPdt /(1. - 2.*G*M/(c**2*r))
    dMdr = 4. * PI * rho *r**2

    return [dMdr,dPdr]

for r in arange(0.0, .003, .0005):
    print r , densityFromPressure(r)







    

