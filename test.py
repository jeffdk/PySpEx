#!/usr/bin/env python
# 
# Jeff Kaplan 
#
# Referene & notation is that of Kidder et al.
#  Blach hole evolutions by spectral methods
#

import sys, os.path

from matplotlib import pyplot as mpl

from JeffSpec import *

from OneDTestFunctions import *

colorList=['k','r','g','b']
lineStyleList=['--','..',':','-.']

functionUsed=ExtendedHeavisideLambdaOffset
#functionUsed = JumpOffset
#functionUsed= SmoothFuncOffset
#functionUsed =  SinOffset
normUsed=LInfinityNorm
normLabel="Infin"
otherPlots=1
calcNormInteractive=1
ColocPoints=[]
biglist=[]
offsetlist=[]
offsets=[]
values=[]

resolutions = [9,10,11]
numResolutions =  len(resolutions)
maxRes = resolutions[numResolutions-1]

legendList=[]
examplePlot=[]
offsetToPlotAt=0.4
offsetStart=-1.0
offsetFinish=1.0
offsetDelta=.01
gotPlotFks=0
normList=[]
offsets=[offsetStart,offsetFinish]

normVsOffsetList=[]

N=30
print dft_matrix(N)
print len(dft_matrix(N))

Psi0=[]
Pi0 =[]
Phi0=[]

points =  []
cheb_xis= []



for i in range(N+1):
    x = -cos(pi*i/N)
    cheb_xis.append(-cos(pi*i/N))
    points.append(functionUsed(-cos(pi*i/N) , 0 ))
    Psi0.append(0.0)
    Pi0.append( 10.*x*(1-x*x)*exp( -26.0 * (cos(pi*i/N)**2 ) ) )
    Phi0.append(0.0)


Psi0 = array (Psi0)
Pi0  = array (Pi0)
Phi0 = array (Phi0)
cheb_xis = array(cheb_xis)

#enforce boundry cond
u_plus  = Pi0 + Phi0
u_minus = Pi0 - Phi0

u_plus[0] = 0.0
u_minus[N] = 0.0

Pi0  = (u_plus + u_minus )/2.0
Phi0 = (u_plus - u_minus )/2.0

#dts =  waveEqRhs1D(Psi0, Pi0, Phi0)


#mpl.plot( cheb_xis,Psi0,    cheb_xis, Pi0,    cheb_xis,  Phi0,
 #         cheb_xis, dts[0], cheb_xis, dts[1],  cheb_xis, dts[2] )

#mpl.show()
#dts =  waveEqRhs1DNoBrdy(Psi0, Pi0, Phi0)
#mpl.plot( cheb_xis,Psi0,    cheb_xis, Pi0,    cheb_xis,  Phi0,
  #        cheb_xis, dts[0], cheb_xis, dts[1],  cheb_xis, dts[2] )
#mpl.show()

deltaT = 1.e-2

fout =  open('psi.ygraph', 'w')

u_0 = array([Psi0, Pi0, Phi0 ])

u = forwardEulerStep(waveEqRhs1D, u_0,deltaT)

outputEvery=1.e-2

for t in arange(0,1.50,deltaT):
    
    if ( t /outputEvery) == int(t/outputEvery):

        st = '\"Time='+str(t ) + '\n'
        fout.write( st )
        for i in range(N+1):
            fout.write(str(cheb_xis[i]) +'\t' + str(u[0][i]) +'\t' + str(u[1][i]) +'\t' + str(u[2][i])  +'\n' )

        fout. write('\n')

    #u_new = forwardEulerStep(waveEqRhs1D,u,deltaT)
    u_new = rungeKutta4Step(waveEqRhs1D,u,deltaT)
    #u_new[0] = filterTopk(u_new[0],3)
    #u_new[1] = filterTopk(u_new[1],3)
    #u_new[2] = filterTopk(u_new[2],3)
    #u_plus  =u_new[1]  + u_new[2]
    #u_minus = u_new[1] - u_new[2]

    #u_plus[0] = 0.0
    #u_minus[N] = 0.0

    #u_new[1] = (u_plus + u_minus )/2.0
    #u_new[2]= (u_plus - u_minus )/2.0
    u = u_new

fout.close()

mpl.plot(  cheb_xis,Psi0,  cheb_xis, u[1], );
mpl.show()
#  cheb_xis, Pi0,    cheb_xis,  Phi0,
#          cheb_xis, dts[0], cheb_xis, dts[1],  cheb_xis, dts[2] )


    


#exit()
#print points
#print matrixTransform(dft_matrix(N),points)
#print calc_fks( lambda x: functionUsed(x,0), N)
#print
#print deriv_matrix(N)
#print
#print  matrixTransform(deriv_matrix(N),points)

#xs = arange(-1.,1.,.01)
#ys=[]
#dys=[]
#fd=[0,0]
#count = 0
#for x in xs:
#    ys.append(NthPartialSum(x,N, matrixTransform(dft_matrix(N),points)))
#    dys.append(NthPartialSum(x,N, matrixTransform(dft_matrix(N),matrixTransform(deriv_matrix(N),points))))
#    if count > 1:
#        fd.append((ys[count] - ys[count-2])/.02)
#    count +=1
#mpl.plot(xs,ys,xs,dys,xs,fd)
#mpl.show()
#exit()

