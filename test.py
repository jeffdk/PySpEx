#!/usr/bin/env python
# 
# Jeff Kaplan 
#
# Referene & notation is that of Kidder et al.
#  Blach hole evolutions by spectral methods
#

import sys, os.path
#from scipy import special as sp

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



resesToLoopThrough=maxRes+1
resolutionsToRun=arange(1,resesToLoopThrough)

for N in resolutionsToRun:
    print "Running res: N",N
    
    xs = arange(-1.,1.1,.1)
    Ns = range(N+1)
    
    ColocPoints.append(resolutionsToRun)
    ColocPoints.append([])
    for i in resolutionsToRun:
        k = i
        if k < N:
            ColocPoints[2*N-1].append(-cos(pi*k/N))
        else:
            ColocPoints[2*N-1].append(1.0)
                   ###UNCOMMENT FOR COEFSS(OFFSET)

    for l in range(numResolutions):
        if N == resolutions[l]:
            biglist.append([])
            examplePlot.append([])
            gotPlotFks = 0
            offsets= arange(offsetStart,offsetFinish,offsetDelta)
            offsetlist=[]
            offsetOfN=[]
            values.append([])

            currentNormList=[]

            for i in Ns:
                values[l].append([])
                #print values
                legendList.append("Coeff#%s"%i)
            for offset in offsets:
                func = lambda x: functionUsed(x,offset)
                
                

                fks = calc_fks(func,N)
                biglist[l].append(Ns)
                biglist[l].append((fks))

                currentNormList.append(normUsed(func, fks)[0])


                if offset >= offsetToPlotAt and not gotPlotFks:
                    #print "Got offset to plot at!"
                    gotPlotFks=1
                    #print fks
                    examplePlot[l] = (fks)
                    normList.append(normUsed(func, fks))

                for k in range(N):
                    values[l][k].append(aabs(fks[k]))


                for k in range(N):        
                    offsetlist.append(offsets)
                            # print "ofsets: ", len(offsets)
                    offsetlist.append(values[l][k])

            normVsOffsetList.append(offsets)
            normVsOffsetList.append(currentNormList)

   # print "values: ", len(values[i])
##THIS ENDS COEFFS(OFFSET META-FLAG)
#print "fks: ",biglist
#mpl.plot(*ColocPoints)



if otherPlots:
##COEFFS FOR OFFSET
    mpl.figure(figsize=(9,7))
    #mpl.subplot(221)
    mpl.semilogy(*offsetlist)
    mpl.axis([offsetStart,offsetFinish,1e-7,1])
    mpl.legend(legendList,loc=0)

    mpl.figure(figsize=(9.3,7.2))
    mpl.semilogy(*normVsOffsetList)
    mpl.axis([offsetStart,offsetFinish,1e-5,10])
    mpl.legend(resolutions)




    #mpl.subplot(222)
    #mpl.figure(figsize=(8,6))
    #mpl.plot(*ColocPoints)
    #mpl.grid()

    #for i in range(10):
    #    mpl.axvline(i+1)
    #mpl.axvline(ColocPoints[i])
    #mpl.show()


from pylab import *
from matplotlib.widgets import Slider, Button, RadioButtons
figure(figsize=(16,10))
ax = subplot(121)
subplots_adjust( bottom=0.15)

offset_index = int(len(offsets/2))

print normList
#####Setup power spectri plot ###
powerSpectri=[]
psLegend=[]
for k in range(numResolutions): 

    psLegend.append("Points: %s\n  L%sDiff: %s\n  +/-: %s\n" %    (resolutions[k],normLabel,  normList[k][0],  normList[k][1]))

    color = colorList[k]
    
    curXs =range(resolutions[k]+1) 
    curYs = aabs(biglist[k][2*offset_index-1])
    
    powerSpectri.append(curXs)
    powerSpectri.append(curYs)
    powerSpectri.append( color)
    #print k,  powerSpectri[2*k]
    #print k, powerSpectri[2*k+1]
    
    
#print powerSpectri
#print

fig1 =semilogy(*powerSpectri) 

axis([0, resolutions[numResolutions-1], 1e-6, 1])
legend(psLegend, loc=0)


#####Setup function plot ###
ax = subplot(122)
xs= arange(-1.0,1.0,.01)

funcPlots =[]
funcPlots.append(xs)
funcPlots.append(functionUsed(xs,offsetToPlotAt))
funcPlots.append('purple')
funcLegend=["func"]
for k in range(numResolutions):
    funcPlots.append(xs)
    funcPlots.append(NthPartialSum(xs,resolutions[k],examplePlot[k]))
    funcLegend.append(resolutions[k])
    color = colorList[k]
    funcPlots.append(color)

fig2 =plot(*funcPlots)
axis([-1,1, -1.2, 1.2])

##plot colocation points for res's
k=numResolutions-1
currentRes = resolutions[k]
for i in range(currentRes+1):
    point = (-cos(pi*i/currentRes))
    #print "line at ",point,"\tres: ",currentRes,lineStyleList[k], colorList[k]
    axvline(point,linestyle=lineStyleList[k], color=colorList[k])


axcolor = 'lightgoldenrodyellow'
axfreq = axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)


sOffset = Slider(axfreq, 'Offset', offsetStart, offsetFinish, valinit=offsetToPlotAt)

##update all plot values given new offset

def update(val):
    Offset = sOffset.val
    psLegend=[]
    offset_index=index_from_value(Offset,offsets)
    func = lambda x: functionUsed(x,Offset)

    fig2[0].set_ydata(functionUsed(xs,Offset))
    for l in range(numResolutions):
        
        freqs=biglist[l][2*offset_index-1]
        
        if calcNormInteractive:
            norm = normUsed(func, freqs)
            psLegend.append("Points: %s\n  L2diff: %s\n  +/-: %s\n" % 
                            (resolutions[l], norm[0], norm[1]))
    
        fig1[l].set_ydata(aabs(freqs))
        
        
        fig2[l+1].set_ydata(NthPartialSum(xs,resolutions[l],freqs))
    if calcNormInteractive:
        ax = subplot(121)
    legend(psLegend,loc=0)
    draw()
sOffset.on_changed(update)


resetax = axes([0.8, 0.025, 0.1, 0.04])
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')
def reset(event):
    sOffset.reset()
   
button.on_clicked(reset)
"""
rax = axes([0.025, 0.5, 0.15, 0.15], axisbg=axcolor)
radio = RadioButtons(rax, ('red', 'blue', 'green'), active=0)
def colorfunc(label):
    l.set_color(label)
    draw()
radio.on_clicked(colorfunc)
"""
show()


"""
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

mpl.rcParams['legend.fontsize'] = 10

fig = plt.figure()
ax = Axes3D(fig)
theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)
z = np.linspace(-2, 2, 100)
r = z**2 + 1
x = r * np.sin(theta)
y = r * np.cos(theta)
ax.semilogy(range(N+1), aabs(biglist), z, label='parametric curve')
ax.legend()

plt.show()
"""


#mpl.plot(xs,Ha(xs),xs,NthPartialSum(xs, N,fks))
#mpl.semilogy(range(N+1),aabs(fks))
#mpl.show()
