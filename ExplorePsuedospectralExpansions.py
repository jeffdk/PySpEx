#!/usr/bin/env python
# 
# Jeff Kaplan 
#
# Explore the effect of shifting a specified function in the interval
# [-1, 1] on a pseudospectral expansion of the function.


import sys, os.path

from matplotlib import pyplot as mpl

from JeffSpec import *

from OneDTestFunctions import *


colorList=['k','r','g','b']
lineStyleList=['--','-.',':','..']

functionUsed=ExtendedHeavisideLambdaOffset
#functionUsed=C1FunctionOddOffset
#functionUsed = JumpOffset
#functionUsed= SmoothFuncOffset
#functionUsed =  SinOffset
normUsed=LInfinityNorm
normLabel="Inf"
#normUsed=L2Norm
#normLabel="2"
otherPlots=1
plotCoefficientPowerVsOffset=0
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
offsetToPlotAt=0.0
offsetStart=-1.0
offsetFinish=1.0
offsetDelta=.01
gotPlotFks=0
normList=[]
offsets=[offsetStart,offsetFinish]

normVsOffsetList=[]


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
            normVsOffsetList.append(colorList[N+numResolutions-1-maxRes])

   # print "values: ", len(values[i])
##THIS ENDS COEFFS(OFFSET META-FLAG)
#print "fks: ",biglist
#mpl.plot(*ColocPoints)

from pylab import *
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib import rc, rcParams
#rc('text', usetex=True)
rc('axes', labelsize=20) #24
rc('axes', titlesize=24) #24
rc('legend', fontsize=18) #16
rcParams['figure.subplot.left'] = 0.10
rcParams['figure.subplot.right'] = 0.96
rcParams['figure.subplot.bottom'] = 0.10
rcParams['figure.subplot.top'] = 0.95
rcParams['figure.figsize'] = 10, 8
rc('xtick', labelsize=16)
rc('ytick', labelsize=16)
if otherPlots:
##COEFFS FOR OFFSET

    if plotCoefficientPowerVsOffset:
        mpl.figure(figsize=(9,7))
        #mpl.subplot(221)
        mpl.semilogy(*offsetlist)
        mpl.axis([offsetStart,offsetFinish,1e-7,1])
        mpl.legend(legendList,loc=0)

    figure(figsize=(9.3,7.2))
    semilogy(*normVsOffsetList)
    for k in range(numResolutions):
        currentRes = resolutions[k]
        for i in range(currentRes+1):
            point = (-cos(pi*i/currentRes))
        #print "line at ",point,"\tres: ",currentRes,lineStyleList[k], colorList[k]
            axvline(point,linestyle=lineStyleList[k], color=colorList[k])

    axis([offsetStart,offsetFinish,1e-3,1e-0])
    legend(["N = %s" % res for res in resolutions], loc=8)
    ylabel("L%s Norm of truncation error" % normLabel)
    xlabel("Location of non-smooth point")
    title("Global error for N'th order expansion")

    #mpl.subplot(222)
    #mpl.figure(figsize=(8,6))
    #mpl.plot(*ColocPoints)
    #mpl.grid()

    #for i in range(10):
    #    mpl.axvline(i+1)
    #mpl.axvline(ColocPoints[i])
    #mpl.show()

figure(figsize=(16,10))
ax = subplot(121)
subplots_adjust( bottom=0.15)

offset_index = int(len(offsets/2))

print normList
#####Setup power spectri plot ###
powerSpectri=[]
psLegend=[]
for k in range(numResolutions): 

    psLegend.append("Order: %s\n  L%s error: %s\n  +/-: %s\n" %    (resolutions[k], normLabel,  normList[k][0],  normList[k][1]))

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
ax.set_xlabel("n (Coefficient #)")
ax.set_ylabel("Value of n'th coefficient")
ax.set_title("Spectral space")


axis([0, resolutions[numResolutions-1], 1e-6, 1])
legend(psLegend, loc=0)


#####Setup function plot ###
ax = subplot(122)
xs= arange(-1.0,1.0,.01)

funcPlots =[]
funcPlots.append(xs)
funcPlots.append(functionUsed(xs,offsetToPlotAt))
funcPlots.append('purple')
funcLegend=["f(x)"]
pltList = []
for k in range(numResolutions):
    funcPlots.append(xs)
    funcPlots.append(NthPartialSum(xs,resolutions[k],examplePlot[k]))
    funcLegend.append("N = %s" %resolutions[k])
    color = colorList[k]
    funcPlots.append(color)

    ##plot colocation points for res's
    #k=numResolutions-1
    currentRes = resolutions[k]
    for i in range(currentRes+1):
        point = (-cos(pi*i/currentRes))
        #print "line at ",point,"\tres: ",currentRes,lineStyleList[k], colorList[k]
        axvline(point,linestyle=lineStyleList[k], color=colorList[k])
fig2 = plot(*funcPlots)
ax.set_xlabel("x")
ax.set_ylabel("f(x) = %s" % functionUsed.__name__)
ax.set_title("Physical space")
legend(fig2, funcLegend, loc=8)
axis([-1,1, -1.2, 1.2])

axcolor = 'lightgoldenrodyellow'
axfreq = axes([0.15, 0.05, 0.55, 0.03], axisbg=axcolor)


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
            psLegend.append("Order: %s\n  L%s error: %s\n  +/-: %s\n" % 
                            (resolutions[l],normLabel,  norm[0], norm[1]))
    
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
