#!/usr/bin/env python
from numpy import *
from matplotlib import pyplot as mpl
#mpl.rcParams['text.usetex'] = True
#mpl.rc('text',usetex=True)
#mpl.rc('font', family='serif')
#font ={'size':24}
#mpl.rc('font',**font)
#mpl.rc('axes', labelsize='large')

import sys

filename = "D0res800x800.dat"

file = open(filename,'r')

print "First Line:"
print file.readline()

line=file.readline()

ns,nu=line.split()

print "Number of radial points, ns: ", ns
print "Number of mu, i.e. cos(theta) points, nu: ", nu

line=file.readline()
radequat,npoly=line.split()

line=file.readline()
comega,rotfunc=line.split()

print "Radequat (avg radius?): ", radequat

#variables defined in file:  r,cosp,ed,alpha,rho,gamma,omega,angv
#r,cosp,ed,alpha,rho,gamma,omega,angv=[],[],[],[],[],[],[],[]
#appearently ed = rhob

data=[]

rhoEquator=[]
rsEquator=[]

numLines=0
for line in file:
    data.append(line.split())
    #if mu = 1.0, theta = Pi/2 
    #print data[numLines][1]
    if float(data[numLines][1])==1.0 and float(data[numLines][2]) != 0.0:
        rsEquator.append(float(data[numLines][0]))
        rhoEquator.append(float(data[numLines][2]))
    numLines+=1

print "TotalLines: ", numLines, "\t ns*nu: ", int(ns)*int(nu)
#print data

#print rsEquator

# mpl.title("Rho from CST file")
# mpl.semilogy(rsEquator,rhoEquator)
# mpl.show()

print range (5, 8)

mpl.title("Rho from CST file")
mpl.plot(rsEquator,rhoEquator)
mpl.show()

avgLength=5

avgdRhos=[]

for i in range(len(rsEquator)):
    istart = i - avgLength/2
    if istart < 0:
        istart = 0
    if istart +avgLength > len(rsEquator) -2:
        istart = len(rsEquator) -2 - avgLength
    
    sum = 0.0
    xrange = rsEquator[istart + avgLength+1] -rsEquator[istart]
    weightSum = 0.0
    for j in range(istart, istart + avgLength+1):
        dx = rsEquator[j+1] - rsEquator[j] 
        trapzoid = (rhoEquator[j+1] + rhoEquator[j])/2.0
        weight = 1./(abs(j-i) + 1.)
        weightSum+=weight
        sum += trapzoid*dx
    sum = sum / xrange

    avgdRhos.append(sum)



mpl.title("Rho from CST file")
mpl.plot(rsEquator, array(avgdRhos)- array(rhoEquator))
mpl.show()
