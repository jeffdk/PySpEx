#!/usr/bin/env python
from numpy import *
from matplotlib import pyplot as mpl
import sys


def getHeaderFromComment(line):
    result = line.partition('=')
    return result[2].strip()


class table:
    headers=[]
    data=[]

    def __init__():
        print "creating table"

filesToRead=[]
filesToRead=['Lev0_GhCe.dat']

dataList=[]

timeConversionFactor = 4.926e-6
densityConversionFactor= 6.173e17
pi = arccos(0.)*2.

#Do the files we are reading contain headers?
headerFlag=True
#If not specify them in headers
headers=[]
headers=[]

#dump first line if file contains a pre-header header
dumpFirstLineFlag=True

columnsToPlot=(0,9)

plotList=[]

colorList=['k','r','g','b','brown','k','r','g','b','k','r','g','b']

##normalize by starting values for use in density plots
normalizeByStartingValues=True

#ONLY USE 1 AT A TIME
constraints=True
maxdense=0
restmass=0
maxdensediff=0
convergenceConstraints=True

GrGridRadians=6*pi
GrAngularPoints=array([6,6,6])+array([0,1,2])
GrAngularDh=1.0/GrAngularPoints *GrGridRadians

HydroGridKm=7.5*1.5*1.5;
HydroGridPoints=3*(array([15,15,15])+array([0,1,2])*6)
HydroDh = 1.0/HydroGridPoints * HydroGridKm

print "HydroDh ",  HydroDh
print "GrAngDh ",  GrAngularDh
legendList=['Vlow res','Low res','Mid res','High res', 'Vhigh res']

if constraints:
    filesToRead=['Lev0_GhCe.dat','Lev1_GhCe.dat','Lev2_GhCe.dat']
    headerFlag=True
    columnsToPlot=(0,6)
    legendList=['Low res','Mid res','High res']
    
    if convergenceConstraints:
            columnsToPlot=(0,7)
            legendList=['|Mid - Low|','|High - Mid|']
            columnsToPlot=(0,7)
                        

if maxdense:
    filesToRead=['DensestPoint0.dat','DensestPoint1.dat','DensestPoint2.dat','DensestPoint3.dat']
    filesToRead=['DensestPoint0.dat','DensestPoint1.dat','DensestPoint2.dat']
    filesToRead=['Lev0_DensestPoint.dat','Lev1_DensestPoint.dat','Lev2_DensestPoint.dat']#,'DensestPoint3.dat']
    columnsToPlot=(0,4)
    legendList=['Vlow res','Low res','Mid res','High res', 'Vhigh res']
    legendList=['Low res','Mid res','High res','Vhigh res']
    legendList=['Low res','Mid res','High res']


if restmass:
    filesToRead=['Lev0_RestMass.dat','Lev1_RestMass.dat','Lev2_RestMass.dat']#,'RestMass1.dat','RestMass2.dat']
    columnsToPlot=(0,2)
    legendList=['|Mid - Low|','|High - Mid|']
 
if maxdensediff:
    filesToRead=['DensestPoint0.dat','DensestPoint1.dat','DensestPoint2.dat','DensestPoint3.dat','DensestPoint4.dat']
    filesToRead=['DensestPoint0.dat','DensestPoint1.dat','DensestPoint2.dat']
    filesToRead=['Lev0_DensestPoint.dat','Lev1_DensestPoint.dat','Lev2_DensestPoint.dat']
    columnsToPlot=(0,9)
    legendList=['|Vlow - Low|','|Mid - Low|','|High - Mid|','|Vhigh - High|']
    legendList=['|Mid - Low|','|High - Mid|']
#    legendList=['|Mid - Low|','|High - Mid|','|Vhigh-High|']
#    filesToRead=['DensestPoint0.dat','DensestPoint1.dat','DensestPoint2.dat','DensestPoint3.dat']

fileCount=0
for i in filesToRead:
    currentFile = open(i,'r')
    dataList.append([])
    if dumpFirstLineFlag:
        currentFile.readline()

    lineCount =0
    for line in currentFile:
        if line[0]=='#':
            header_i=getHeaderFromComment(line)
            print header_i
            headers.append(header_i)
            ##if a header we don't process line as normal data
            # so continue
            
            continue
        
        
        splittedLine=line.split()

        if normalizeByStartingValues and (maxdense or maxdensediff):
            densityConversionFactor = 1.0/float(splittedLine[4])
            normalizeByStartingValues=0

#        if len(splittedLine) != len(headers):
#            print "WARNING ERROR LENGTH OF HEADERS NOT SAME AS DATA COLUMNS!"
#            sys.exit()

        for k in range(len(splittedLine)):
            dataList[fileCount].append([])
            #ASSUMING 0 is time column
            if k ==0:
                dataList[fileCount][k].append(float(splittedLine[k])*timeConversionFactor)
            ## if we are working with density, convert to cgs
            elif k == 4 and maxdense:
                dataList[fileCount][k].append(float(splittedLine[k])*densityConversionFactor -1.0)
            else:
               # print fileCount, k
                dataList[fileCount][k].append(float(splittedLine[k])) 
            
            
                
    
        ##This normalizes column 2 by dividing by column 6 
        # and storing in new column 7
        if constraints:
            dataList[fileCount].append([])
            dataList[fileCount][6].append(float(splittedLine[1])/float(splittedLine[5]))

            if convergenceConstraints and fileCount != 0: 
                dataList[fileCount].append([])
                previousFileConstraint =  dataList[fileCount-1][1][len( dataList[fileCount][1])-1]
                dataList[fileCount][7].append(log(float(splittedLine[1])/previousFileConstraint)/
                                           #   log(GrAngularDh[fileCount]/GrAngularDh[fileCount-1]))
                                              log(HydroDh[fileCount]/HydroDh[fileCount-1]) )

        if restmass and fileCount != 0:
            dataList[fileCount].append([])
            previousFileMass =  dataList[fileCount-1][1][len( dataList[fileCount][1])-1]
            dataList[fileCount][2].append(abs(float(splittedLine[1]) - previousFileMass ) )
       

        if maxdensediff and fileCount != 0:
            dataList[fileCount].append([])
            previousFileDensity =  float(dataList[fileCount-1][4][len( dataList[fileCount][4])-1])
            print previousFileDensity, float(splittedLine[4])
            dataList[fileCount][9].append(abs(float(splittedLine[4]) - previousFileDensity )*densityConversionFactor )


        lineCount+=1


    currentFile.close()
    times = dataList[fileCount][columnsToPlot[0]]
#    if restmass:
#        maxpoints =lineCount
#        times=dataList[fileCount][columnsToPlot[0]]
    print len(times), len(dataList[fileCount][columnsToPlot[1]]), fileCount
    #if restmass, we dont want to plot the first guy!
    if (restmass or maxdensediff or convergenceConstraints) and fileCount ==0:
        fileCount=fileCount+1
        continue
        
    
    plotList.append(times)
    plotList.append(dataList[fileCount][columnsToPlot[1]])
    plotList.append(colorList[fileCount])
    fileCount=fileCount+1
    
#print headers
#print dataList

font ={'size':16}
mpl.rc('font',**font)

if constraints:
    mpl.grid(True)
    mpl.title("Generalized harmonic constraint violation for an isolated TOV star\n")
    
    mpl.xlabel(' \n Time elapsed (s)')
    mpl.ylabel('Normalized constraint violation')

    mpl.semilogy(*plotList)
    #mpl.axis([0,.020,8e-4,1.3e-1])
    mpl.legend(legendList,loc=0)
    mpl.show()

if convergenceConstraints:
    mpl.grid(True)
    mpl.title("Convergence rate  for constraints \n")
    
    mpl.xlabel(' \n Time elapsed (s)')
    mpl.ylabel('Normalized constraint violation')

    mpl.plot(*plotList)
    #mpl.axis([0,.020,8e-4,1.3e-1])
    mpl.legend(legendList,loc=0)
    mpl.show()


if maxdense:
    mpl.title("MaxDense: WENO3 fixed FLUID resolution\n")
    mpl.grid(True)
    a=mpl.axes()
    #mpl.setp(a,xticks=[0.0,0.001,0.002,0.003,0.004,0.005],yticks=[0.9985,0.999,0.9995,1.0,1.0005,1.001,1.0015])
    mpl.xlabel('\n Time elapsed (s)')
    mpl.ylabel('Max density / Initial Value -1.0')
    mpl.plot(*plotList)
    mpl.legend(legendList,loc=0)
    mpl.show()

if maxdensediff:
    mpl.title("Convergence of the maximum density for an isolated TOV star\n")
    mpl.grid(True)
    a=mpl.axes()
    #mpl.setp(a,xticks=[0.0,0.004,0.008,0.012,0.016,0.020])
    mpl.xlabel('\n Time elapsed (s)')
    mpl.ylabel('Difference in maximum density/Initial density - 1.0')
    mpl.plot(*plotList)
    mpl.legend(legendList,loc=0)
    mpl.show()

if restmass:
   
    mpl.title("Convergence of rest mass for an isolated 1.625$M_{sun}$ TOV star\n")
    mpl.grid(True)
    #print plotList
    mpl.semilogy(*plotList)
    a=mpl.axes()
    mpl.setp(a,xticks=[0.002,0.004,0.006,0.008])
    mpl.xlabel('\n Time elapsed (s)')
    mpl.ylabel('Difference in Rest Mass $M_{sun}$')
    mpl.legend(legendList,loc=0)
    mpl.show()
