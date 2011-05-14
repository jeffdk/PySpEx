#!/usr/bin/env python
#
#  Goal:  Read in specific ascii tab/space delimited file
#         and put it in an array of tabularFunctions to have
#         a function vs time

from dataFunction import *


#array of dataFunctions managed properly
class dataFunctionVsTime:

    dataFuncs=[]
    times=[]
    
    numTimes=0
    
    #inputs, outputs & points for each dataFunction
    numInputs=0
    numOutputs=0
    numPoints=0


    def __init__(self,numInputs=1,numOutputs=1,numPoints=2,numTimes=1):
    
        self.numInputs  =numInputs
        self.numOutputs =numOutputs
        self.numPoints  =numPoints
        self.numTimes   =numTimes
        print "Creating dataFunction..."
        print "Dim of inputs to dataFuncVsTime: ", self.numInputs
        print "Dim of outputs from dataFuncVsTime: ", self.numOutputs

    def set_numPoints(self, numPoints):
        self.numPoints=numPoints

    def set_numTimes(self, numTimes):
        self.numTimes=numTimes

    #
    def readYgraph(self, filename):
        
        filehandle=open(filename,'r')
        
        #go for the whole file!
        while True:
            line = filehandle.readline()
             #if first line is whitespace, continue until it is not
            while line.isspace(): 
                line = filehandle.readline()
            #if we are at a timestamp line, get the timestamp and read the
            # data at that time into a function vs time
            if len(line) != 0:
                if line[0]== '"' or line[1]=='"' or line[2]=='"':
                    self.times.append( float(line.split('=')[1]))
                    
                    data = dataFunction(self.numInputs,self.numOutputs,self.numPoints)
                    data.readFuncDataFromFile(filehandle)
                    self.dataFuncs.append(data)
                    

            if len(line) == 0:
                break
                
        print "DONE!"



#fvt=dataFunctionVsTime()


#fvt.readYgraph("rho.xg")

    
        
        
#print len(fvt.dataFuncs)
#print "Times in fvt:", len(fvt.times)
#print fvt.times
#for i in range(len(fvt.dataFuncs)):
    #print fvt.dataFuncs[i].data
