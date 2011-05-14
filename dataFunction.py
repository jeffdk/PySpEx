#!/usr/bin/env python
#
#  Goal:  Read in any ascii tab/space delimited file
#         put in a format it is easy to do data manipulations on
#
# TO DO LIST: 
#   
#  1)  SEE LINE 27: add check (size points)/numInputs = size
#                  MUST ALSO DO THIS IN CHECK IN setPoints/setData (that is normalize by numInputs/numOutputs)
# 
#  2)  Add manual choice of columns for input/output data
#      This includes having #columns not equal to numInputs+numOutputs
#
from numpy import *

import sys
#creates a function out of tabular data
#
# inputs index an array of size numOutputs
#
class dataFunction:

    numInputs=0
    numOutputs=0
    
    numPoints=0

    ###must add check (size points)/numInputs = size
    points = array([0.0])
    data= array([0.0])

    def __init__(self,numInputs=1,numOutputs=1,numPoints=2):

        self.numInputs  =numInputs
        self.numOutputs =numOutputs
        self.numPoints  =numPoints 
        #print "Creating dataFunction..."
        #print "Dim of inputs to dataFunction: ", self.numInputs
        #print "Dim of outputs from dataFunction: ", self.numOutputs

    def __str__(self):
        return "<dataFunction object instance>"

    def set_numPoints(self, numPoints):
        self.numPoints=numPoints
        
    def setPoints(self, points):
        ##
        #Data points must be numpy array of floats of dim 1 and correct length!
        ##
        if type(points) != type(array([])):
            print "Input to dataFunction.setPoints not a numpy array: ",  type(points)
            print " WARNING: POINTS NOT SET!"
            return -1
        if points.dtype != dtype('float64'):
            print "Input array to dataFunction.setPoints has wrong data type: ", points.dtype
            print " WARNING: POINTS NOT SET!"
            return -1
        if points.ndim != 1:
            print "Input array to setPoints is not 1-dimesional as dataFunction expects: ", points.ndim
            print " WARNING: POINTS NOT SET!"
            return -1
        if self.numPoints*self.numInputs != points.size:
            print "Input array to dataFunction.setPoints has wrong number of points: ", points.size
            print "This dataFunction has it's numPoints set to: ", self.numPoints
            print " WARNING: POINTS NOT SET!"
            return -1
            
        self.points=points
    #done with setPoints
    def setData(self, data):
        ##
        #Data data must be numpy array of floats of dim 1 and correct length!
        ##
        if type(data) != type(array([])):
            print "Input to dataFunction.setData not a numpy array: ",  type(data)
            print " WARNING: DATA NOT SET!"
            return -1
        if data.dtype != dtype('float64'):
            print "Input array to dataFunction.setData has wrong data type: ", data.dtype
            print " WARNING: DATA NOT SET!"
            return -1
        if data.ndim != 1:
            print "Input array to setData is not 1-dimesional as dataFunction expects: ", data.ndim
            print " WARNING: DATA NOT SET!"
            return -1
        if self.numPoints*self.numOutputs != data.size:
            print "Input array to dataFunction.setData has wrong number of dat: ", data.size
            print "This dataFunction has it's numData set to: ", self.numData
            print " WARNING: DATA NOT SET!"
            return -1
            
        self.data=data
    #Done with setData

    # NOTE: Filehandle must ALREADY be past all header data
    #   This function can ONLY HANDLE WHITESPACE AND WHITESPACE
    #    DELIMINTED DATA!!!
    # 
    # dumps any whitespace at start, then reads in column data
    #  and stops reading at next line of whitespace
    def readFuncDataFromFile(self, filehandle, inputsColumnList=0, outputsColumnList=0):
        minCols = self.numInputs+self.numOutputs
        #default arguments means first cols are inputs
        if inputsColumnList==0:
            inputsColumnList=range(self.numInputs)
        # and then next cols are outputs
        if outputsColumnList==0:
            outputsColumnList=range(self.numInputs, minCols)
        
       
        if type(filehandle) != type(open("dumb.txt",'w')):
            print "Must pass a file-handle to readFromFile, not a: ", type(filehandle)
            print " WARNING: DATA NOT READ FROM FILE in dataFunction.readFromFile#"
            return -1
       
        #read first line
        line = filehandle.readline()
       
        #if first line is whitespace, continue until it is not
        while line.isspace(): 
            line = filehandle.readline()
        
        
        #now split the line, it must have more columns (entries)
        # than number of total inputs + outputs
        lineData=line.split()
        if len(lineData) < minCols:
            print "Must have numCols in data file = numInputs+numOutputs: "
            print "\t columns in data file: {0}".format(len(lineData))
            print "\t numInputs+numOutputs: {0}".format(minCols)
            print "FATAL ERROR IN readFuncDataFromFile!"
            sys.exit()
        #print lineData

        #one for inputs one for outputs
        tempData=[[],[]]
        
        
        for j in inputsColumnList:
            tempData[0].append( float(lineData[j]) )
       
        #print inputsColumnList 
        #print outputsColumnList
        for j in outputsColumnList:
            tempData[1].append( float(lineData[j]) )
       
        # already have first line
        numberOfLines=1
        while True:
            line=filehandle.readline()
           
            if line.isspace() or not line:
                break
            lineData=line.split()
            numberOfLines+=1
            
            for i in inputsColumnList:
                #print lineData, i
                tempData[0].append(float(lineData[i])) 
            for i in outputsColumnList:
                tempData[1].append(float(lineData[i])) 
        #print numberOfLines
        
        self.set_numPoints(numberOfLines)

        self.setPoints(array(tempData[0]))
        self.setData(array(tempData[1]))
    
    #currently only supports 1D input/outputs!!!!
    #also assumes gridpoints line up!!!! Will be a catastrophic failure otherwise!
    #(though one may have 2x, 3x, etc. points other)
    def __sub__(self, rightArg):
        
        #if we are adding the same thing
        if str(self)==str(rightArg):
            #answer is only defined on a range where both functoins
            # have values
            pointMax=min(self.points[self.numPoints-1],
                          rightArg.points[rightArg.numPoints-1]
                          )
            numPointsMax=max(self.numPoints,rightArg.numPoints)
            #print pointMax, numPointsMax
            
            #compute which one has least resolution (highest grid spacing)
            #, then use that for the grid-spacing
            selfRange= self.points[self.numPoints-1] - self.points[0]
            rightArgRange= rightArg.points[rightArg.numPoints-1] - rightArg.points[0]
            selfSpacing=selfRange/(self.numPoints - 1)
            rightArgSpacing=rightArgRange/(rightArg.numPoints -1)

           # print selfSpacing, rightArgSpacing
            
            
            if selfSpacing > rightArgSpacing:
                grid=self.points
            else:
                grid=rightArg.points
                

            i=0
            resultData=[]
            resultGrid=[]
            currentPoint = -9e99
            selfIndex=0
            rightArgIndex=0
            while currentPoint < pointMax and i < len(grid):
                currentPoint=grid[i]
                currentSelfPoint=self.points[selfIndex]
                currentSelfValue=self.data[selfIndex]
                currentRightArgPoint=rightArg.points[rightArgIndex]
                currentRightArgValue=rightArg.data[rightArgIndex]

                print currentSelfPoint, currentPoint, currentRightArgPoint
                if round(currentSelfPoint,5)==round(currentPoint,5) and round(currentPoint,5)==round(currentRightArgPoint,5):
                    #print "all equal" 
                    result=currentSelfValue-currentRightArgValue
                    resultGrid.append(currentPoint)
                    resultData.append(result)
                    i+=1
                    selfIndex+=1
                    rightArgIndex+=1
                
                elif currentSelfPoint < currentPoint:
                    #print "self less than: ", currentSelfPoint,selfIndex,i
                    selfIndex+=1
                elif currentRightArgPoint < currentPoint:
                    #print "arg less than: ", currentRightArgPoint, rightArgIndex,i
                    rightArgIndex+=1

            resultFunc=dataFunction()
            resultFunc.set_numPoints(i)
            resultFunc.setData(array(resultData))
            resultFunc.setPoints(array(resultGrid))
        #or, we define subtraction by a int or float as subtract
        # that value from every point in the function
        elif type(rightArg)==type(1) or type (rightArg)==type(1.1):
            resultFunc=dataFunction()
            resultFunc.set_numPoints(self.numPoints)
            resultFunc.setPoints(self.points)
            #numpy subtract float/int from array is how we want it
            resultFunc.setData( self.data - rightArg)
        else:
            print "SUBTRACTION of ", type(rightArg), " from dataFunction NOT DEFINED!!!"
            return

        return resultFunc
    
           # for i in range(numPointsMax):
                #print self.points[i]
                #print rightArg.points[i]
            #    print "no"
              
        
