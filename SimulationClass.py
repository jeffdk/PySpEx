#!/usr/bin/env python
#
#
from numpy import *

from dataFunction import *

class DatSet:

    title=''
    simulationName =''

    Data=dataFunction()

    def __init__(self, filename):
        
        filehandle=open(filename, 'r')

        
        

class Simulation:

    #Resolutions; should be tuples of points in a linear direction
    #First item: GrResolution
    #second item: HydroResolution
    # GrResolution could be a set of Npoints (Ex.  Nr, Ntheta, Nphi
    #     for spherical shells)
    Levs=[]
    
    DataSets=[]


    
