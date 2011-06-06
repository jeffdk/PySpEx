#!/usr/bin/env python
#
#
from numpy import *

class DatSet:
    something="nothing"

class Simulation:

    #Resolutions; should be tuples of points in a linear direction
    #First item: GrResolution
    #second item: HydroResolution
    # GrResolution could be a set of Npoints (Ex.  Nr, Ntheta, Nphi
    #     for spherical shells)
    
    ##The following 4 must have the same length
    Levs=[]
    DatSets=[]          #This is the filename prefix
    DatSetsLabels={}
    Resolutions=[]

    Legends=[]

    Directory='./'

    Name='sim'

    timeConversionFactor = 4.926e-6
    densityConversionFactor= 6.173e17
    

    
    def __init__():
        Levs=['Lev0,Lev1,Lev2']
        DatSets=['DensestPoint']          #This is the filename prefix
        DatSetsLabels={'DensestPoint':1}
        Resolutions=[45,45+6*3,45+2*6*3]
        
        Legends=[]
        
        Directory='./data/'
        
        Name='sim'


