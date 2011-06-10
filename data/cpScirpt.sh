#!/bin/bash


for i in 0 1 2 
do
    #portionList= ( 'DensestPoint.dat' 'GhCe.dat' 'RestMass.dat' )
    
    postFix=pd1/Run/
    base=/home/jeff/tovStars/stdTOVFiles/longerGaugesL
    scp zwicky:${base}${i}fL${i}${postFix}DensestPoint.dat ./Lev${i}_DensestPoint.dat
    scp zwicky:${base}${i}fL${i}${postFix}RestMass.dat ./Lev${i}_RestMass.dat
    scp zwicky:${base}${i}fL${i}${postFix}/Constraints/GhCe.dat ./Lev${i}_GhCe.dat
done