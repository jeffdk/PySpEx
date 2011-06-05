#!/bin/bash


for i in 0 1 2 
do
    portionList={'DensestPoint.dat'}
    postFix=pd1/Run/
    base=/home/jeff/tovStars/stdTOVFiles/tenMSsL
    scp zwicky:${base}${i}fL${i}${postFix}DensestPoint.dat ./Lev${i}_DensestPoint.dat
done