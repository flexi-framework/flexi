#! /usr/bin/env python
#************************************************************************************
#
# Author:       Thomas Bolemann
# Institution:  Inst. of Aero- and Gasdynamics, University of Stuttgart
# Date:         07.07.2016
#
# Description:  This script will compute a the number of cores for a
#               given number of elements, polynomial degree and desired load/core 
#               in a sense that for equal element distribution there will be no
#               single cores having more elements then the average, only less.
# 
#************************************************************************************

from numpy import *

nElems=86712
N=7
nDOFTarget=6000.0
dim=3


nCoresMax=5000.0
nCoresPerNode=24

nDOF=nElems*(N+1)**dim

print "============================================="
print "Number of DOF " + str(nDOF)

nCores=min(nCoresMax,nDOF/nDOFTarget)
nElemsPerCoreIdeal=ceil(nElems / nCores)
nElemsPerNodeIdeal=nElemsPerCoreIdeal*nCoresPerNode
nNodes=ceil(nElems/(nElemsPerNodeIdeal))
nCores=nNodes*nCoresPerNode
nElemsPerCore = nElems / nCores
Imbalance=(nElemsPerCoreIdeal - nElemsPerCore) / nElemsPerCoreIdeal

print "============================================="
print "Approaching from bottom:"
print "Suggested number of nodes " + str(int(nNodes)) + " (" + str(int(nCores)) + " cores)"
print "results in " + str(nElems / nCores) + " elements "
print "and " + str(nElems / nCores *(N+1)**dim) + " DOF per core "
print "Efficiency loss: " + str(100*Imbalance) + " %"
print "============================================="

nCores=min(nCoresMax,nDOF/nDOFTarget)
nElemsPerCoreIdeal=floor(nElems / nCores)
nElemsPerNodeIdeal=nElemsPerCoreIdeal*nCoresPerNode
nNodes=ceil(nElems/(nElemsPerNodeIdeal))
nCores=nNodes*nCoresPerNode
nElemsPerCore = nElems / nCores
Imbalance=(nElemsPerCoreIdeal - nElemsPerCore) / nElemsPerCoreIdeal

print "Approaching from top:"
print "Suggested number of nodes " + str(int(nNodes)) + " (" + str(int(nCores)) + " cores)"
print "results in " + str(nElems / nCores) + " elements "
print "and " + str(nElems / nCores *(N+1)**dim) + " DOF per core "
print "Efficiency loss: " + str(100*Imbalance) + " %"
print "============================================="

