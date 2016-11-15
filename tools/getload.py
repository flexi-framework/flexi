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


nCoresMax=9999999
nCoresPerNode=24
tol=0.05

nDOF=nElems*(N+1)**3

print "============================================="
print "Number of DOF " + str(nDOF)


nCores=min(nCoresMax,nDOF/nDOFTarget)
nNodesTmp=floor(nCores/nCoresPerNode)
diff=ceil(nElems / nCores) - (nElems / nCores)
while diff>tol:
    nNodesTmp=nNodesTmp-1
    nCores=nNodesTmp*nCoresPerNode
    diff=ceil(nElems / nCores) - (nElems / nCores)

print "============================================="
print "Approaching from bottom:"
print "Suggested number of nodes " + str(nNodesTmp)
print "results in " + str(nElems / nCores) + " elements "
print "and " + str(nElems / nCores *(N+1)**3) + " DOF per core "
print "============================================="

nCores=min(nCoresMax,nDOF/nDOFTarget)
nNodesTmp=floor(nCores/nCoresPerNode)
diff=ceil(nElems / nCores) - (nElems / nCores)
while diff>tol:
    nNodesTmp=nNodesTmp+1
    nCores=nNodesTmp*nCoresPerNode
    diff=ceil(nElems / nCores) - (nElems / nCores)

print "Approaching from top:"
print "Suggested number of nodes " + str(nNodesTmp)
print "results in " + str(nElems / nCores) + " elements "
print "and " + str(nElems / nCores *(N+1)**3) + " DOF per core."
print "============================================="

