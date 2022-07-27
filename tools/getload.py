#! /usr/bin/env python
#************************************************************************************
#
# Description:  This script will compute a the number of cores for a
#               given number of elements, polynomial degree and desired load/core
#               in a sense that for equal element distribution there will be no
#               single cores having more elements then the average, only less.
#
#************************************************************************************

nElems        = 86712
N             = 7
nDOFTarget    = 6000.0
dim           = 3

nCoresMax     = 5000.0
nCoresPerNode = 24

# END OF USER INPUT BLOCK ***********************************************************

from numpy     import *
from termcolor import colored

nDOF = nElems*(N+1)**dim

print("┌────────────────────────────────────────────")
print("│ Number of DOF: {}".format(nDOF))

nCores             = min(nCoresMax,nDOF/nDOFTarget)
nElemsPerCoreIdeal = ceil(nElems / nCores)
nElemsPerNodeIdeal = nElemsPerCoreIdeal*nCoresPerNode
nNodes             = int(ceil(nElems/(nElemsPerNodeIdeal)))
nCores             = nNodes*nCoresPerNode
nElemsPerCore      = nElems / nCores
Imbalance          = (nElemsPerCoreIdeal - nElemsPerCore) / nElemsPerCoreIdeal

print("├────────────────────────────────────────────")
print("│ " + colored("Approaching from bottom:",'white',attrs=['bold']))
print("│ Suggested number of nodes : {:9d} ({:d} cores)".format(nNodes,nCores))
print("├── resulting elements/core : {:9.3f}".format(nElems / nCores))
print("├── resulting DOF     /core : {:9.3f}".format(nElems / nCores *(N+1)**dim))
print("├── efficiency loss         : {:9.3f}%".format(100*Imbalance))

nCores             = min(nCoresMax,nDOF/nDOFTarget)
nElemsPerCoreIdeal = floor(nElems / nCores)
nElemsPerNodeIdeal = nElemsPerCoreIdeal*nCoresPerNode
nNodes             = int(ceil(nElems/(nElemsPerNodeIdeal)))
nCores             = nNodes*nCoresPerNode
nElemsPerCore      = nElems / nCores
Imbalance          = (nElemsPerCoreIdeal - nElemsPerCore) / nElemsPerCoreIdeal

print("├────────────────────────────────────────────")
print("│ " + colored("Approaching from top:",'white',attrs=['bold']))
print("│ Suggested number of nodes : {:9d} ({:d} cores)".format(nNodes,nCores))
print("├── resulting elements/core : {:9.3f}".format(nElems / nCores))
print("├── resulting DOF     /core : {:9.3f}".format(nElems / nCores *(N+1)**dim))
print("├── efficiency loss         : {:9.3f}%".format(100*Imbalance))
print("└────────────────────────────────────────────")
