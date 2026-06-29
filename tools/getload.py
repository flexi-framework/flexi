#! /usr/bin/env python
# ************************************************************************************
#
# Description:  This script will compute a the number of cores for a
#               given number of elements, polynomial degree and desired load/core
#               in a sense that for equal element distribution there will be no
#               single cores having more elements then the average, only less.
#
# ************************************************************************************
import numpy as np
from termcolor import colored

# ************************************************************************************

nElems        = 86712
N             = 7
nDOFTarget    = 6000.0
dim           = 3

nCoresMax     = 5000.0
nCoresPerNode = 24

# END OF USER INPUT BLOCK ***********************************************************

nDOF = nElems*(N+1)**dim

print("┌────────────────────────────────────────────")
print(f"│ Number of DOF: {nDOF}")

nCores             = min(nCoresMax, nDOF/nDOFTarget)
nElemsPerCoreIdeal = np.ceil(nElems / nCores)
nElemsPerNodeIdeal = nElemsPerCoreIdeal*nCoresPerNode
nNodes             = int(np.ceil(nElems/(nElemsPerNodeIdeal)))
nCores             = nNodes*nCoresPerNode
nElemsPerCore      = nElems / nCores
Imbalance          = (nElemsPerCoreIdeal - nElemsPerCore) / nElemsPerCoreIdeal

print("├────────────────────────────────────────────")
print("│ " + colored("Approaching from bottom:", 'white', attrs=['bold']))
print(f"│ Suggested number of nodes : {nNodes:9d} ({nCores:d} cores)")
print(f"├── resulting elements/core : {nElems / nCores:9.3f}")
print(f"├── resulting DOF     /core : {nElems / nCores *(N+1)**dim:9.3f}")
print(f"├── efficiency loss         : {100*Imbalance:9.3f}%")

nCores             = min(nCoresMax, nDOF/nDOFTarget)
nElemsPerCoreIdeal = np.floor(nElems / nCores)
nElemsPerNodeIdeal = nElemsPerCoreIdeal*nCoresPerNode
nNodes             = int(np.ceil(nElems/(nElemsPerNodeIdeal)))
nCores             = nNodes*nCoresPerNode
nElemsPerCore      = nElems / nCores
Imbalance          = (nElemsPerCoreIdeal - nElemsPerCore) / nElemsPerCoreIdeal

print("├────────────────────────────────────────────")
print("│ " + colored("Approaching from top:", 'white', attrs=['bold']))
print(f"│ Suggested number of nodes : {nNodes:9d} ({nCores:d} cores)")
print(f"├── resulting elements/core : {nElems / nCores:9.3f}")
print(f"├── resulting DOF     /core : {nElems / nCores *(N+1)**dim:9.3f}")
print(f"├── efficiency loss         : {100*Imbalance:9.3f}%")
print("└────────────────────────────────────────────")
