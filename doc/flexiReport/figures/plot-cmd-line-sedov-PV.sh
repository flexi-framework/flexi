#!/bin/sh

python3 ../tools/xyPlotting/multiXY.py -t ./plotInputFiles/ParaVExported/flexi-PV-sedov.tsv -t ./plotInputFiles/exact/exact-PV-sedov.dat -o "sedov-PV.pdf"
python3 ../tools/norm/solnNorm.py -t ./plotInputFiles/ParaVExported/flexi-PV-sedov.tsv -t ./plotInputFiles/exact/exact-PV-sedov.dat -o "sedov-L2Norm-PV.dat"

