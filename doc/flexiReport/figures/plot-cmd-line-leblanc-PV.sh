#!/bin/sh

python3 ../tools/xyPlotting/multiXY.py -t ./plotInputFiles/ParaVExported/flexi-PV-leblanc.csv -t ./plotInputFiles/exact/exact-PV-leblanc.dat -o "leblanc-PV.pdf"
python3 ../tools/norm/solnNorm.py -t ./plotInputFiles/ParaVExported/flexi-PV-leblanc.csv -t ./plotInputFiles/exact/exact-PV-leblanc.dat -o "leblanc-L2Norm-PV.dat"

