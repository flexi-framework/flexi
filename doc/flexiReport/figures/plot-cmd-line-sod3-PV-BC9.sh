#!/bin/sh

python3 ../tools/xyPlotting/multiXY.py -t ./plotInputFiles/ParaVExported/flexi-sod3-BC9-best.tsv -t ./plotInputFiles/exact/exact-PV-sod3-BC9-mod.dat -o "sod3-BC9-PV.pdf"
python3 ../tools/norm/solnNorm.py -t ./plotInputFiles/ParaVExported/flexi-sod3-BC9-best.tsv -t ./plotInputFiles/exact/exact-PV-sod3-BC9-mod.dat -o "sod3-L2Norm-BC9-PV.dat"
