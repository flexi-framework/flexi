#!/bin/sh

python3 ../tools/xyPlotting/multiXY.py -t ./plotInputFiles/ParaVExported/flexi-sod2-BC2.tsv -t ./plotInputFiles/exact/exact-PV-sod2-BC9.dat -o "sod2-BC2-PV.pdf"
python3 ../tools/norm/solnNorm.py -t ./plotInputFiles/ParaVExported/flexi-sod2-BC2.tsv -t ./plotInputFiles/exact/exact-PV-sod2-BC9.dat -o "sod2-L2Norm-BC2-PV.dat"
