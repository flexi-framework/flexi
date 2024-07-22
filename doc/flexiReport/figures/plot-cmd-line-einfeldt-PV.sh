#!/bin/sh

python3 ../tools/xyPlotting/multiXY.py -t ./plotInputFiles/ParaVExported/flexi-PV-einfeldt.tsv -t ./plotInputFiles/exact/exact-PV-einfeldt123.dat -o "einfeldt-PV.pdf"
python3 ../tools/norm/solnNorm.py -t ./plotInputFiles/ParaVExported/flexi-PV-einfeldt.tsv -t ./plotInputFiles/exact/exact-PV-einfeldt123.dat -o "einfeldt-L2Norm-PV.dat"
