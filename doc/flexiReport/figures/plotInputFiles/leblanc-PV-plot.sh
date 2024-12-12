#!/bin/sh

python3 /home/rod/Research/flexi/doc/flexiReport/tools/xyPlotting/multiXY.py -t ./plotInputFiles/ParaVExported/flexi-PV-einfeldt.tsv -t ./plotInputFiles/exact/exact-PV-einfeldt123.dat -o "PVeinfeldt-multi.pdf"
python3 ../tools/xyPlotting/solnNorm.py -t ./plotInputFiles/ParaVExported/flexi-PV-einfeldt.tsv -t ./plotInputFiles/exact/exact-PV-einfeldt123.dat -o "PVeinfeldt-L2Norm.dat"
