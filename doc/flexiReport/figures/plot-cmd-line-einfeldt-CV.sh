#!/bin/sh

python3 ../tools/xyPlotting/CVmultiXY.py -t ./plotInputFiles/ParaVExported/flexi-CV-einfeldt.tsv -t ./plotInputFiles/exact/exact-CV-einfeldt123.dat -o "einfeldt-CV.pdf"
python3 ../tools/norm/solnNorm.py -t ./plotInputFiles/ParaVExported/flexi-CV-einfeldt.tsv -t ./plotInputFiles/exact/exact-CV-einfeldt123.dat -o "einfeldt-L2Norm-CV.dat"

