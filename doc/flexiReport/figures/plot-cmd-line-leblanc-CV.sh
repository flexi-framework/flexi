#!/bin/sh

python3 ../tools/xyPlotting/CVmultiXY.py -t ./plotInputFiles/ParaVExported/flexi-CV-leblanc.tsv -t ./plotInputFiles/exact/exact-CV-leblanc.dat -o "leblanc-CV.pdf"
python3 ../tools/norm/solnNorm.py -t ./plotInputFiles/ParaVExported/flexi-CV-leblanc.tsv -t ./plotInputFiles/exact/exact-CV-leblanc.dat -o "leblanc-L2Norm-CV.dat"

