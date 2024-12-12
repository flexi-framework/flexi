#!/bin/sh

python3 ../tools/xyPlotting/CVmultiXY.py -t ./plotInputFiles/ParaVExported/flexi-CV-sedov.tsv -t ./plotInputFiles/exact/exact-CV-sedov.dat -o "sedov-CV.pdf"
python3 ../tools/norm/solnNorm.py -t ./plotInputFiles/ParaVExported/flexi-CV-sedov.tsv -t ./plotInputFiles/exact/exact-CV-sedov.dat -o "sedov-L2Norm-CV.dat"

