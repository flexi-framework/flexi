#!/bin/sh

python3 ../tools/xyPlotting/multiXY.py -t ./plotInputFiles/ParaVExported/flexi-sod1-g-5o3-t-178o1000-x0-3o4.tsv -t ./plotInputFiles/exact/exact-PV-sod1-BC9.dat -o "sod1-BC9-PV.pdf"
python3 ../tools/norm/solnNorm.py -t ./plotInputFiles/ParaVExported/flexi-sod1-g-5o3-t-178o1000-x0-3o4.tsv -t ./plotInputFiles/exact/exact-PV-sod1-BC9.dat -o "sod1-BC9-L2norm-PV.dat"
