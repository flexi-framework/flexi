#!/bin/sh

python3 ../../../tools/ToroExact/toro_exact.py -n 1001 -d nodal -g 1.6666667 -p user --left "[1.0,0.0,1.0]" --right "[0.125,0.0,0.1]" -b "[0.0,1.0]" -x 0.5 -t 0.2746 -o "exact-PV-sod2-BC9.dat"


