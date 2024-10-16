#!/bin/sh

python3 ../../../tools/ToroExact/toro_exact.py -n 101 -d nodal -g 1.6666667 -p user --left "[1.0,0.0,1.0]" --right "[1.0,0.0,0.1]" -b "[0.0,1.0]" -x 0.75 -t 0.1780 -VT "CV" -o "exact-CV-sod1-BC9.dat"
