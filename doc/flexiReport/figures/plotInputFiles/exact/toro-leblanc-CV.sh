#!/bin/sh

python3 ../../../tools/ToroExact/toro_exact.py -n 1001 -d nodal -g 1.6666667 -p user --left "[1.0,0.0,0.066666667]" --right "[0.001,0.0,0.000000000066666667]" -b "[0.0,9.0]" -x 3.0 -t 6.0 -VT "CV" -o "exact-CV-leblanc.dat"
