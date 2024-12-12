#!/bin/sh

python3 ../../../tools/ToroExact/toro_exact.py -n 1001 -d nodal -g 1.4 -p user --left "[1.0,-2.0,0.4]" --right "[1.0,2.0,0.4]" -b "[0.0,1.0]" -x 0.5 -t 0.05 -VT "PV" -o "exact-PV-einfeldt123.dat"
