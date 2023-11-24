#!/bin/sh

python3 ../../../tools/ToroExact/toro_exact.py -n 1001 -d nodal -g 1.6666667 -p user --left "[1.0,0.0,66.6667]" --right "[0.2,0.0,0.06667]" -b "[0.0,1.0]" -x 0.3 -t 0.01812 -VT "PV" -o "exact-PV-sod3-BC9-mod.dat"
#python3 /usr/local/packages/ToroExact/toro_exact.py -n 101 -d nodal -g 1.6666667 -p user --left "[1.0,0.0,66.6667]" --right "[0.2,0.0,0.06667]" -b "[0.0,1.0]" -x 0.3 -t 0.01812 -VT "CV" 
#python3 /usr/local/packages/ToroExact/toro_exact.py -n 101 -d nodal -g 1.6666667 -p user --left "[1.0,0.0,66.6667]" --right "[0.2,0.0,0.06667]" -b "[0.0,1.0]" -x 0.3 -t 0.01812 -VT "CV" -o "exact-CV-sod3-BC9-mod-101.dat"

