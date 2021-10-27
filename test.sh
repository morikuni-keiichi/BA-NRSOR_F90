#!/bin/bash
make
./main --at=1 --nin=50 --omg=1.0 --tol=1.0e-8 --omax=800 --rmax=0 -v --directory=RANDL7/
gnuplot plot.plt