#!/bin/bash
make
./main --at=1 --nin=50 --omg=1.0d0 --tol=1e-8 --omax=800 --rmax=0 --output_mode=0 --fi=RANDL7/
gnuplot plot.plt