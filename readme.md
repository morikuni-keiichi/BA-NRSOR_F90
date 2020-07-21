# Fortran codes for the BA-GMRES method preconditioned by the NR-SOR inner iterations

by Keiichi Morikuni and Ken Hayami

Latest update: July 21, 2020


This software is currently released under the GNU General Public License
http://www.gnu.org/copyleft/gpl.html

## Copyright 

2013-2020 Keiichi Morikuni [morikuni.keiichi.fw@u.tsukuba.ac.jp](mailto:morikuni.keiichi.fw@u.tsukuba.ac.jp)
URL: http://researchmap.jp/KeiichiMorikuni/

## Maintenace

Keiichi Morikuni <[morikuni.keiichi.fw@u.tsukuba.ac.jp](mailto:morikuni.keiichi.fw@u.tsukuba.ac.jp)>
Assistant Professor, Ph.D.
Faculty of Engineering, Information and Systems
University of Tsukuba
1-1-1 Tennodai, Tsukuba, Ibaraki 305-8573, Japan

## Support

**The Graduate University for Advanced Studies (SOKENDAI)**
Shonan Village, Hayama, Kanagawa 240-0193 Japan

If you use this code in research for publication, please cite the papers

1. Keiichi Morikuni and Ken Hayami, Inner-iteration Krylov subspace methods for least squares problems, SIAM Journal on Matrix Analysis and Applications, Volume 34, Number 1, pages 1-22, 2013. DOI: 10.1137/110828472
2. Keiichi Morikuni and Ken Hayami, Convergence of inner-iteration GMRES methods for rank-deficient least squares problems, SIAM Journal on Matrix Analysis and Applications, Volume 36, Number 1, pages 225-250, 2015. DOI: 10.1137/130946009

For the commercial use, please make a contact to
Keiichi Morikuni [morikuni.keiichi.fw@u.tsukuba.ac.jp](mailto:morikuni.keiichi.fw@u.tsukuba.ac.jp)

BA-GMRES preconditioned by NR-SOR inner iterations involves the following files:

main.f90
solver.f90
sub.f90
func.f90
globvar.f90 
plot.plt
makefile
readme.md

A test matrix called RANDL7 in the compressed column storage (CCS) format
is given in directory RANDL7. The values of the NR-SOR inner iteration
parameters can be automatically tuned at each restart or specified by you.

To compile the codes and run the program, proceed as follows:

$ make
$ ./main --at=1 --nin=50 --omg=1.0 --tol=1.0e-8 --omax=800 --rmax=0 —output_mode=0 —fi=RANDL7/

Then the program outputs the approximate solution data solution.dat,
the result data info.dat, and the relative residual norm history data reshis.dat.
In addition, some specific data is output in log.csv.

Here,
--at=
This option enables to automatically determine the values of the number of inner
iterations and the relaxation parameter. A value for this option must be provided;
possible values are
0: turn off the automatic parameter tuning
1: turn on it

--nin=
This option determines the

- number of inner iterations for the setting —at=0 and

- maximum number of inner iterations for the setting --at=n, n>0 and the actual number of inner iterations are automatically determined,

which must be a nonnegative integer.

--omg=
This option determines the value of the relaxation parameter for --at=0; otherwise the value provided is used as the initial value of the relaxation parameter for the automatic parameter tuning

--tol=
This option determines the threshold in terms of the relative residual norm for terminating the iterations. Typically, the value is less than one.

--omax=
This option determines the maximum number of outer iterations.

--rmax=
This option determines the restart cycle. The restart is turned off for —rmax=0.

--output_mode=
This option enables a detailed output display.

0: turn off the detailed output display.
1: turn on the detailed output display.

--fi=
This option determines the directory name in which the matrix data used is contained.
The directory name must be the relative one.

Please provide feedback if you have any questions or suggestions.
[morikuni.keiichi.fw@u.tsukuba.ac.jp](mailto:morikuni.keiichi.fw@u.tsukuba.ac.jp)


### References

[1] Keiichi Morikuni and Ken Hayami,
Inner-iteration Krylov subspace methods for least squares problems,
SIAM Journal on Matrix Analysis and Applications, Volume 34, Number 1, pages 1-22, 2013.
DOI: 10.1137/110828472
[2] Keiichi Morikuni and Ken Hayami,
Convergence of inner-iteration GMRES methods for rank-deficient least squares problems,
SIAM Journal on Matrix Analysis and Applications, Volume 36, Number 1, pages 225-250, 2015.
DOI: 10.1137/130946009
