# Fortran codes for the BA-GMRES method preconditioned by the NR-SOR inner iterations

by Keiichi Morikuni and Ken Hayami

Latest update: July 21, 2020


## License

This software is currently released under the [GNU General Public License](http://www.gnu.org/copyleft/gpl.html).

If you use this code in research for publication, please cite the papers

1. Keiichi Morikuni and Ken Hayami, Inner-iteration Krylov subspace methods for least squares problems, SIAM Journal on Matrix Analysis and Applications, Volume 34, Number 1, pages 1-22, 2013. DOI: [10.1137/110828472](https://doi.org/10.1137/110828472)
2. Keiichi Morikuni and Ken Hayami, Convergence of inner-iteration GMRES methods for rank-deficient least squares problems, SIAM Journal on Matrix Analysis and Applications, Volume 36, Number 1, pages 225-250, 2015. DOI: [10.1137/130946009](https://doi.org/10.1137/130946009)

For commercial use, please make a contact to
Keiichi Morikuni [morikuni.keiichi.fw@u.tsukuba.ac.jp](mailto:morikuni.keiichi.fw@u.tsukuba.ac.jp)



## Organization

The codes consist of the following files:

* main.f90
* solver.f90
* sub.f90
* func.f90
* floating_point_kinds.f90
* plot.plt
* makefile
* readme.md
* test.sh

A test matrix called RANDL7 in the compressed column storage (CCS) format
is given in directory RANDL7. The values of the NR-SOR inner-iteration
parameters can be automatically tuned at each restart or specified by you.

## Usage

To simply run the program, execute the following:

```
$ sh test.sh
```

To compile the codes, execute the following:

```
$ make
```

Change the Fortran compiler given in the Makefile file if necessary.

To simply run the program with the default values of parameters on the test matrix RANDL7, execute the following:

```
$ ./main --fi=RANDL7/
```
Then the program outputs the approximate solution data solution.dat, the result data info.dat, and the relative residual norm history data reshis.dat.

To run the program with specific values of parameters on the test matrix RANDL7, execute the following:

```
$ ./main --at=1 --nin=50 --omg=1.0 --tol=1.0e-8 --omax=800 --rmax=0 —output_mode=0 --fi=RANDL7/
```

Some specific data is output in log.csv.


- `--at=`
This option enables to automatically determine the values of the number of inner
iterations and the relaxation parameter. A value for this option must be provided;
possible values are
-- `0`: turn off the automatic parameter tuning
-- `1`: turn on it.

- `--nin=`: Nonnegative integer
This option determines the  
	- number of inner iterations for the setting for `--a=0`
	- maximum number of inner iterations for `--at=n`, n > 0 and the actual number of inner iterations are automatically determined.

- `--omg=`: This option determines the value of the relaxation parameter for `--at=0`; otherwise the value provided is used as the initial value of the relaxation parameter for the automatic parameter tuning.

- `--tol=`: This option determines the threshold in terms of the relative residual norm for terminating the iterations. Typically, the value is less than one.

- `--omax=`: This option determines the maximum number of outer iterations.

- `--rmax=`: This option determines the number of restart cycles. The restart is turned off for `—rmax=0`.

- `--output_mode=`: This option enables a detailed output display.  
	- `--output_mode=0`: Turn off the detailed output display  
	- `--output_mode=1`: Turn on the detailed output display.

- `--fi=`: This option determines the directory name in which the matrix data used is contained. 
The directory name must be the relative one.

## Contacts

Please provide feedback to [morikuni.keiichi.fw@u.tsukuba.ac.jp](mailto:morikuni.keiichi.fw@u.tsukuba.ac.jp) if you have any questions or suggestions.

Keiichi Morikuni, Ph.D.  

Affiliation: Faculty of Engineering, Information and Systems, University of Tsukuba  
Postal Address: 1-1-1 Tennodai, Tsukuba, Ibaraki 305-8573, Japan

Homepage URL: [http://researchmap.jp/KeiichiMorikuni/](http://researchmap.jp/KeiichiMorikuni/)

## Support

The Graduate University for Advanced Studies (SOKENDAI), Shonan Village, Hayama, Kanagawa 240-0193 Japan


### References

1. Keiichi Morikuni and Ken Hayami, Inner-iteration Krylov subspace methods for least squares problems, SIAM Journal on Matrix Analysis and Applications, Volume 34, Number 1, pages 1-22, 2013. DOI: [10.1137/110828472](https://doi.org/10.1137/110828472)
2. Keiichi Morikuni and Ken Hayami, Convergence of inner-iteration GMRES methods for rank-deficient least squares problems, SIAM Journal on Matrix Analysis and Applications, Volume 36, Number 1, pages 225-250, 2015. DOI: [10.1137/130946009](https://doi.org/10.1137/130946009)

