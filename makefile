PROG		= main

# GNU fortran
FC 		= gfortran
# normal ver.
FFLAGS		= -O4 # -pg
# debug ver.
#FFLAGS		= -g -O0 -Wall -Wtabs -Wintrinsics-std -Wintrinsic-shadow -fbounds-check -O -Wuninitialized -fbacktrace

# Intel Fortran option
#FC		= ifort
# normal ver.
#FFLAGS		=
# IEEE 754 ver.
#FFLAGS		= -fp-model precise -fimf-arch-consistency=true 
# debug ver.
#FFLAGS		= -CB -traceback -g -check uninit -warn all -check all -std


COMMON_MOD 	= globvar.f90 func.f90 solver.f90
OBJCTS		= globvar.o func.o sub.o solver.o main.o

.SUFFIXES	:
.SUFFIXES	: .o .f90

.f90.o:
	${FC} -c $< ${FFLAGS} 

${PROG}	:	${OBJCTS}
	${FC} -o $@ ${OBJCTS} ${FFLAGS} >& err.d

${OBJCTS}: ${COMMON_MOD}

clean:
	rm -f ${PROG} err.d *.o *.mod *genmod.*

allclean:
	rm -f ${PROG} ${PROG}.exe err.d *.o *.mod *genmod.f90 *.mod0 solution.dat reshis.dat log.csv info.dat
