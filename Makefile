# Fortran Compiler
FC  = gfortran
# Compilation flags
F90FLAGS = -funderscoring -O0 -g -Wall -Wtabs -fbounds-check -ffpe-trap='zero'
# Linking with locally-installed LAPACK and BLAS libraries
# Mac OS
FLINKER = -framework Accelerate
# Linux
# FLINKER = -L/usr/local/lib -llapack -lblas

# Main compilation
all: clean elopph

# .f90: free-form Fortran 95/90 source files
%.o: %.f90
	$(FC) -c $(F90FLAGS) $< -o $@ # $< input files; $@ output files

%.o: %.f
	$(FC) -c $(F90FLAGS) $< -o $@

MODS = globvar.o
OBJS = libMath.o swntpar.o swntstruct.o
TUBE = libswntElec.o libswntPhon.o libswntOpt.f90

elopph: $(OBJS) $(MODS) $(TUBE)
	$(FC) $(OBJS) $(MODS) $(TUBE) \
	elopphtube.f90 -o elopphtube.out $(FLINKER)

clean:
	rm -f *.a *.exe *.mod *.o *.out *~ *.amk *.default
