# Fortran Compiler
FC  = gfortran
# Compilation flags
F90FLAGS = -funderscoring -O0 -g -Wall -Wtabs -fbounds-check
# Linking with locally-installed LAPACK and BLAS libraries
# Mac OS
FLINKER = -framework Accelerate
# Linux
# FLINKER = -L/usr/local/lib -llapack -lblas

all: clean cntabsorpt

# .f90: free-form Fortran 95/90 source files
%.o: %.f90
	$(FC) -c $(F90FLAGS) $< -o $@ # $< input files; $@ output files

%.o: %.f
	$(FC) -c $(F90FLAGS) $< -o $@

OBJS = libMath.o swntpar.o swntstruct.o 
TUBE = libswntElec.o libswntOpt.o libDrudeOpt.o libswntSTB.o

# Main compilation
cntabsorpt: $(OBJS) $(TUBE)
	$(FC) $(OBJS) $(TUBE) \
	$(F90FLAGS) cntabsorpt.f90 -o cntabsorpt.out $(FLINKER)

clean:
	rm -f *.a *.exe *.mod *.o *.out *~ *.amk *.default
