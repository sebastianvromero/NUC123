FORTRAN = gfortran
FLAGS = -fdefault-real-8
F77FLAGS = $(FLAGS) -ffixed-line-length-none
F90FLAGS = $(FLAGS) -fimplicit-none -cpp

modules = bessel.o phys_parameters.o comp_parameters.o variables.o \
          scalpot.o universe.o rate_integrator.o reactions.o quadpack.o

bbn.exe: bbn.f90 driver.o $(modules)
	$(FORTRAN) $(F90FLAGS) bbn.f90 -o bbn.exe driver.o $(modules)

nuc123.exe: NUC123.f driver.o $(modules)
	$(FORTRAN) $(F77FLAGS) NUC123.f -o nuc123.exe driver.o $(modules)

driver.o: driver.f $(modules)
	$(FORTRAN) $(F77FLAGS) -c driver.f quadpack.f

# Default rule for modules (if they have no dependencies)
%.o: %.f90
	$(FORTRAN) $(F90FLAGS) -c $<

%.o: %.F90
	$(FORTRAN) $(F90FLAGS) -c $<


.PHONY: clean
clean:
	-rm -f bbn.exe nuc123.exe *.o *.mod

