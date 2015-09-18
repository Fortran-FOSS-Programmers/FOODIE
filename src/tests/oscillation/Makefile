# FOODiE Makefile

INCLUDE = ../../lib
LIB     = ../../lib

.SUFFIXES: .f90 .o

.f90.o:
	$(FC) -c $(FCFLAGS) -I$(INCLUDE) $<

oscillation: oscillation.o
	$(FC) oscillation.f90 type_oscillation.o $(FCFLAGS) -I$(INCLUDE) -L$(LIB) -lfoodie -o $@

oscillation.o: oscillation.f90 type_oscillation.o
type_oscillation.o: type_oscillation.f90

.PHONY:
clean:
	rm -f oscillation *.o *.mod
