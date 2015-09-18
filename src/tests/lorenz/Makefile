# FOODiE Makefile

INCLUDE = ../../lib
LIB     = ../../lib

.SUFFIXES: .f90 .o

.f90.o:
	$(FC) -c $(FCFLAGS) -I$(INCLUDE) $<

lorenz: lorenz.o
	$(FC) lorenz.f90 type_lorenz.o $(FCFLAGS) -I$(INCLUDE) -L$(LIB) -lfoodie -o $@

lorenz.o: lorenz.f90 type_lorenz.o
type_lorenz.o: type_lorenz.f90

.PHONY:
clean:
	rm -f lorenz *.o *.mod
