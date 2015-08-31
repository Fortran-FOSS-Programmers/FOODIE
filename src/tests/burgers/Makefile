# FOODiE Makefile

INCLUDE = ../../lib
LIB     = ../../lib

.SUFFIXES: .f90 .o

.f90.o:
	$(FC) -c $(FCFLAGS) -I$(INCLUDE) $<

burgers: burgers.o
	$(FC) burgers.f90 type_burgers.o $(LIB)/*.o $(FCFLAGS) -I$(INCLUDE) -o $@

burgers.o: burgers.f90 type_burgers.o
type_burgers.o: type_burgers.f90

.PHONY:
clean:
	rm -f burgers *.o *.mod
