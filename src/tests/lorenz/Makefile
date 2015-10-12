# FOODiE Makefile

INCLUDE = ../../lib
EXTERNAL = ../../../external/FLAP
LIB     = ../../lib

OBJS = lorenz.o type_lorenz.o

.SUFFIXES: .f90 .o

.f90.o:
	$(FC) -c $(FCFLAGS) -I$(INCLUDE) -I$(EXTERNAL) $<

lorenz: lorenz.o
	$(FC) $(OBJS) $(FCFLAGS) -L$(EXTERNAL) -lexternal -L$(LIB) -lfoodie -o $@

lorenz.o: lorenz.f90 type_lorenz.o
type_lorenz.o: type_lorenz.f90

.PHONY:
clean:
	rm -f lorenz *.o *.mod
