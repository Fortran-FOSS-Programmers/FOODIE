# FOODIE Makefile

INCLUDE  = ../../lib
EXTERNAL = ../../../external/FLAP
LIB      = ../../lib

OBJS = oscillation.o type_oscillation.o

.SUFFIXES: .f90 .o

.f90.o:
	$(FC) -c $(FCFLAGS) -I$(INCLUDE) -I$(EXTERNAL) $<

oscillation: oscillation.o
	$(FC) $(OBJS) $(FCFLAGS) -L$(EXTERNAL) -lexternal -L$(LIB) -lfoodie -o $@

oscillation.o: oscillation.f90 type_oscillation.o
type_oscillation.o: type_oscillation.f90

.PHONY:
clean:
	rm -f oscillation *.o *.mod
