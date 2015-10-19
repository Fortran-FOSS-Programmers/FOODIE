# FOODIE Makefile

INCLUDE  = ../../lib
EXTERNAL = ../../../external/FLAP
LIB      = ../../lib

OBJS = burgers.o type_burgers.o

.SUFFIXES: .f90 .o

.f90.o:
	$(FC) -c $(FCFLAGS) -I$(INCLUDE) -I$(EXTERNAL) $<

burgers: burgers.o
	$(FC) $(OBJS) $(FCFLAGS) -L$(EXTERNAL) -lexternal -L$(LIB) -lfoodie -o $@

burgers.o: burgers.f90 type_burgers.o
type_burgers.o: type_burgers.f90

.PHONY:
clean:
	rm -f burgers *.o *.mod
