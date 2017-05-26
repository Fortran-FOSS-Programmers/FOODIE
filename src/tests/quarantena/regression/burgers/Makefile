# FOODIE Makefile

OBJS = burgers.o type_burgers.o

.SUFFIXES: .f90 .o

.f90.o:
	$(FC) -c $(FCFLAGS) -I$(INCLUDE) -I$(FLAPDIR) -I$(PYPLOTDIR)/src $<

burgers: burgers.o
	$(FC) $(OBJS) $(FCFLAGS) -L$(FLAPDIR) -lflap -L$(PYPLOTDIR) -lpyplot -L$(LIB) -lfoodie -o $@

burgers.o: burgers.f90 type_burgers.o
type_burgers.o: type_burgers.f90

.PHONY:
clean:
	rm -f burgers *.o *.mod
