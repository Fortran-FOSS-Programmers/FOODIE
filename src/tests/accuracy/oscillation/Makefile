# FOODIE Makefile

OBJS = oscillation.o oscillation_t.o

.SUFFIXES: .f90 .o

.f90.o:
	$(FC) -c $(FCFLAGS) -I$(INCLUDE) -I$(FLAPDIR) -I$(PYPLOTDIR)/src $<

oscillation: oscillation.o
	$(FC) $(OBJS) $(FCFLAGS) -L$(FLAPDIR) -lflap -L$(PYPLOTDIR) -lpyplot -L$(LIB) -lfoodie -o $@

oscillation.o: oscillation.f90 oscillation_t.o
oscillation_t.o: oscillation_t.f90

.PHONY:
clean:
	rm -f oscillation *.o *.mod
