# FOODIE Makefile

OBJS = euler-1D.o type_euler-1D.o

.SUFFIXES: .f90 .o

.f90.o:
	$(FC) -c $(FCFLAGS) -I$(INCLUDE) -I$(FLAPDIR) -I$(PYPLOTDIR)/src -I$(WENOOFDIR)/src/lib $<

euler-1D: euler-1D.o
	$(FC) $(OBJS) $(FCFLAGS) -L$(FLAPDIR) -lflap -L$(PYPLOTDIR) -lpyplot -L$(WENOOFDIR) -lwenoof -L$(LIB) -lfoodie -o $@

euler-1D.o: euler-1D.f90 type_euler-1D.o
type_euler-1D.o: type_euler-1D.f90

.PHONY:
clean:
	rm -f euler-1D *.o *.mod
