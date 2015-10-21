# FOODIE Makefile

INCLUDE = ../../lib
LIB     = ../../lib

.SUFFIXES: .f90 .o

.f90.o:
	$(FC) -c $(FCFLAGS) -I$(INCLUDE) $<

euler-1D: euler-1D.o
	$(FC) euler-1D.f90 type_euler-1D.o $(FCFLAGS) -I$(INCLUDE) -L$(LIB) -lfoodie -o $@

euler-1D.o: euler-1D.f90 type_euler-1D.o
type_euler-1D.o: type_euler-1D.f90

.PHONY:
clean:
	rm -f euler-1D *.o *.mod
