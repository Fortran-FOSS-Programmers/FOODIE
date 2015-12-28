# WENOOF Makefile

OBJS = IR_Precision.o wenoof.o type_weno_interpolator.o type_weno_interpolator_upwind.o 

.SUFFIXES: .f90 .o

.f90.o:
	$(FC) -c $(FCFLAGS) $<

wenoof.o: wenoof.f90 IR_Precision.o type_weno_interpolator.o type_weno_interpolator_upwind.o
type_weno_interpolator.o: type_weno_interpolator.f90 IR_Precision.o
type_weno_interpolator_upwind.o: type_weno_interpolator_upwind.f90 type_weno_interpolator.o IR_Precision.o
IR_Precision.o: IR_Precision.f90

.PHONY:
clean:
	rm -f *.o *.mod
