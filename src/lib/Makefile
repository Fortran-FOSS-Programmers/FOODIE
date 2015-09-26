# FOODiE Makefile

.SUFFIXES: .F90 .f90 .o
.PHONY: all foodie clean

OBJS = type_integrand.o \
       foodie_integrator_adams_bashforth.o \
       foodie_integrator_euler_explicit.o \
       foodie_integrator_leapfrog.o \
       foodie_integrator_low_storage_runge_kutta.o \
       foodie_integrator_tvd_runge_kutta.o \
       Data_Type_Command_Line_Interface.o \
       IR_Precision.o \
       pyplot_module.o

.f90.o:
	$(FC) -c $(FCFLAGS) $<

.F90.o:
	$(FC) -c $(FCFLAGS) $<

all:
	ln -sfv ../../external/IR_Precision/src/IR_Precision.f90 .
	ln -sfv ../../external/pyplot-fortran/src/pyplot_module.f90 .
	ln -sfv ../../external/FLAP/src/lib/Data_Type_Command_Line_Interface.F90 .
	$(MAKE) foodie
	ar ruv libfoodie.a $(OBJS)

foodie: foodie.o

foodie.o: foodie.f90 $(OBJS)
foodie_integrator_adams_bashforth.o: foodie_integrator_adams_bashforth.f90 type_integrand.o IR_Precision.o
foodie_integrator_euler_explicit.o: foodie_integrator_euler_explicit.f90 type_integrand.o IR_Precision.o
foodie_integrator_leapfrog.o: foodie_integrator_leapfrog.f90 type_integrand.o IR_Precision.o
foodie_integrator_low_storage_runge_kutta.o: foodie_integrator_low_storage_runge_kutta.f90 type_integrand.o IR_Precision.o
foodie_integrator_tvd_runge_kutta.o: foodie_integrator_tvd_runge_kutta.f90 type_integrand.o IR_Precision.o
type_integrand.o: type_integrand.f90 IR_Precision.o
IR_precision.o: IR_Precision.f90
pyplot_module.o: pyplot_module.f90
Data_Type_Command_Line_Interface.o: Data_Type_Command_Line_Interface.F90

clean:
	rm -f *.o *.mod *.a
