# FOODiE Makefile

.SUFFIXES: .f90 .o
.PHONY: all foodie clean

.f90.o:
	$(FC) -c $(FCFLAGS) $<

all:
	ln -sfv ../tests/pyplot_module.f90 .
	$(MAKE) foodie

foodie: foodie.o

foodie.o: foodie.f90 type_integrand.o foodie_integrator_euler_explicit.o foodie_integrator_tvd_runge_kutta.o pyplot_module.o 
foodie_integrator_euler_explicit.o: foodie_integrator_euler_explicit.f90 type_integrand.o IR_Precision.o
#foodie_integrator_euler_implicit.o: foodie_integrator_euler_implicit.f90 type_integrand.o IR_Precision.o
foodie_integrator_tvd_runge_kutta.o: foodie_integrator_tvd_runge_kutta.f90 type_integrand.o IR_Precision.o
type_integrand.o: type_integrand.f90 IR_Precision.o
IR_precision.o: IR_Precision.f90
pyplot_module.o: pyplot_module.f90

clean:
	rm -f *.o *.mod
