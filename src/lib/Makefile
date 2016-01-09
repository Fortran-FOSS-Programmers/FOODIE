# FOODIE Makefile

.SUFFIXES: .F90 .f90 .o
.PHONY: all foodie clean

OBJS = foodie_adt_integrand.o \
       foodie_integrator_adams_bashforth.o \
       foodie_integrator_adams_bashforth_moulton.o \
       foodie_integrator_adams_moulton.o \
       foodie_integrator_backward_differentiation_formula.o \
       foodie_integrator_euler_explicit.o \
       foodie_integrator_leapfrog.o \
       foodie_integrator_low_storage_runge_kutta.o \
       foodie_integrator_embedded_runge_kutta.o \
       foodie_integrator_tvd_runge_kutta.o \
       foodie_kinds.o \
       foodie_utils.o

.f90.o:
	$(FC) -c $(FCFLAGS) $<

.F90.o:
	$(FC) -c $(FCFLAGS) $<

all:
	$(MAKE) foodie
	ar ruv libfoodie.a $(OBJS)

foodie: foodie.o

foodie.o: foodie.f90 $(OBJS)
foodie_integrator_adams_bashforth.o: foodie_integrator_adams_bashforth.f90 foodie_adt_integrand.o foodie_kinds.o foodie_utils.o
foodie_integrator_adams_bashforth_moulton.o: foodie_integrator_adams_bashforth_moulton.f90 foodie_integrator_adams_bashforth.o foodie_integrator_adams_moulton.o foodie_adt_integrand.o foodie_kinds.o foodie_utils.o
foodie_integrator_adams_moulton.o: foodie_integrator_adams_moulton.f90 foodie_adt_integrand.o foodie_kinds.o foodie_utils.o
foodie_integrator_backward_differentiation_formula.o: foodie_integrator_backward_differentiation_formula.f90 foodie_adt_integrand.o foodie_kinds.o
foodie_integrator_euler_explicit.o: foodie_integrator_euler_explicit.f90 foodie_adt_integrand.o foodie_kinds.o
foodie_integrator_leapfrog.o: foodie_integrator_leapfrog.f90 foodie_adt_integrand.o foodie_kinds.o
foodie_integrator_embedded_runge_kutta.o: foodie_integrator_embedded_runge_kutta.f90 foodie_adt_integrand.o foodie_kinds.o foodie_utils.o
foodie_integrator_low_storage_runge_kutta.o: foodie_integrator_low_storage_runge_kutta.f90 foodie_adt_integrand.o foodie_kinds.o foodie_utils.o
foodie_integrator_tvd_runge_kutta.o: foodie_integrator_tvd_runge_kutta.f90 foodie_adt_integrand.o foodie_kinds.o
foodie_adt_integrand.o: foodie_adt_integrand.f90 foodie_kinds.o
foodie_kinds.o: foodie_kinds.f90
foodie_utils.o: foodie_utils.f90

clean:
	rm -f *.o *.mod *.a
