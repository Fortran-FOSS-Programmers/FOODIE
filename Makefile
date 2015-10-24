#
# FOODIE Makefile
#
########################################################################
# Compiler and flags

### GNU ###
FC      = gfortran
FCFLAGS = "-cpp -g -O0 -C -fbacktrace"

### Intel ###
#FC      = ifort
#FCFLAGS = "-cpp -g -O0 -C -traceback -assume realloc_lhs"

### Cray ###
#FC      = ftn
#FCFLAGS = "-O 0 -e Z -g"

OPTSC = $(FCFLAGS)" -c"

########################################################################
# Paths

LIBDIR  = $(shell pwd)/src/lib
FLAPDIR = $(shell pwd)/external/FLAP

########################################################################
# Targets

.PHONY: all external foodie tests lorenz oscillation burgers clean

all: tests

foodie:
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS) --directory=src/lib

external:
	$(MAKE) FC=$(FC) OPTSC=$(OPTSC) --directory=external/FLAP
	ar ruv external/FLAP/libexternal.a external/FLAP/tests/obj/data_type_command_line_interface.o external/FLAP/tests/obj/ir_precision.o

tests: accuracy parallel regression

accuracy: oscillation
parallel: euler-1D-openmp euler-1D-openmp-no-foodie
regression: lorenz burgers

oscillation: foodie external
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS) INCLUDE=$(LIBDIR) LIB=$(LIBDIR) EXTERNAL=$(FLAPDIR) --directory=src/tests/accuracy/$@
	cp src/tests/accuracy/$@/$@ .

lorenz: foodie external
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS) INCLUDE=$(LIBDIR) LIB=$(LIBDIR) EXTERNAL=$(FLAPDIR) --directory=src/tests/regression/$@
	cp src/tests/regression/$@/$@ .

burgers: foodie external
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS) INCLUDE=$(LIBDIR) LIB=$(LIBDIR) EXTERNAL=$(FLAPDIR) --directory=src/tests/regression/$@
	cp src/tests/regression/$@/$@ .

euler-1D-openmp:
	@echo "euler-1D-openmp build not implemented in Makefile yet; use FoBIS instead"

euler-1D-openmp-no-foodie:
	@echo "euler-1D-openmp-no-foodie Build not implemented in Makefile yet; use FoBIS instead"

clean:
	rm -vf lorenz oscillation burgers
	rm -vf external/FLAP/libexternal.a
	$(MAKE) --directory=external/FLAP clean
	$(MAKE) --directory=src/lib clean
	$(MAKE) --directory=src/tests/accuracy/oscillation clean
	$(MAKE) --directory=src/tests/regression/burgers clean
	$(MAKE) --directory=src/tests/regression/lorenz clean
