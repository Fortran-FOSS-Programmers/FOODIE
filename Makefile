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
# Targets

.PHONY: all external foodie tests lorenz oscillation burgers clean

all: tests

foodie:
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS) --directory=src/lib

external:
	$(MAKE) FC=$(FC) OPTSC=$(OPTSC) --directory=external/FLAP
	ar ruv external/FLAP/libexternal.a external/FLAP/tests/obj/data_type_command_line_interface.o external/FLAP/tests/obj/ir_precision.o

tests: lorenz oscillation burgers

lorenz: foodie external
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS) --directory=src/tests/$@
	cp src/tests/$@/$@ .

oscillation: foodie external
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS) --directory=src/tests/$@
	cp src/tests/$@/$@ .

burgers: foodie external
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS) --directory=src/tests/$@
	cp src/tests/$@/$@ .

clean:
	rm -vf lorenz oscillation burgers
	rm -vf external/FLAP/libexternal.a
	$(MAKE) --directory=external/FLAP clean
	$(MAKE) --directory=src/lib clean
	$(MAKE) --directory=src/tests/lorenz clean
	$(MAKE) --directory=src/tests/oscillation clean
	$(MAKE) --directory=src/tests/burgers clean
