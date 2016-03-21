#
# FOODIE Makefile
#
########################################################################
# Compiler and flags

### IBM XL ###
FC      = bgxlf2008_r
FCFLAGS = -qsuffix=cpp=f90

### GNU ###
#FC      = gfortran
#FCFLAGS = -cpp -g -O0 -C -fbacktrace

### Intel ###
#FC      = ifort
#FCFLAGS = -cpp -g -O0 -C -traceback -assume realloc_lhs

### Cray ###
#FC      = ftn
#FCFLAGS = -O 0 -e Z -g

FCFLAGS_PASS = "$(FCFLAGS)"
OPTSC = $(FCFLAGS_PASS)" -c"

########################################################################
# Paths

PWD       = $(shell pwd)
LIBDIR    = $(PWD)/src/lib
FLAPDIR   = $(PWD)/src/third_party/FLAP
PYPLOTDIR = $(PWD)/src/third_party/pyplot-fortran
WENOOFDIR = $(PWD)/src/third_party/WenOOF

########################################################################
# Targets

.PHONY: all external flap pyplot wenoof foodie tests lorenz oscillation burgers clean

all: tests

foodie:
	@echo "Building $@"
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS_PASS) --directory=src/lib

########################################################################
# External (third-party) libraries

external: flap wenoof pyplot

flap:
	@echo "Building $@ (third-party library)"
	$(MAKE) FC=$(FC) OPTSC=$(OPTSC) --directory=$(FLAPDIR)
	ar ruv $(FLAPDIR)/libflap.a $(FLAPDIR)/tests/obj/data_type_command_line_interface.o $(FLAPDIR)/tests/obj/ir_precision.o

pyplot:
	@echo "Building $@ (third-party library)"
	cd $(PYPLOTDIR)/src && $(FC) -c $(FCFLAGS) pyplot_module.f90
	ar ruv $(PYPLOTDIR)/libpyplot.a $(PYPLOTDIR)/src/pyplot_module.o

wenoof:
	@echo "Building $@ (third-party library)"
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS_PASS) --directory=$(WENOOFDIR)/src/lib
	ar ruv $(WENOOFDIR)/libwenoof.a $(WENOOFDIR)/src/lib/*.o

########################################################################

tests: accuracy parallel regression

accuracy: oscillation
parallel: euler-1D-openmp euler-1D-openmp-no-foodie
regression: lorenz burgers euler-1D

oscillation: foodie external
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS_PASS) INCLUDE=$(LIBDIR) LIB=$(LIBDIR) \
                FLAPDIR=$(FLAPDIR) PYPLOTDIR=$(PYPLOTDIR) --directory=src/tests/accuracy/$@
	cp src/tests/accuracy/$@/$@ .

lorenz: foodie external
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS_PASS) INCLUDE=$(LIBDIR) LIB=$(LIBDIR) \
                FLAPDIR=$(FLAPDIR) PYPLOTDIR=$(PYPLOTDIR) --directory=src/tests/regression/$@
	cp src/tests/regression/$@/$@ .

burgers: foodie external
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS_PASS) INCLUDE=$(LIBDIR) LIB=$(LIBDIR) \
                FLAPDIR=$(FLAPDIR) PYPLOTDIR=$(PYPLOTDIR) --directory=src/tests/regression/$@
	cp src/tests/regression/$@/$@ .

euler-1D: foodie external
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS_PASS) INCLUDE=$(LIBDIR) LIB=$(LIBDIR) \
                FLAPDIR=$(FLAPDIR) PYPLOTDIR=$(PYPLOTDIR) WENOOFDIR=$(WENOOFDIR) --directory=src/tests/regression/$@
	cp src/tests/regression/$@/$@ .

euler-1D-openmp:
	@echo "euler-1D-openmp build not implemented in Makefile yet; use FoBIS instead"

euler-1D-openmp-no-foodie:
	@echo "euler-1D-openmp-no-foodie Build not implemented in Makefile yet; use FoBIS instead"

clean:
	rm -vf lorenz oscillation burgers euler-1D
	rm -vf $(FLAPDIR)/libflap.a
	rm -vf $(WENOOFDIR)/libwenoof.a
	rm -vf $(PYPLOTDIR)/libpyplot.a
	$(MAKE) --directory=$(FLAPDIR) clean
	$(MAKE) --directory=$(WENOOFDIR)/src/lib clean
	$(MAKE) --directory=$(LIBDIR) clean
	$(MAKE) --directory=src/tests/accuracy/oscillation clean
	$(MAKE) --directory=src/tests/regression/burgers clean
	$(MAKE) --directory=src/tests/regression/lorenz clean
	$(MAKE) --directory=src/tests/regression/euler-1D clean
