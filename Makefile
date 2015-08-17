# FOODiE Makefile

### GNU ###
FC      = gfortran
FCFLAGS = "-cpp -g -O0 -C -fbacktrace"

### Intel ###
#FC      = ifort
#FCFLAGS = "-cpp -g -O0 -C -traceback -assume realloc_lhs"

### Cray ###
#FC      = ftn
#FCFLAGS = "-O 0 -e Z -g"


.PHONY: all foodie tests lorenz oscillation burgers clean

all: tests

foodie:
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS) --directory=src/lib

tests: lorenz oscillation burgers

lorenz: foodie
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS) --directory=src/tests/$@
	cp src/tests/$@/$@ .

oscillation: foodie
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS) --directory=src/tests/$@
	cp src/tests/$@/$@ .

burgers: foodie
	$(MAKE) FC=$(FC) FCFLAGS=$(FCFLAGS) --directory=src/tests/$@
	cp src/tests/$@/$@ .

clean:
	rm -vf lorenz oscillation burgers
	$(MAKE) --directory=src/lib clean
	$(MAKE) --directory=src/tests/lorenz clean
	$(MAKE) --directory=src/tests/burgers clean
