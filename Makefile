# Makefile for Time-Dependent Schrödinger Equation Solver
FC ?= gfortran
FFLAGS = -O3 -Wall -fdefault-real-8
OMPFLAGS = -fopenmp
PYTHON ?= python3

.PHONY: all clean run plot parallel

all: schr armonico

schr: schr.f95
	$(FC) $(FFLAGS) -o $@ schr.f95

armonico: armonico.f95
	$(FC) $(FFLAGS) -o $@ armonico.f95

parallel: multi/multOpen.f95
	$(FC) $(FFLAGS) $(OMPFLAGS) -o multi/multOpen multi/multOpen.f95

run: schr
	./schr

plot:
	$(PYTHON) representar.py

clean:
	rm -f schr armonico mult schr2 multi/multi multi/multOpen *.o *.mod *.exe fort.* onda_resultado.dat

