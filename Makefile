# Makefile for the High Performance Computing programming project,
# Academic Year 2022/2023.
#
# Available targets:
#
# - sph: builds the non-GUI version (default)
#
# - sph.gui: builds the GUI version
#
# - all: builds both the GUI and non-GUI versions
#
# - clean: clean up
#
# Last modified on 2022-11-27 by Moreno Marzolla

EXE:=sph sph-omp.c
CFLAGS+=-std=c99 -Wall -Wpedantic
LDLIBS=-lm
OPM_FLAG=-fopenmp

.PHONY: clean

all: $(EXE)

sph: sph.c
	gcc sph.c $(CFLAGS) $(LDLIBS) -o bin/sph

omp: sph-omp.c
	gcc sph-omp.c $(CFLAGS) $(LDLIBS) $(OPM_FLAG) -o bin/sph-omp

clean:
	\rm -f $(EXE) *.o *~
