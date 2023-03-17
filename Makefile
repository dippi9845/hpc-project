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

EXE:=sph sph-omp
CFLAGS+=-std=c99 -Wall -Wpedantic
LDLIBS=-lm
OPM_FLAG=-fopenmp
SIMD_FLAG=-march=native -O2 -ftree-vectorize -fopt-info-vec -funsafe-math-optimizations

.PHONY: clean

all: $(EXE)

sph: sph.c
	gcc sph.c $(CFLAGS) $(LDLIBS) -o bin/sph

omp: sph-omp.c
	gcc sph-omp.c $(CFLAGS) $(LDLIBS) $(OPM_FLAG) -o bin/sph-omp

simd: sph-simd.c
	gcc sph-simd.c $(CFLAGS) $(LDLIBS) $(SIMD_FLAG) -o bin/sph-simd

omp-simd: sph-omp-simd.c
	gcc sph-omp-simd.c $(CFLAGS) $(LDLIBS) $(SIMD_FLAG) $(OPM_FLAG) -o bin/sph-omp-simd

simd-acc: sph-simd-acc.c
	gcc sph-simd-acc.c $(CFLAGS) $(LDLIBS) $(SIMD_FLAG) $(OPM_FLAG) -o bin/sph-simd-acc

omp-simd-acc: sph-omp-simd-acc.c
	gcc sph-omp-simd-acc.c $(CFLAGS) $(LDLIBS) $(SIMD_FLAG) $(OPM_FLAG) -o bin/sph-omp-simd-acc

clean:
	\rm -f $(EXE) *.o *~
