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


CFLAGS+=-std=c99 -Wall -Wpedantic
LDLIBS=-lm
OPM_FLAG=-fopenmp
SIMD_FLAG=-march=native -O2 -ftree-vectorize -fopt-info-vec -funsafe-math-optimizations

.PHONY: clean

all: sph omp cuda

sph: sph.c
	gcc sph.c $(CFLAGS) $(LDLIBS) -o bin/sph

omp: sph-omp.c
	gcc sph-omp.c $(CFLAGS) $(LDLIBS) $(OPM_FLAG) -o bin/sph-omp

omp-dynamic: sph-omp-dynamic.c
	gcc sph-omp-dynamic.c $(CFLAGS) $(LDLIBS) $(OPM_FLAG) -o bin/sph-omp-dynamic

simd: sph-simd.c
	gcc sph-simd.c $(CFLAGS) $(LDLIBS) $(SIMD_FLAG) -o bin/sph-simd

omp-simd: sph-omp-simd.c
	gcc sph-omp-simd.c $(CFLAGS) $(LDLIBS) $(SIMD_FLAG) $(OPM_FLAG) -o bin/sph-omp-simd

cuda: sph-cuda.cu
	nvcc sph-cuda.cu $(LDLIBS) -o bin/sph-cuda

cuda-SoA: sph-cuda-SoA.cu
	nvcc sph-cuda-SoA.cu $(LDLIBS) -o bin/sph-cuda-SoA

cuda-shared: sph-cuda-shared.cu
	nvcc sph-cuda-shared.cu $(LDLIBS) -o bin/sph-cuda-shared

cuda-dbg: sph-cuda.cu
	nvcc -g -G sph-cuda.cu $(LDLIBS) -o bin/sph-cuda

clean:
	\rm -f bin/sph*
