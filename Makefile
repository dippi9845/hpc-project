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
STEP_PERF=-D"STEP_PERFORMANCE"

.PHONY: clean

all: sph omp cuda

sph: sph.c
	gcc sph.c $(CFLAGS) $(LDLIBS) -o bin/sph

sph-step: sph.c
	gcc sph.c $(STEP_PERF) $(CFLAGS) $(LDLIBS) -o bin/sph-step

omp: sph-omp.c
	gcc sph-omp.c $(CFLAGS) $(LDLIBS) $(OPM_FLAG) -o bin/sph-omp

omp-step: sph-omp.c
	gcc sph-omp.c $(STEP_PERF) $(CFLAGS) $(LDLIBS) $(OPM_FLAG) -o bin/sph-omp-step

simd: sph-simd.c
	gcc sph-simd.c $(CFLAGS) $(LDLIBS) $(SIMD_FLAG) -o bin/sph-simd

simd-step: sph-simd.c
	gcc sph-simd.c $(STEP_PERF) $(CFLAGS) $(LDLIBS) $(SIMD_FLAG) -o bin/sph-simd-step

omp-simd: sph-omp-simd.c
	gcc sph-omp-simd.c $(CFLAGS) $(LDLIBS) $(SIMD_FLAG) $(OPM_FLAG) -o bin/sph-omp-simd

omp-simd-step: sph-omp-simd.c
	gcc sph-omp-simd.c $(STEP_PERF) $(CFLAGS) $(LDLIBS) $(SIMD_FLAG) $(OPM_FLAG) -o bin/sph-omp-simd-step

cuda: sph-cuda.cu
	nvcc sph-cuda.cu $(LDLIBS) -o bin/sph-cuda

cuda-step: sph-cuda.cu
	nvcc sph-cuda.cu $(STEP_PERF) $(LDLIBS) -o bin/sph-cuda-step

cuda-SoA: sph-cuda-SoA.cu
	nvcc sph-cuda-SoA.cu $(LDLIBS) -o bin/sph-cuda-SoA

cuda-shared: sph-cuda-shared.cu
	nvcc sph-cuda-shared.cu $(LDLIBS) -o bin/sph-cuda-shared

cuda-dbg: sph-cuda.cu
	nvcc -g -G sph-cuda.cu $(LDLIBS) -o bin/sph-cuda

simd-acc: sph-simd-acc.c
	gcc sph-simd-acc.c $(CFLAGS) $(LDLIBS) $(SIMD_FLAG) $(OPM_FLAG) -o bin/sph-simd-acc

clean:
	\rm -f bin/sph*
