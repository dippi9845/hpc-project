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
SRC_FOLDER=src/
BIN_FOLDER=bin/

.PHONY: clean

all: sph simd omp-dynamic omp-one-thread-pool omp-simd-SoA cuda cuda-SoA cuda-shared

consegna: omp-simd-SoA cuda-shared

sph: ${SRC_FOLDER}sph.c
	gcc ${SRC_FOLDER}sph.c $(CFLAGS) $(LDLIBS) -o ${BIN_FOLDER}sph

simd: ${SRC_FOLDER}sph-simd.c
	gcc ${SRC_FOLDER}sph-simd.c $(CFLAGS) $(LDLIBS) $(SIMD_FLAG) -o ${BIN_FOLDER}sph-simd

omp-dynamic: ${SRC_FOLDER}sph-omp-dynamic.c
	gcc ${SRC_FOLDER}sph-omp-dynamic.c $(CFLAGS) $(LDLIBS) $(OPM_FLAG) -o ${BIN_FOLDER}sph-omp-dynamic

omp-simd-SoA: ${SRC_FOLDER}sph-omp-simd-SoA.c
	gcc ${SRC_FOLDER}sph-omp-simd-SoA.c $(CFLAGS) $(LDLIBS) $(SIMD_FLAG) $(OPM_FLAG) -o ${BIN_FOLDER}sph-omp-simd-SoA

omp-one-thread-pool: ${SRC_FOLDER}sph-omp-one-thread-pool.c
	gcc ${SRC_FOLDER}sph-omp-one-thread-pool.c $(CFLAGS) $(LDLIBS) $(OPM_FLAG) -o ${BIN_FOLDER}sph-omp-one-thread-pool

cuda: ${SRC_FOLDER}sph-cuda.cu
	nvcc ${SRC_FOLDER}sph-cuda.cu $(LDLIBS) -o ${BIN_FOLDER}sph-cuda

cuda-SoA: ${SRC_FOLDER}sph-cuda-SoA.cu
	nvcc ${SRC_FOLDER}sph-cuda-SoA.cu $(LDLIBS) -o ${BIN_FOLDER}sph-cuda-SoA

cuda-shared: ${SRC_FOLDER}sph-cuda-shared.cu
	nvcc ${SRC_FOLDER}sph-cuda-shared.cu $(LDLIBS) -o ${BIN_FOLDER}sph-cuda-shared

clean:
	\rm -f bin/sph*
