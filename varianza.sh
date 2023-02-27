#!/bin/bash
STEPS=100
CORES=1
for (( i=1000; $i<20000; i=$i*2 )); do
    OMP_NUM_THREADS=${CORES} $1 $i ${STEPS} > "varianza_${i}_${CORES}.txt"
done

CORES=12
for (( i=1000; $i<20000; i=$i*2 )); do
    OMP_NUM_THREADS=${CORES} $1 $i ${STEPS} > "varianza_${i}_${CORES}.txt"
done