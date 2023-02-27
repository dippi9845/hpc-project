#!/bin/bash
for (( i=1000; $i<20000; i=$i*2 )); do
    OMP_NUM_THREADS=1 $1 $i 100 > "varianza_${i}_100_1.txt"
done