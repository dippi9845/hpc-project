#!/bin/bash

EXPORT_PATH=/media/dippi/Volume1/hpc_tests/
MAX_THREAD=24
REPETITIONS=8
EXE_PATH=bin/

CURRENT_DIR="${EXPORT_PATH}/sph"
CURRENT_EXE="${EXE_PATH}/sph"

make sph
mkdir CURRENT_DIR

for (( try=1; $try<=$try; i=$try+1 )); do
    echo "[runno i sph: $try]"
done

for (( i=1; $i<=$MAX_THREAD; i=$i+1 )); do
    for (( try=1; $try<=$REPETITIONS; try=$try+1 )); do
        echo "[runno con $i ths la $try volta ]"
    done
done