#!/bin/bash

EXPORT_PATH=/media/dippi/Volume1/hpc_tests/
MAX_THREAD=24
REPETITIONS=8
EXE_PATH=bin
MAX_STEPS=1000
MAX_PARTICLES=20000

CURRENT_DIR="${EXPORT_PATH}sph"
CURRENT_EXE="${EXE_PATH}/sph"

make sph 1>/dev/null 2>&1 
mkdir $CURRENT_DIR

for (( CUR_PAR=2000; $CUR_PAR<=$MAX_PARTICLES; CUR_PAR=$CUR_PAR+2000 )); do # particles [2000, 20000] -> 10
    for (( CUR_STEP=200; $CUR_STEP<=$MAX_STEPS; CUR_STEP=$CUR_STEP+200 )); do # steps [200, 1000] -> 5
        for (( try=0; $try<$REPETITIONS; try=$try+1 )); do #  * 8
            #$CURRENT_EXE $CUR_PAR 
            echo "[par: $CUR_PAR; steps: $CUR_STEP; R $try;]"
        done
    done
done

for (( i=1; $i<=$MAX_THREAD; i=$i+1 )); do
    for (( try=1; $try<=$REPETITIONS; try=$try+1 )); do
        echo "[runno con $i ths la $try volta ]"
    done
    echo
done