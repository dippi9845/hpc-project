#!/bin/bash

run_parallel() {
    EXE_NAME=$1
    MAKE_NAME=$2
    TH=$3
    N=$4
    S=$5

    CURRENT_DIR="${EXPORT_PATH}${EXE_NAME}"
    CURRENT_EXE="${EXE_PATH}/${EXE_NAME}"

    make $MAKE_NAME 1>/dev/null 2>&1 
    mkdir $CURRENT_DIR 2>/dev/null

    TO_PRINT=""
    for (( try=0; $try<$REPETITIONS; try=$try+1 )); do #  * 8 
        OUT=`OMP_NUM_THREADS=${3} $CURRENT_EXE $N $S`
        TO_PRINT="$TO_PRINT;$OUT"
    done
    echo ${TO_PRINT:1} > "$CURRENT_DIR/run_p_${N}_s_${S}_t_${TH}.csv"
}

$MAX_THREAD=24
MAKE_NAME=""
EXE_NAME="" 
STEPS=100

for (( i=0 ; $i<$MAX_THREAD ; i=$i+1 )); do
    N=printf "%.0f\n" $(echo "1000*sqrt(4.713423831*${i})" | bc) 2>/dev/null
    run_parallel $MAKE_NAME $EXE_NAME $i $N $STEPS
done