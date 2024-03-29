#!/bin/bash

test_loop() {
    EXE_NAME=$1
    MAKE_NAME=$2

    CURRENT_DIR="${EXPORT_PATH}${EXE_NAME}"
    CURRENT_EXE="${EXE_PATH}/${EXE_NAME}"

    make $MAKE_NAME
    mkdir $CURRENT_DIR

    for (( CUR_PAR=500; $CUR_PAR<=$MAX_PARTICLES; CUR_PAR=$CUR_PAR+2000 )); do # particles [1000, 20000] -> 20
        echo "[P: $CUR_PAR]"
        for (( CUR_STEP=200; $CUR_STEP<=200; CUR_STEP=$CUR_STEP+50 )); do # steps [100, 100] -> 1
            TO_PRINT=""
            for (( try=0; $try<$REPETITIONS; try=$try+1 )); do #  * 8 
                OUT=`$CURRENT_EXE $CUR_PAR $CUR_STEP`
                TO_PRINT="$TO_PRINT;$OUT"
            done
            echo ${TO_PRINT:1} > "$CURRENT_DIR/run_p_${CUR_PAR}_s_${CUR_STEP}.csv"
        done
    done
}

test_loop_parallel() {
    EXE_NAME=$1
    MAKE_NAME=$2
    TH=$3

    CURRENT_DIR="${EXPORT_PATH}${EXE_NAME}"
    CURRENT_EXE="${EXE_PATH}/${EXE_NAME}"

    make $MAKE_NAME 1>/dev/null 2>&1 
    mkdir $CURRENT_DIR 2>/dev/null

    for (( CUR_PAR=6000; $CUR_PAR<=$MAX_PARTICLES; CUR_PAR=$CUR_PAR+500 )); do # particles [500, 6000] -> 10
        echo "        [P: $CUR_PAR]"
        for (( CUR_STEP=100; $CUR_STEP<=$MAX_STEPS; CUR_STEP=$CUR_STEP+50 )); do # steps [50, 200] -> 4
            TO_PRINT=""
            echo "              [S: $CUR_STEP]"
            for (( try=0; $try<$REPETITIONS; try=$try+1 )); do #  * 8 
                OUT=`OMP_NUM_THREADS=${3} $CURRENT_EXE $CUR_PAR $CUR_STEP`
                #echo "$CURRENT_EXE $CUR_PAR $CUR_STEP"
                TO_PRINT="$TO_PRINT;$OUT"
            done
            echo ${TO_PRINT:1} > "$CURRENT_DIR/run_p_${CUR_PAR}_s_${CUR_STEP}_t_${TH}.csv"
        done
    done
}

EXPORT_PATH=/media/dippi/Volume1/hpc_tests/
MAX_THREAD=24
REPETITIONS=8
EXE_PATH=../bin
MAX_STEPS=100
MAX_PARTICLES=6000

#echo "Inizio Versione seriale"

#test_loop "sph" "sph"

#echo "Finito la versione seriale"

#test_loop "sph-simd" "simd"

#test_loop "sph-simd-acc" "simd-acc"

#echo "Finito la versione simd"

#echo "inizio omp"

#for (( i=1; $i<=$MAX_THREAD; i=$i+1 )); do
#    echo "    [th: $i]" 
#    test_loop_parallel "sph-omp-dynamic" "omp-dynamic" $i
#done

#echo "fine omp"
#echo "inizio omp simd"

#for (( i=1; $i<=$MAX_THREAD; i=$i+1 )); do
#    echo "    [th: $i]" 
#    test_loop_parallel "sph-omp-simd" "omp-simd" $i
#done

#echo "Dimanica"

#for (( i=1; $i<=$MAX_THREAD; i=$i+1 )); do
#    echo "    [th: $i]" 
#    test_loop_parallel "sph-omp-dynamic" "omp-dynamic" $i
#done

#for (( i=1; $i<=$MAX_THREAD; i=$i+1 )); do
#    echo "    [th: $i]" 
#    test_loop_parallel "sph-omp-simd" "omp-simd" $i
#done

#echo "Dimanica simd"

#for (( i=1; $i<=$MAX_THREAD; i=$i+1 )); do
#    echo "    [th: $i]" 
#    test_loop_parallel "sph-omp-simd" "omp-simd" $i
#done

#echo "fine omp simd"

#echo "Inizio cuda"

#test_loop "sph-cuda" "cuda"

#echo "Fine cuda"

#echo "Inizio cuda"

#test_loop "sph-cuda-SoA" "cuda-SoA"

#echo "Fine cuda"

#echo "Inizio cuda"

#test_loop "sph-cuda-shared" "cuda-shared"

#echo "Fine cuda"
