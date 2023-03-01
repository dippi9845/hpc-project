#!/bin/bash

test_loop() {
    EXE_NAME=$1
    MAKE_NAME=$2

    CURRENT_DIR="${EXPORT_PATH}${EXE_NAME}"
    CURRENT_EXE="${EXE_PATH}/${EXE_NAME}"

    make $MAKE_NAME 1>/dev/null 2>&1 
    mkdir $CURRENT_DIR 2>/dev/null

    for (( CUR_PAR=2000; $CUR_PAR<=$MAX_PARTICLES; CUR_PAR=$CUR_PAR+2000 )); do # particles [2000, 20000] -> 10
        echo "[P: $CUR_PAR]"
        for (( CUR_STEP=200; $CUR_STEP<=$MAX_STEPS; CUR_STEP=$CUR_STEP+200 )); do # steps [200, 1000] -> 5
            TO_PRINT=""
            for (( try=0; $try<$REPETITIONS; try=$try+1 )); do #  * 8 
                OUT=`$CURRENT_EXE $CUR_PAR $CUR_STEP`
                TO_PRINT="$TO_PRINT;$OUT"
            done
            echo ${TO_PRINT:1} > "$CURRENT_DIR/run_p_${CUR_PAR}_s_${CUR_STEP}${3}.csv"
        done
    done
}

EXPORT_PATH=/media/dippi/Volume1/hpc_tests/
MAX_THREAD=24
REPETITIONS=8
EXE_PATH=bin
MAX_STEPS=1000
MAX_PARTICLES=20000

#test_loop "sph" "sph"

#echo "Finito la versione seriale"

#test_loop "sph-simd" "simd"

#echo "Finito la versione simd"

echo "inizio omp"

for (( i=1; $i<=$MAX_THREAD; i=$i+1 )); do
    test_loop "sph-omp" "omp" "t_${i}"
done

echo "fine omp"
echo "inizio omp simd"

for (( i=1; $i<=$MAX_THREAD; i=$i+1 )); do
    test_loop "sph-omp" "omp" "t_${i}"
done

echo "fine omp simd"