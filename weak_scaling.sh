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
    for (( try=0; $try<$REPETITIONS; try=$try+1 )); do #  8 
        OUT=`OMP_NUM_THREADS=${3} $CURRENT_EXE $N $S`
        echo -e "\t[NUM: ${try}] $OUT"
        TO_PRINT="$TO_PRINT;$OUT"
    done
    echo ${TO_PRINT:1} > "$CURRENT_DIR/run_p_${N}_s_${S}_t_${TH}.csv"
}

omp_union_particles() {
    OMP="${EXPORT_PATH}sph-omp${1}"
    OMP_OUTPUT="${EXPORT_PATH}omp-weak${1}.csv"
    CUR_STEP=$STEPS

    RAW=""
    for (( THD=1; $THD<=$MAX_THREAD; THD=$THD+1 )); do
        CUR_PAR=`printf "%.0f\n" $(echo "scale=15; ${C_P}*sqrt($THD)" | bc) 2>/dev/null`
        OMP_FILE="${OMP}/run_p_${CUR_PAR}_s_${CUR_STEP}_t_${THD}.csv"
        if [ -e $OMP_FILE ]; then
            RAW="${RAW};${THD}"
        fi
    done
    RAW="${RAW:1}\n"
    
    for (( try=1; $try<=$REPETITIONS; try=$try+1 )); do
        ROW=""
        for (( THD=1; $THD<=$MAX_THREAD; THD=$THD+1 )); do
            CUR_PAR=`printf "%.0f\n" $(echo "scale=15; ${C_P}*sqrt($THD)" | bc) 2>/dev/null`
            OMP_FILE="${OMP}/run_p_${CUR_PAR}_s_${CUR_STEP}_t_${THD}.csv"
            if [ -e $OMP_FILE ]; then
                timeing=`cut -d";" -f $try $OMP_FILE`
                ROW="${ROW};${timeing}"
            fi
        done
        RAW="${RAW}${ROW:1}\n"
    done

    echo -n -e $RAW > $OMP_OUTPUT
}

omp_union_steps() {
    OMP="${EXPORT_PATH}sph-omp${1}"
    OMP_OUTPUT="${EXPORT_PATH}omp-weak${1}.csv"
    CUR_PAR=$PARTICLES

    RAW=""
    for (( THD=1; $THD<=$MAX_THREAD; THD=$THD+1 )); do
        CUR_STEP=`printf "%.0f\n" $(echo "scale=15; ${C_S}*$THD" | bc) 2>/dev/null`
        OMP_FILE="${OMP}/run_p_${CUR_PAR}_s_${CUR_STEP}_t_${THD}.csv"
        echo "[S: $CUR_STEP][$OMP_FILE]"
        if [ -e $OMP_FILE ]; then
            RAW="${RAW};${THD}"
        fi
    done
    RAW="${RAW:1}\n"
    
    for (( try=1; $try<=$REPETITIONS; try=$try+1 )); do
        ROW=""
        for (( THD=1; $THD<=$MAX_THREAD; THD=$THD+1 )); do
            CUR_STEP=`printf "%.0f\n" $(echo "scale=15; ${C_S}*$THD" | bc) 2>/dev/null`
            OMP_FILE="${OMP}/run_p_${CUR_PAR}_s_${CUR_STEP}_t_${THD}.csv"
            if [ -e $OMP_FILE ]; then
                timeing=`cut -d";" -f $try $OMP_FILE`
                ROW="${ROW};${timeing}"
            fi
        done
        RAW="${RAW}${ROW:1}\n"
    done

    echo -n -e $RAW > $OMP_OUTPUT
}

REPETITIONS=8
MAX_THREAD=24
MAKE_NAME="omp"
EXE_NAME="sph-omp"
STEPS=100
PARTICLES=1500
EXE_PATH=bin
EXPORT_PATH=/media/dippi/Volume1/hpc_tests/
C_P="1500"
C_S="100"

for (( i=1 ; $i<=$MAX_THREAD ; i=$i+1 )); do
    echo "[TH: ${i}]"
    N=`printf "%.0f\n" $(echo "scale=15; ${C_P}*sqrt($i)" | bc) 2>/dev/null`
    echo "[PARTS: $N]"
    #run_parallel $EXE_NAME $MAKE_NAME $i $N $STEPS
done

#omp_union_particles

for (( i=1 ; $i<=$MAX_THREAD ; i=$i+1 )); do
    echo "[TH: ${i}]"
    S=`printf "%.0f\n" $(echo "scale=15; ${C_S}*$i" | bc) 2>/dev/null`
    echo "[STEPS: $S]"
    #run_parallel $EXE_NAME $MAKE_NAME $i $PARTICLES $S
done

omp_union_steps