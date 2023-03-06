#!/bin/bash
# comando per pescare i campi: cut -d";" -f FIELD_NUM FILE

serial_union() {

    SERIAL="${EXPORT_PATH}sph"
    SERIAL_SIMD="${EXPORT_PATH}sph-simd"
    SERIAL_OUTPUT="${EXPORT_PATH}serial-simd.csv"

    RAW=""
    for (( CUR_PAR=500; $CUR_PAR<=$MAX_PARTICLES; CUR_PAR=$CUR_PAR+500 )); do # particles [2000, 20000] -> 10
        for (( CUR_STEP=200; $CUR_STEP<=$MAX_STEPS; CUR_STEP=$CUR_STEP+50 )); do # steps [200, 200] -> 1
            SERIAL_FILE="${SERIAL}/run_p_${CUR_PAR}_s_${CUR_STEP}.csv"
            SERIAL_SIMD_FILE="${SERIAL_SIMD}/run_p_${CUR_PAR}_s_${CUR_STEP}.csv"
            if [ -e $SERIAL_FILE -a -e $SERIAL_SIMD_FILE ]; then
                RAW="${RAW};Seriale Particelle: ${CUR_PAR};Seriale SIMD Particelle: ${CUR_PAR}"
            fi
        done
    done
    RAW="${RAW:1}\n" # salto il primo ;

    for (( try=1; $try<=$REPETITIONS; try=$try+1 )); do
        ROW=""
        for (( CUR_PAR=500; $CUR_PAR<=$MAX_PARTICLES; CUR_PAR=$CUR_PAR+500 )); do # particles [2000, 20000] -> 10
            for (( CUR_STEP=200; $CUR_STEP<=$MAX_STEPS; CUR_STEP=$CUR_STEP+50 )); do # steps [200, 200] -> 1
                SERIAL_FILE="${SERIAL}/run_p_${CUR_PAR}_s_${CUR_STEP}.csv"
                SERIAL_SIMD_FILE="${SERIAL_SIMD}/run_p_${CUR_PAR}_s_${CUR_STEP}.csv"
                if [ -e $SERIAL_FILE -a -e $SERIAL_SIMD_FILE ]; then
                    timeing1=`cut -d";" -f $try $SERIAL_FILE`
                    timeing2=`cut -d";" -f $try $SERIAL_SIMD_FILE`
                    ROW="${ROW};${timeing1};${timeing2}"
                fi
            done
        done
        RAW="${RAW}${ROW:1}\n"
    done

    echo -n -e $RAW > $SERIAL_OUTPUT

}

omp_union() {
    OMP="${EXPORT_PATH}sph-omp${1}"
    OMP_OUTPUT="${EXPORT_PATH}omp${1}.csv"

    RAW=""
    for (( CUR_PAR=500; $CUR_PAR<=$MAX_PARTICLES; CUR_PAR=$CUR_PAR+500 )); do # particles [2000, 20000] -> 10
        for (( CUR_STEP=50; $CUR_STEP<=$MAX_STEPS; CUR_STEP=$CUR_STEP+50 )); do # steps [200, 200] -> 1
            for (( THD=1; $THD<=$MAX_THREAD; THD=$THD+1 )); do # per tutti i thread
                OMP_FILE="${OMP}/run_p_${CUR_PAR}_s_${CUR_STEP}_t_${THD}.csv"
                if [ -e $OMP_FILE ]; then
                    RAW="${RAW};Particelle: ${CUR_PAR} Step: ${CUR_STEP} TH ${THD}"
                fi
            done
        done
    done
    RAW="${RAW:1}\n" # salto il primo ;
    
    for (( try=1; $try<=$REPETITIONS; try=$try+1 )); do
        ROW=""
        for (( CUR_PAR=500; $CUR_PAR<=$MAX_PARTICLES; CUR_PAR=$CUR_PAR+500 )); do # particles [2000, 20000] -> 10
            for (( CUR_STEP=50; $CUR_STEP<=$MAX_STEPS; CUR_STEP=$CUR_STEP+50 )); do # steps [200, 200] -> 1
                for (( THD=1; $THD<=$MAX_THREAD; THD=$THD+1 )); do # per tutti i thread
                    OMP_FILE="${OMP}/run_p_${CUR_PAR}_s_${CUR_STEP}_t_${THD}.csv"
                    if [ -e $OMP_FILE ]; then
                        timeing=`cut -d";" -f $try $OMP_FILE`
                        ROW="${ROW};${timeing}"
                    fi
                done
            done
        done
        RAW="${RAW}${ROW:1}\n"
    done

    #echo -n -e $RAW

    echo -n -e $RAW > $OMP_OUTPUT
}

EXPORT_PATH=/media/dippi/Volume1/hpc_tests/
MAX_THREAD=24
REPETITIONS=8
EXE_PATH=bin
MAX_STEPS=200
MAX_PARTICLES=6000

# unisci i tempi simd e seriali
#serial_union

# unisci i tempi paralleli di omp
#omp_union

# unisci i tempi paralleli di omp-simd
#omp_union "-simd"
