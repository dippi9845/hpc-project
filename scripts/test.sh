#!/bin/bash

COMMAND=$1
for (( i=100; i<=2000; i=$i+100 )); do
    TESTING=/media/dippi/Volume1/hpc_tests/testing/test-${i}.txt
    $COMMAND ${i} 1000 > $TESTING
    #OUTPUT=`diff /media/dippi/Volume1/hpc_tests/test-real-${i}-.txt $TESTING`
    diff /media/dippi/Volume1/hpc_tests/test-real-${i}.txt $TESTING

    if [[ $? != 0 ]]; then
        #echo $OUTPUT
        echo "Test non passato per $i"
        #exit 255
    fi
done