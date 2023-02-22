#!/bin/bash
for ((i=100; i<=10000; i=$i+100)); do
    for ((j=50; j<=1000; j=$j+50)); do
        bin/sph $i $j > /media/dippi/Volume1/hpc_tests/test-real-$i-$j.txt
    done
done
