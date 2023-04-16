#!/bin/bash
for ((i=100; i<=1000; i=$i+100)); do
    bin/sph $i 1000 > /media/dippi/Volume1/hpc_tests/test-real-$i.txt
    echo "Done $i"
done
