for ((i=100; i<=20000; i=$i+100)); do
    for ((j=50; j<=1000; j=$j+50)); do
        ./sph $i $j > test-$i-$j.txt
    done
done