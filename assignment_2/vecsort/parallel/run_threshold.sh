#!/bin/bash
array=(2 8 32 128 512 1024 4096 8192 16384 32786 65536)
#array=(32786 65536)
for T in ${array[@]}
do
counter=1
    echo $T >> result_threshold_guided.txt
    while [ $counter -le 3 ]
    do
    prun -1 -np 1 ./parallel  -r -l 1000 -n 10 -x 100000 -s 42 -p 32 -c 1 -t $T  >> result_threshold_guided.txt
    ((counter++))
    done
done
