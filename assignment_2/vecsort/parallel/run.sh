#!/bin/bash
array=(1 2 5 10 20 50 100 500 1000 2000) 
for chunk in ${array[@]}
do
counter=1
    echo $chunk >> result4096_chunk_dynamic.txt
    while [ $counter -le 3 ]
    do
    prun -1 -np 1 ./parallel  -r -l 1000 -n 10 -x 100000 -s 42 -p 32 -c $chunk -t 4096 >> result4096_chunk_dynamic.txt
    ((counter++))
    done
done



