#!/bin/bash
array=(64 128 512 1024 4096 8192 16384 32786 65536)
f() { array=("${BASH_ARGV[@]}"); } # reverse arr

shopt -s extdebug
f "${array[@]}"
shopt -u extdebug

for T in ${array[@]}
do
counter=1
    echo $T >> result_threshold.txt
    while [ $counter -le 3 ]
    do
    prun -1 -np 1 ./parallel  -r -l 10000000 -p 32 -s 42 -t $T  >> result_threshold.txt
    ((counter++))
    done
done