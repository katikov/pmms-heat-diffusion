#!bin/bash

array=(10000 100000 1000000 10000000 100000000) 
for M in a d r
do
for T in ${array[@]}
do
counter=1
    echo $M $T >> result.txt
    while [ $counter -le 3 ]
    do
    prun -1 -np 1 ./sequential  -$M -l $T -p 32 -s 42 >> result.txt
    ((counter++))
    done
done
done



