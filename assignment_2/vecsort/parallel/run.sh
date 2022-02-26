#!bin/bash
array1=(2 4 8 16 32)
array2=(4096 8192 16384) 
for P in ${array1[@]}
do
for T in ${array2[@]}
do
counter=1
k=2
    echo $P $T >> result_combine.txt
    while [ $counter -le 3 ]
    do
    prun -1 -np 1 ./parallel  -r -l 1000 -n 10 -x 100000 -s 42 -p $P -o $P/k -t $T  >> result_combine.txt
    ((counter++))
    done
done
done
