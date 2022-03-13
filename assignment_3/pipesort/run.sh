#bin/bash
#array1=(1 2 4 8 16 32 50 64 100 200 500)
array1=(700 1000 2000 4000 8000)
array2=(4090)
for P in ${array1[@]}
do
for T in ${array2[@]}
do
counter=1
    echo $P $T >> result_test.txt
    while [ $counter -le 3 ]
    do
    prun -1 -np 1 ./pipesort  -l $T -s 42 -b $P   >> result_test.txt
    ((counter++))
    done
done
done
