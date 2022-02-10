#!/bin/bash

sheet=$1
myarray=$(grep Lane $sheet | sed "s/,/ /g")
my_array=($myarray)
value='index'

for i in "${!my_array[@]}"; do
   if [[ "${my_array[$i]}" = "${value}" ]]; then
       index=${i}
   fi
done

for line in $(cat $sheet | grep -v Lane | grep -v "\[Data\]"); do
    linearray=($(echo $line | sed "s/,,/,null,/g" | sed "s/^,/null,/g" | sed "s/,/ /g" ))
    seq=${linearray[$index]}
    seqlen=$(expr length $seq)
done

echo $seqlen

