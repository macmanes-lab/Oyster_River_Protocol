#!/bin/bash

#./numbers.sh
#this script finds SCOs found in each assembly ()

rm -f list list2 buscos.txt

for table in $(find . -name full_table*); do 
    echo $table >> list;
done

grep Complete $(sed -n 1p list) \
| cut -f1 \
| grep -f - $(sed -n 2p list) \
| grep Complete \
| cut -f1 \
| grep -f - $(sed -n 3p list) \
| grep Complete \
| cut -f1 \
| sort -R \
| head >> list2


for file in $(find ./ -name full_table* 2> /dev/null)
do
    for busco in $(cat list2)
    do
        echo $(basename $file) $(grep $busco $file) | awk '{print $1 "\t" $2 "\t" $4}' | tee -a buscos.txt
    done
done

rm list list2
