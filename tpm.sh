#!/bin/bash

cat buscos.txt | sed -n 1,10p | cut -f2
echo $1
for i in $(sed -n 1,10p buscos.txt | cut -f3)
do
    grep $i $1/quant.sf | cut -f4
done

echo $2
for i in $(sed -n 11,20p buscos.txt | cut -f3)
do
    grep $i $2/quant.sf | cut -f4
done

echo $3
for i in $(sed -n 21,30p buscos.txt | cut -f3)
do
    grep $i $3/quant.sf | cut -f4
done

echo $4
for i in $(sed -n 31,40p buscos.txt | cut -f3)
do
    grep $i $4/quant.sf | cut -f4
done

rm buscos.txt
