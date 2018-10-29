#!/bin/bash

#sh ./tpm.sh DRR2345097

rm -f list list2 buscos.txt

for table in $(find . -name full_table_*$1*); do
    echo $table;
done | sort -k1 >> list

grep Complete $(sed -n 1p list) \
| cut -f1 \
| grep -f - $(sed -n 2p list) \
| grep Complete \
| cut -f1 \
| grep -f - $(sed -n 3p list) \
| grep Complete \
| cut -f1 \
| grep -f - $(sed -n 4p list) \
| grep Complete \
| cut -f1 \
| grep -f - $(sed -n 5p list) \
| grep Complete \
| cut -f1 \
| sort -R \
| head >> list2


for file in $(find ./ -name full_table_*$1* 2> /dev/null)
do
    for busco in $(cat list2)
    do
        echo $(basename $file) $(grep $busco $file) | awk '{print $1 "\t" $2 "\t" $4}'
    done
done | sort -k1 >> buscos.txt


x=$(wc -l list2 | cut -d " " -f1)
x1=$(expr $x + 1)
x2=$(echo "$x * 2" | bc)
x3=$(expr $x2 + 1)
x4=$(echo "$x * 3" | bc)
x5=$(expr $x4 + 1)
x6=$(echo "$x * 4" | bc)
x7=$(expr $x6 + 1)
x8=$(echo "$x * 5" | bc)

sed -n 1,"$x"p buscos.txt | cut -f2
echo "orthomerged"
for i in $(sed -n 1,"$x"p buscos.txt | cut -f3)
do
    grep -w "$i" "$(find . -name salmon_orthomerged_"$1")"/quant.sf | cut -f4
done

echo "Shannon"
for i in $(sed -n "$x1","$x2"p buscos.txt | cut -f3)
do
    grep -w $i $(find . -name salmon_shannon_"$1")/quant.sf | cut -f4
done

echo "SPAdes55"
for i in $(sed -n "$x3","$x4"p buscos.txt | cut -f3)
do
    grep -w $i $(find . -name salmon_sp55_"$1")/quant.sf | cut -f4
done

echo "SPAdes75"
for i in $(sed -n "$x5","$x6"p buscos.txt | cut -f3)
do
    grep -w $i $(find . -name salmon_sp75_"$1")/quant.sf | cut -f4
done

echo "Trinity"
for i in $(sed -n "$x7","$x8"p buscos.txt | cut -f3)
do
    grep -w $i $(find . -name salmon_trin_"$1")/quant.sf | cut -f4
done

#rm buscos.txt list list2
