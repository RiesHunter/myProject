#!/bin/bash

pwd

for file in *.fasta
do
f=${file%%.fasta}
## HA
sed -n '/>*HA/p' ./${f}.fasta > name
sed "s/HA/HA|${f}/g" name > name.temp
cat name.temp >> segmented_compiled-HA
sed -n '/>*HA/,/>/p' ./${f}.fasta > seq
sed '$ d' seq > seq.temp
sed '1d' seq.temp > seq
cat seq >> segmented_compiled-HA
echo -ne "${f}\r"
echo ""
done

mv segmented_compiled-HA segmented_compiled-HA.fasta

#mkdir segmented_compiled_fastas
#mv segmented_compiled-* segmented_compiled_fastas

rm seq
rm name
rm name.temp
rm seq.temp