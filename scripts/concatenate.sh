#!/bin/bash

pwd

for file in *.fasta
do
f=${file%%.fasta}

## NP
sed -n '/>*NP/p' ./${f}.fasta > name
sed "s/NP/NP|${f}/g" name > name.temp
cat name.temp >> segmented_compiled-NP
sed -n '/>*NP/,/>/p' ./${f}.fasta > seq
sed '$ d' seq > seq.temp
sed '1d' seq.temp > seq
cat seq >> segmented_compiled-NP
echo -ne "[#--------] ${f}\r"

## NS
sed -n '/>*NS/p' ./${f}.fasta > name
sed "s/NS/NS|${f}/g" name > name.temp
cat name.temp >> segmented_compiled-NS
sed -n '/>*NS/,/>/p' ./${f}.fasta > seq
#sed '$ d' seq > seq.temp        #removed because it's last in fasta
sed '1d' seq.temp > seq
cat seq >> segmented_compiled-NS
echo -ne "[##-------] ${f}\r"

## MP
sed -n '/>*M/p' ./${f}.fasta > name
sed "s/M/M|${f}/g" name > name.temp
cat name.temp >> segmented_compiled-M
sed -n '/>*M/,/>/p' ./${f}.fasta > seq
sed '$ d' seq > seq.temp
sed '1d' seq.temp > seq
cat seq >> segmented_compiled-M
echo -ne "[###------] ${f}\r"

## PA
sed -n '/>*PA/p' ./${f}.fasta > name
sed "s/PA/PA|${f}/g" name > name.temp
cat name.temp >> segmented_compiled-PA
sed -n '/>*PA/,/>/p' ./${f}.fasta > seq
sed '$ d' seq > seq.temp
sed '1d' seq.temp > seq
cat seq >> segmented_compiled-PA
echo -ne "[####-----] ${f}\r"

## PB1
sed -n '/>*PB1/p' ./${f}.fasta > name
sed "s/PB1/PB1|${f}/g" name > name.temp
cat name.temp >> segmented_compiled-PB1
sed -n '/>*PB1/,/>/p' ./${f}.fasta > seq
sed '$ d' seq > seq.temp
sed '1d' seq.temp > seq
cat seq >> segmented_compiled-PB1
echo -ne "[#####----] ${f}\r"

## PB2
sed -n '/>*PB2/p' ./${f}.fasta > name
sed "s/PB2/PB2|${f}/g" name > name.temp
cat name.temp >> segmented_compiled-PB2
sed -n '/>*PB2/,/>/p' ./${f}.fasta > seq
sed '$ d' seq > seq.temp
sed '1d' seq.temp > seq
cat seq >> segmented_compiled-PB2
echo -ne "[######---] ${f}\r"

## NA
sed -n '/>*NR/p' ./${f}.fasta > name
sed "s/NR/NR|${f}/g" name > name.temp
cat name.temp >> segmented_compiled-NA
sed -n '/>*NR/,/>/p' ./${f}.fasta > seq
sed '$ d' seq > seq.temp
sed '1d' seq.temp > seq
cat seq >> segmented_compiled-NR
echo -ne "[#######--] ${f}\r"

## HA
sed -n '/>*HA/p' ./${f}.fasta > name
sed "s/HA/HA|${f}/g" name > name.temp
cat name.temp >> segmented_compiled-HA
sed -n '/>*HA/,/>/p' ./${f}.fasta > seq
sed '$ d' seq > seq.temp
sed '1d' seq.temp > seq
cat seq >> segmented_compiled-HA
echo -ne "[########-] ${f}\r"

## Whole Genome
echo ">Whole_Genome|${f}" > name.temp
sed '/>/d' ./${f}.fasta > seq
cat name.temp >> segmented_compiled-Whole_Genome
cat seq >> segmented_compiled-Whole_Genome
echo -ne "[#########] ${f}\r"

echo ""
done


mv segmented_compiled-NP segmented_compiled-NP.fasta
mv segmented_compiled-NS segmented_compiled-NS.fasta
mv segmented_compiled-M segmented_compiled-M.fasta
mv segmented_compiled-PA segmented_compiled-PA.fasta
mv segmented_compiled-PB2 segmented_compiled-PB2.fasta
mv segmented_compiled-PB1 segmented_compiled-PB1.fasta
mv segmented_compiled-NR segmented_compiled-NR.fasta
mv segmented_compiled-HA segmented_compiled-HA.fasta
mv segmented_compiled-Whole_Genome segmented_compiled-Whole_Genome.fasta

#mkdir segmented_compiled_fastas
#mv segmented_compiled-* segmented_compiled_fastas

rm seq
rm name
rm name.temp
rm seq.temp