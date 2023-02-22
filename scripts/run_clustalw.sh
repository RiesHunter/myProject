#!/bin/bash

## TIME START
SECONDS=0

#data_dir="/home/rieshunter/GitHub/myProject/data/aligned"
cwd=$(pwd)
echo $cwd

for file in segmented_compiled-*.fasta
do
f=${file%%.fasta}
echo $f
clustalw \
  -ALIGN \
  -INFILE=${cwd}/${f}.fasta \
  -OUTFILE=${cwd}/clustalw_${f}.fasta \
  -OUTPUT=FASTA
done

mkdir clustalw
mv clustalw* clustalw
mv *.dnd clustalw

## REPORT TIME             
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
echo "----- Done -----"
echo ""