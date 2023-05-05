#!/bin/bash

## TIME START
SECONDS=0

cwd=$(pwd)
echo $cwd

for file in segmented_compiled-*.fasta
do
f=${file%%.fasta}
echo $f
clustalw \
  -ALIGN \
  -INFILE=./${f}.fasta \
  -OUTFILE=./clustalw_${f}.fasta \
  -OUTPUT=FASTA
rm ${f}.dnd
done

## REPORT TIME             
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
echo "----- Done -----"
echo ""
