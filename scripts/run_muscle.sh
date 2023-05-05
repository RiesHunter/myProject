#!/bin/bash

## TIME START
SECONDS=0

cwd=$(pwd)
echo $cwd

for file in segmented_compiled-*.fasta
do
f=${file%%.fasta}
echo $f
muscle \
  -in ./${f}.fasta \
  -out ./muscle_${f}.fasta
done

## REPORT TIME             
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
echo "----- Done -----"
echo ""
