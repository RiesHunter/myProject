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
t_coffee ${cwd}/${f}.fasta
mv ${cwd}/${f}.aln ${cwd}/tcoffee_MSA_${f}.aln
mv ${cwd}/${f}.dnd ${cwd}/tcoffee_GUIDE_TREE_${f}.dnd
mv ${cwd}/${f}.html ${cwd}/tcoffee_MSA_${f}.html
done

mkdir tcoffee; mv tcoffee* tcoffee

## REPORT TIME             
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
echo "----- Done -----"
echo ""