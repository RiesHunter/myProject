# Code for semester project
### Description of data
Source: https://elifesciences.org/articles/35962#content and https://www.ncbi.nlm.nih.gov/bioproject/PRJNA412631
Who: HIVE Cohort for respiratory virus surveillance
What: Flu sequences from 200 people in a test-negative study
Where: JT McCrone et al. at UMich
When: 2010-2015
Why: We know flu evolution is evolutionarily stochastic within-host but deterministic globally. Maybe the answer lies between-host?

Aims:
Evolution w/in hosts through serial time points
Evolution b/w hosts through household pairing

### Stage of data
Raw SRA reads in need of QC and alignment

### Download data
Processed fastas retrieved from GitHub page:
https://github.com/lauringlab/Host_level_IAV_evolution

### Create concatenated fasta file for each segment
Created "concatenate.sh", which effectively grabs each segment from each sample and concatenates it to one file, per gene. Here's small sample:
```shell
sed -n '/>*NP/p' ./${f}.fasta > name
sed "s/NP/NP|${f}/g" name > name.temp
cat name.temp >> segmented_compiled-NP
sed -n '/>*NP/,/>/p' ./${f}.fasta > seq
sed '$ d' seq > seq.temp
sed '1d' seq.temp > seq
cat seq >> segmented_compiled-NP
echo -ne "[#--------] ${f}\r"
```

### Align HK_1 with ClustalW and Muscle
ClustalW will be a good "Vanilla" test that is most strong in its use of weight-based scoring. Because so many of my sequences will be largely identical, this could be an important component that other software may not have. ClustalW also outputs a .dnd file which can be used to make a dendrogram. Might save me a step later!
Muscle will be a nice experiment on progressive alignment, which should be an improvement in accuracy and speed (yes, almost 3x faster!) in comparison to ClustalW. While this is faster and perhaps more accurate, I'm unsure at this point how to test accuracy in comparison to other programs. How do we know which is more accurate?

Created "run_clustalw.sh":
```shell
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
```

Benchmarking 10 flu genomes:
 - Muscle: 1 minute and 14 seconds
 - Clustalw: 3 minutes and 11 seconds
 - Tcoffee: [soon-to-be-clocked]

### Utilize ape, phangorn, and adegenet to construct trees
In R, I imported the muscle alignment for HA
The alignment was imported with both "fasta2DNAbin" and "read.phyDat"
The "fasta2DNAbin" command was used for a UPGMA distance-based trees
  Pairwise distances were calculated using JC69
  Distance matrices were converted to trees with UPGMA and NJ
The "read.phyDat" command was used for a parsimony-based tree
  Parsimony of all trees were calculated:
    ```shell
      > parsimony(treeUPGMA, phydat_muscle_HA)
    [1] 97
    > parsimony(treeNJ, phydat_muscle_HA)
    [1] 90
    > parsimony(treeRatchet, phydat_muscle_HA)
    [1] 90
    ```
  treeRatchet is the parsimony tree
    A parsimony ratchet (100 iterations) was applied to allow for quick and relatively accurate maximum parsimony calculations
    acctran was applied to allow for polytomies (this dataset has a bunch of them!)
    di2multi collapses branches of zero length

Chosen algorithms:
  JC69
    - Incredibly simple
    - Assumes equal base freq.
    - Assumes equal base mutation rates
    - Only parameter is µ, the overall substitution rate

Figures:
![treeNJ](figures/treeNJ.tiff)
![treeUPGMA](figures/treeUPGMA.tiff)
![treeRatchet](figures/treeRatchet.tiff)

