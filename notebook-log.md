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
Fastas to be aligned.

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

#### Write nexus file in R
 - Used write.nexus(), which is a function in APE

### RAxML-NG for maximum likelihood
RAxML-NG is purportedly a fast, modular ML program. It touts its speed and precision, but seems to fall short of IQ-Tree in tree inference accuracy. Still, it generally finds the best-scoring tree overall! Its speed and general precision makes it the optimal choice for my preliminary ML calculations.
Assumptions: mutational processes same at each branch; all sites evolve independently of eachother
(this is unlikely to be suitable for influenza in general, but this analysis will be a suitable first step toward a tree)

#### download
Downloaded raxml-ng_v1.1.0_linux_x86_64.zip from https://github.com/amkozlov/raxml-ng/releases/tag/1.1.0, as I had issues with the make step after pulling from GitHub 
```shell
## moved .zip file to ~/Code/
mkdir raxml-ng
mv raxml-ng*.zip raxml-ng/
unzip raxml-ng*.zip

./raxml-ng -v
 #(base) rieshunter@hries:~/Code/raxml-ng$ ./raxml-ng -v

 #RAxML-NG v. 1.1.0 released on 29.11.2021 by The Exelixis Lab.
 #Developed by: Alexey M. Kozlov and Alexandros Stamatakis.
 #Contributors: Diego Darriba, Tomas Flouri, Benoit Morel, Sarah Lutteropp, Ben Bettisworth.
 #Latest version: https://github.com/amkozlov/raxml-ng
 #Questions/problems/suggestions? Please visit: https://groups.google.com/forum/#!forum/raxml

 #System: Intel(R) Core(TM) i5-7360U CPU @ 2.30GHz, 2 cores, 7 GB RAM

## added to PATH, so I can call in myProject
# nano .bashrc
# PATH="$PATH:/home/rieshunter/Code/raxml-ng"
## reloaded terminal session
raxml-ng -v
 #<pre>(base) <font color="#8AE234"><b>rieshunter@hries</b></font>:<font color="#729FCF"><b>~/Code/myProject</b></font>$ raxml-ng -v

 #RAxML-NG v. 1.1.0 released on 29.11.2021 by The Exelixis Lab.
 #Developed by: Alexey M. Kozlov and Alexandros Stamatakis.
 #Contributors: Diego Darriba, Tomas Flouri, Benoit Morel, Sarah Lutteropp, Ben Bettisworth.
 #Latest version: https://github.com/amkozlov/raxml-ng
 #Questions/problems/suggestions? Please visit: https://groups.google.com/forum/#!forum/raxml

 #System: Intel(R) Core(TM) i5-7360U CPU @ 2.30GHz, 2 cores, 7 GB RAM
```

#### execute
```shell
#(base) rieshunter@hries:~/Code/myProject/data/HK_1/parsed_fa/muscle/raxml-ng$ 
raxml-ng --check --msa ../muscle_segmented_compiled-HA.fasta --model GTR+G
mv ../*raxml.* .
raxml-ng --check --msa *.phy --model GTR+G
 # in both logs, we see identical consensus-level sequences, which are expected for rapid flu transmission within a community. Acute infections transmit little diversity.

## time to infer the tree
raxml-ng --msa *.phy --prefix HK_1_HA --model GTR+G --seed 920
 # (base) rieshunter@hries:~/Code/myProject/data/HK_1/parsed_fa/muscle/raxml-ng$ ls
 #total 236K
 #-rw-rw-r-- 1 rieshunter rieshunter  130 Mar 30 20:13 HK_1_HA.raxml.bestModel
 #-rw-rw-r-- 1 rieshunter rieshunter 1.4K Mar 30 20:13 HK_1_HA.raxml.bestTree
 #-rw-rw-r-- 1 rieshunter rieshunter 1.3K Mar 30 20:13 HK_1_HA.raxml.bestTreeCollapsed
 #-rw-rw-r-- 1 rieshunter rieshunter  14K Mar 30 20:13 HK_1_HA.raxml.log
 #-rw-rw-r-- 1 rieshunter rieshunter  28K Mar 30 20:13 HK_1_HA.raxml.mlTrees
 #-rw-rw-r-- 1 rieshunter rieshunter 3.8K Mar 30 20:13 HK_1_HA.raxml.rba
 #-rw-rw-r-- 1 rieshunter rieshunter  28K Mar 30 20:13 HK_1_HA.raxml.startTree
 #-rw-rw-r-- 1 rieshunter rieshunter 4.7K Mar 30 20:13 muscle_segmented_compiled-HA.fasta.raxml.log
 #-rw-rw-r-- 1 rieshunter rieshunter  54K Mar 30 20:13 muscle_segmented_compiled-HA.fasta.raxml.reduced.phy
 #-rw-rw-r-- 1 rieshunter rieshunter 1.4K Mar 30 20:13 muscle_segmented_compiled-HA.fasta.raxml.reduced.phy.raxml.log
```

### MrBayes
 - Created shell script "Bayesian.sh" in ~/scrips/

#### Download (linux)
`sudo apt install mrbayes`
```
mb -v
## returns the following:
##  Version:    3.2.7a
##  Features:   Beagle readline
##  Host type:  x86_64-pc-linux-gnu (CPU: x86_64)
##  Compiler:   gnu 11.2.0
```

#### MrBayes Block
 - "mb_block.txt" in ~/scripts/
 - See "https://github.com/crsl4/phylogenetics-class/blob/master/exercises/notebook-log.md" for example used here

#### Append block to .nex and run
```
cat HK_1.nex mb_block.txt > HK_1-mb.nex

mb HK_1-mb.nex

```