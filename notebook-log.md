# Code for semester project
### Description of data
Source: https://elifesciences.org/articles/35962#content and https://www.ncbi.nlm.nih.gov/bioproject/PRJNA412631
Who: HIVE Cohort for respiratory virus surveillance
What: Flu sequences from 200 people in a test-negative study
Where: JT McCrone et al. at UMich
When: 2010-2015
Why: We know flu evolution is evolutionarily stochastic within-host but deterministic globally. Maybe the answer lies 
between-host?

Aims:
Evolution w/in hosts through serial time points
Evolution b/w hosts through household pairing

### Stage of data
Fastas to be aligned.

### Download data
Processed fastas retrieved from GitHub page:
https://github.com/lauringlab/Host_level_IAV_evolution

### Create concatenated fasta file for each segment
Created "concatenate.sh", which effectively grabs each segment from each sample and concatenates it to one file, per gene. 
Here's small sample:
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
ClustalW will be a good "Vanilla" test that is most strong in its use of weight-based scoring. Because so many of my 
sequences will be largely identical, this could be an important component that other software may not have. ClustalW also 
outputs a .dnd file which can be used to make a dendrogram. Might save me a step later!
Muscle will be a nice experiment on progressive alignment, which should be an improvement in accuracy and speed (yes, almost 
3x faster!) in comparison to ClustalW. While this is faster and perhaps more accurate, I'm unsure at this point how to test 
accuracy in comparison to other programs. How do we know which is more accurate?

Installed clustalw with `brew install brewsci/bio/clustal-w`
clustalw version: CLUSTAL 2.1 Multiple Sequence Alignments

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
    A parsimony ratchet (100 iterations) was applied to allow for quick and relatively accurate maximum parsimony 
calculations
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

### RAxML-NG for maximum likelihood
RAxML-NG is purportedly a fast, modular ML program. It touts its speed and precision, but seems to fall short of IQ-Tree in 
tree inference accuracy. Still, it generally finds the best-scoring tree overall! Its speed and general precision makes it 
the optimal choice for my preliminary ML calculations.
Assumptions: mutational processes same at each branch; all sites evolve independently of eachother
(this is unlikely to be suitable for influenza in general, but this analysis will be a suitable first step toward a tree)

#### download
Downloaded raxml-ng from brew with `brew install brewsci/bio/raxml-ng` on MacOS
```shell
./raxml-ng -v
 RAxML-NG v. 1.1-master released on 29.11.2021 by The Exelixis Lab.
 Developed by: Alexey M. Kozlov and Alexandros Stamatakis.
 Contributors: Diego Darriba, Tomas Flouri, Benoit Morel, Sarah Lutteropp, Ben Bettisworth.
 Latest version: https://github.com/amkozlov/raxml-ng
 Questions/problems/suggestions? Please visit: https://groups.google.com/forum/#!forum/raxml
 
 System: Intel(R) Core(TM) i5-8279U CPU @ 2.40GHz, 4 cores, 16 GB RAM
```

#### execute
```shell
## remove spaces in taxa names for cali09
sed -i.bak 's/parsed .*/parsed/g' muscle_segmented_compiled-HA.fasta
 # creates .bak backup of original

## remove Genbank and spaces from taxa names from perth
sed -i.bak 's/ GenBank.*//g' muscle_segmented_compiled-HA.fasta
 
#(base) rieshunter@hries:~/Code/myProject/data/HK_1/parsed_fa/muscle/raxml-ng$ 
raxml-ng --check --msa ../muscle_segmented_compiled-HA.fasta --model GTR+G
mv ../*raxml.* .
raxml-ng --check --msa *.phy --model GTR+G
 # in both logs, we see identical consensus-level sequences, which are expected for rapid flu transmission within a 
community. Acute infections transmit little diversity.

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

### MrBayes for bayesian inference
#### Rationale
MrBayes is a program for Bayesian inference using a Metropolis-Coupled Markov Chain Monte Carlo (MC-MCMC). Using posterior 
probabilities of trees, we can make phylogenetic inferences for a group of species. To approximate this, MrBayes utilizes 
MC-MCMC on phylogenetic trees. Importantly, we can adjust models and assumptions easily with MrBayes, allowing for the 
assessment of diverse parameters relatively easily. The major strength of MrBayes is its computational efficiency, largely 
due to Bayesian inference. Additionally, we can apply prior knowledge to a dataset to output a distribution of a 
parameter—not just a point estimate. The major weaknesses of MrBayes are largely in the user—the user's priors, models, and 
assumptions of both. Further, some may argue that the incorporation of beliefs on a dataset is inherently flawed. Still, 
MrBayes is clearly a useful tool, amongst many, to analyze phylogenetic relationships.

#### My assumptions
Proportions:  Beta
Branch len.:  log-normal
State freq.:  Dirichlet

#### Download (Mac)
```shell
brew tap brewsci/bio
brew install mrbayes

mb -v
## returns the following:
##  MrBayes, Bayesian Analysis of Phylogeny
##  
##  Version:   3.2.7a
##  Features:  SSE Beagle MPI
##  Host type: x86_64-apple-darwin22.1.0 (CPU: x86_64)
##  Compiler:  clang 14.0.0
```

#### MrBayes Block
 - "mb_block.txt" in ~/scripts/
 - See "https://github.com/crsl4/phylogenetics-class/blob/master/exercises/notebook-log.md" for example used here
 - Modified to use ML-determined outgroup MH8714_A

#### Create Nexus file (Mac)
 - Used run_clustalW_Nexus.sh script in ~/scripts/ on segmented_compiled files
 - Same code as run_clustalW.sh used earlier but outputs Nexus file
 - Copied file from working ~/data/HK_1/parsed_fa/clustalw/ folder to ~/results/ folder
```shell
mv clustalw_segmented_compiled-HA.nexus ./HK_1.nex
```
 - Trimmed tip labels:
```shell
sed -i.bu 's/HA//g' HK_1.nex
sed -i.bu 's/.removed.parsed//g' HK_1.nex
sed -i.bu 's/|//g' HK_1.nex
rm HK_1.nex.bu #intermediate file for OS X
```

#### Append block to .nex and run
```shell
cat HK_1.nex ../scripts/mb_block.txt > HK_1-mb.nex
mb HK_1-mb.nex
```
 - Performed one run at 10,000 runs, outputting every 100
    - See Terminal output: "TerminalOutput.txt"

#### Analyzed in Tracer v1.7.2
 - Imported HK_1-mb.nex.p to Tracer
 - First run at 10,000 runs didn't appear to reach convergence completely. See screenshots in ~/results/ngen_10k

### ASTRAL for coalescence
#### Rationale
ASTRAL is a program for multi-species coalescence, outputting a species tree from a set of gene trees. I do not plan to 
utilize ASTRAL or any other coalescence program in my work, so I will use the test data set outlined in astral-tutorial.md.
Main strengths: CLI and incorporation of priors is a pro. Statistics seem to be an 
improvement from ML phylogenies when gene tree discordance is not incredily high or low.
Assumptions: Substitution model priors are usually imperfect biologically but sufficient nonetheless.

#### Download and execute
```shell
git clone https://github.com/smirarab/ASTRAL.git
unzip Astral.5.7.8.zip 
cd Astral
java -jar astral.5.7.8.jar -i test_data/song_mammals.424.gene.tre 
```

Output:
```
================== ASTRAL ===================== 

This is ASTRAL version 5.7.8
Gene trees are treated as unrooted
424 trees read from test_data/song_mammals.424.gene.tre
index0
All output trees will be *arbitrarily* rooted at Chicken

======== Running the main analysis
Number of taxa: 37 (37 species)
Taxa: [Chicken, Marmoset, Orangutan, Human, Chimpanzee, Gorilla, Macaque, Galagos, Mouse_Lemur, Tree_Shrew, Mouse, Rat, 
Kangaroo_Rat, Guinea_Pig, Squirrel, Tarsier, Rabbit, Pika, Microbat, Megabat, Horse, Dolphin, Cow, Alpaca, Pig, Dog, Cat, 
Shrew, Hedgehog, Lesser_Hedgehog_Tenrec, Hyrax, Elephant, Sloth, Armadillos, Platypus, Opossum, Wallaby]
Taxon occupancy: {Rat=424, Tarsier=424, Dolphin=424, Rabbit=424, Macaque=424, Pika=424, Alpaca=424, Shrew=424, Sloth=424, 
Tree_Shrew=424, Kangaroo_Rat=424, Armadillos=424, Chimpanzee=424, Horse=424, Dog=424, Human=424, Lesser_Hedgehog_Tenrec=424, 
Microbat=424, Platypus=424, Wallaby=424, Cow=424, Pig=424, Marmoset=424, Megabat=424, Hedgehog=424, Mouse=424, 
Guinea_Pig=424, Mouse_Lemur=424, Cat=424, Hyrax=424, Elephant=424, Chicken=424, Orangutan=424, Opossum=424, Galagos=424, 
Squirrel=424, Gorilla=424}
Number of gene trees: 424
0 trees have missing taxa
Calculating quartet distance matrix (for completion of X)
Species tree distances calculated ...
Building set of clusters (X) from gene trees 
------------------------------
gradient0: 1933
Number of Clusters after addition by distance: 1933
calculating extra bipartitions to be added at level 1 ...
Adding to X using resolutions of greedy consensus ...
Limit for sigma of degrees:975
polytomy size limit : 4
discarded polytomies:  [3, 3, 4, 4]
Threshold 0.0:
Threshold 0.01:
Threshold 0.02:
Threshold 0.05:
Threshold 0.1:
Threshold 0.2:
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 1933
Threshold 0.3333333333333333:
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 1933
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 1933
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 1933
max k is :0
Number of Clusters after addition by greedy: 1933
gradient0 in heuristiic: 1933
partitions formed in 0.733 secs
Dynamic Programming starting after 0.733 secs
Using tree-based weight calculation.
Using polytree-based weight calculation.
Polytree max score: 28003080
Polytree building time: 0.202 seconds.
Number of quartet trees in the gene trees: 28003080
Size of largest cluster: 37
Greedy score: 24862814
estimationFactor: 1.1263037241078182
Sub-optimal score: 25489533
Total Number of elements weighted: 3372
Normalized score (portion of input quartet trees satisfied before correcting for multiple individuals): 0.9115752624354179
Optimization score: 25526915
Optimal tree inferred in 1.417 secs.
(Chimpanzee,(Human,(Gorilla,(Orangutan,(Macaque,(Marmoset,(Tarsier,((Galagos,Mouse_Lemur),((Tree_Shrew,((Rabbit,Pika),(Squirrel,(Guinea_Pig,(Kangaroo_Rat,(Mouse,Rat)))))),((((Chicken,Platypus),(Opossum,Wallaby)),((Sloth,Armadillos),(Lesser_Hedgehog_Tenrec,(Hyrax,Elephant)))),((Shrew,Hedgehog),((Microbat,Megabat),((Horse,(Dog,Cat)),(Alpaca,(Pig,(Dolphin,Cow))))))))))))))));
Final quartet score is: 25526915
Final normalized quartet score is: 0.9115752624354179
Extended species tree:
(Chicken,(Platypus,((Opossum,Wallaby)1:4.9534768802562965,(((Sloth,Armadillos)1:4.3332364705044455,(Lesser_Hedgehog_Tenrec,(Hyrax,Elephant)1:1.5072311535707152)1:3.0324267685203896)1:0.11752571582479492,(((Shrew,Hedgehog)1:1.0128082767511968,((Microbat,Megabat)1:1.5529600021930556,((Horse,(Dog,Cat)1:2.9057840368910477)0.9:0.0641150813890718,(Alpaca,(Pig,(Dolphin,Cow)1:1.3549566469016843)1:0.6392263251934307)1:3.5727535641858488)1:0.11521870556469156)1:0.43083761402314513)1:1.9714179086025256,((Tree_Shrew,((Rabbit,Pika)1:3.449399483480025,(Squirrel,(Guinea_Pig,(Kangaroo_Rat,(Mouse,Rat)1:5.045850200387328)1:0.6153942169908233)1:0.14915704136865612)1:1.6583346999346553)1:0.843253312273408)0.91:0.08610422555073367,((Galagos,Mouse_Lemur)1:2.4173938646853443,(Tarsier,(Marmoset,(Macaque,(Orangutan,(Gorilla,(Human,Chimpanzee)1:0.6434676443273446)1:2.3363022618905545)1:2.5809965452562977)1:2.687856981495734)1:4.347341076686)1:0.6511829814609134)1:2.1216694610241857)1:2.3363739622649566)1:1.5631449093243042)1:3.3544501191225256)1:0.9262330413691204));
(Chicken,(Platypus,((Opossum,Wallaby)1:4.9534768802562965,(((Sloth,Armadillos)1:4.3332364705044455,(Lesser_Hedgehog_Tenrec,(Hyrax,Elephant)1:1.5072311535707152)1:3.0324267685203896)1:0.11752571582479492,(((Shrew,Hedgehog)1:1.0128082767511968,((Microbat,Megabat)1:1.5529600021930556,((Horse,(Dog,Cat)1:2.9057840368910477)0.9:0.0641150813890718,(Alpaca,(Pig,(Dolphin,Cow)1:1.3549566469016843)1:0.6392263251934307)1:3.5727535641858488)1:0.11521870556469156)1:0.43083761402314513)1:1.9714179086025256,((Tree_Shrew,((Rabbit,Pika)1:3.449399483480025,(Squirrel,(Guinea_Pig,(Kangaroo_Rat,(Mouse,Rat)1:5.045850200387328)1:0.6153942169908233)1:0.14915704136865612)1:1.6583346999346553)1:0.843253312273408)0.91:0.08610422555073367,((Galagos,Mouse_Lemur)1:2.4173938646853443,(Tarsier,(Marmoset,(Macaque,(Orangutan,(Gorilla,(Human,Chimpanzee)1:0.6434676443273446)1:2.3363022618905545)1:2.5809965452562977)1:2.687856981495734)1:4.347341076686)1:0.6511829814609134)1:2.1216694610241857)1:2.3363739622649566)1:1.5631449093243042)1:3.3544501191225256)1:0.9262330413691204):0.0); 
Weight calculation took 0.369462758 secs
ASTRAL finished in 1.941 secs
```
