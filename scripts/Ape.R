## Clear Global Environment
rm(list = ls())

#### Session prep ####
## Install packages and load libraries as required
if(!require(ape)){
  install.packages("ape", dependencies = T)
  library(ape)
}
if(!require(adegenet)){
  install.packages("adegenet", dependencies = T)
  library(adegenet)
}
if(!require(phangorn)){
  install.packages("phangorn", dependencies = T)
  library(phangorn)
}

#### Import data ####
dir_s <- paste("/Users/rieshunter/Documents/GitHub/myProject/data/HK_1/parsed_fa")
dir_muscle <- paste("/Users/rieshunter/Documents/GitHub/myProject/data/HK_1/parsed_fa/muscle")
setwd(dir_muscle); getwd(); dir()
dna_muscle_HA <- fasta2DNAbin(file="muscle_segmented_compiled-HA.fasta")
phydat_muscle_HA <- read.phyDat(file="muscle_segmented_compiled-HA.fasta", format="fasta")

#### Analyze data ####
## Neighbor-joining distance-based tree
dm  <- dist.ml(dna_muscle_HA, model = "JC69")
treeUPGMA  <- upgma(dm)
treeNJ  <- NJ(dm)

## Parsimony-based tree
parsimony(treeUPGMA, phydat_muscle_HA)
parsimony(treeNJ, phydat_muscle_HA)
treeRatchet  <- pratchet(phydat_muscle_HA, trace = 0, minit=100)
parsimony(treeRatchet, phydat_muscle_HA)
treeRatchet  <- acctran(treeRatchet, phydat_muscle_HA)
treeRatchet  <- di2multi(treeRatchet)
if(inherits(treeRatchet, "multiPhylo")){
  treeRatchet <- unique(treeRatchet)
}

#### Plot data ####
## Neighbor-joining distance-based trees
plot(treeNJ, cex=.4, main = "NJ JC69 distance-based tree")
plot(treeUPGMA, cex=.4, main = "UPGMA JC69 distance-based tree")

## Parsimony-based tree
plotBS(midpoint(treeRatchet), 
       type="phylogram", 
       main = "UPGMA JC69 parsimony-based tree",
       cex=.4)
add.scale.bar()