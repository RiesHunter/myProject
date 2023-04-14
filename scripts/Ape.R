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
dir_s <- paste("/Users/rieshunter/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/_Admin/Thesis/Coursework/Phylo/myProject/results")
dir_muscle <- paste("/Users/rieshunter/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/_Admin/Thesis/Coursework/Phylo/myProject/data/HK_1/parsed_fa/muscle")
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

#### Root ML tree ####
dir <- paste("/Users/rieshunter/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/_Admin/Thesis/Coursework/Phylo/myProject/data/HK_1/parsed_fa/muscle/raxml-ng")
setwd(dir); getwd(); dir()
## .bestTree
bestTree <- paste("HK_1_HA.raxml.bestTree")
tr_bestTree <- read.tree(text = c('((HA|MH8714_A.removed.parsed:0.000572,(((((((HA|HS1566_B.removed.parsed:0.001146,(HA|HS1579_B.removed.parsed:0.015078,(HA|PC1A.removed.parsed:0.000596,HA|PC2B.removed.parsed:0.001123):0.001087):0.001817):0.000001,(HA|MH8924_A.removed.parsed:0.001128,(HA|HS1364_B.removed.parsed:0.000001,HA|HS1341_B.removed.parsed:0.000572):0.000934):0.000232):0.000001,HA|HS1409_B.removed.parsed:0.000001):0.000001,(HA|HS1345_B.removed.parsed:0.000572,HA|HS1417_B.removed.parsed:0.000572):0.000001):0.000001,(((HA|HS1425_B.removed.parsed:0.000572,HA|HS1357_B.removed.parsed:0.000001):0.001727,HA|HS1465_B.removed.parsed:0.001146):0.000577,(HA|HS1294_B.removed.parsed:0.000572,(((HA|MH8589_A.removed.parsed:0.001731,(HA|MH8453_B.removed.parsed:0.001145,HA|HS1518_B.removed.parsed:0.001720):0.000001):0.000001,HA|HS1562_B.removed.parsed:0.000001):0.000001,HA|HS1457_B.removed.parsed:0.001146):0.001722):0.000001):0.000001):0.000001,(HA|HS1573_B.removed.parsed:0.001146,HA|HS1532_B.removed.parsed:0.001722):0.000001):0.000001,(HA|MH9045_A.removed.parsed:0.001743,HA|MH9547_B.removed.parsed:0.001151):0.002302):0.000572):0.000001,((HA|MH8922_A.removed.parsed:0.001719,(HA|HS1530_B.removed.parsed:0.000001,(HA|HS1377_B.removed.parsed:0.001145,HA|HS1484_B.removed.parsed:0.001145):0.000001):0.000572):0.000001,(HA|HS1427_B.removed.parsed:0.000001,HA|HS1455_B.removed.parsed:0.000572):0.001147):0.000001,HA|HS1402_B.removed.parsed:0.000001);'))

## clean labels
tr_bestTree$tip.label <- gsub(".removed.parsed", "", tr_bestTree$tip.label)
tr_bestTree$tip.label <- gsub("HA\\|", "", tr_bestTree$tip.label)
is.rooted(tr_bestTree)

tr_bestTree_root <- root(tr_bestTree, outgroup = 1, resolve.root = TRUE)
 # rooting to MH8714_A, the first node
is.rooted(tr_bestTree_root)
plot(tr_bestTree_root)

## tree to nexus format
setwd(dir_s); getwd(); dir()
write.nexus(tr_bestTree_root, file = "HK_1.nex")

## double-check that it looks good
plot(read.nexus("HK_1.nex"))
