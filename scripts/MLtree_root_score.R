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
dir_HK <- paste("/Users/rieshunter/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/_Admin/Thesis/Coursework/Phylo/myProject/data/HK")
dir_p <- paste("/Users/rieshunter/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/_Admin/Thesis/Coursework/Phylo/myProject/data/perth")
dir_c <- paste("/Users/rieshunter/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/_Admin/Thesis/Coursework/Phylo/myProject/data/cali09")

#### Functions ####
clean_tips <- function(x) {
  x$tip.label <- gsub(".removed.parsed", "", x$tip.label)
  x$tip.label <- gsub("HA\\|", "", x$tip.label)
  is.rooted(x)
  return(x)
}
root_tree <- function (x) {
  root(x, outgroup = 1, resolve.root = TRUE)
  is.rooted(x)
  return(x)
}

#### Root ML tree ####
### HK
## import
setwd(dir_HK); getwd(); dir()

### cali09
## import
setwd(dir_c); getwd(); dir()

### perth
## import
setwd(dir_p); getwd(); dir()
#ClustalW through RAxML-NG
bestTree <- paste("HK_1_HA.raxml.bestTree")
tr_bestTree <- read.tree(text = c('((HA|MH8714_A.removed.parsed:0.000572,(((((((HA|HS1566_B.removed.parsed:0.001146,(HA|HS1579_B.removed.parsed:0.015078,(HA|PC1A.removed.parsed:0.000596,HA|PC2B.removed.parsed:0.001123):0.001087):0.001817):0.000001,(HA|MH8924_A.removed.parsed:0.001128,(HA|HS1364_B.removed.parsed:0.000001,HA|HS1341_B.removed.parsed:0.000572):0.000934):0.000232):0.000001,HA|HS1409_B.removed.parsed:0.000001):0.000001,(HA|HS1345_B.removed.parsed:0.000572,HA|HS1417_B.removed.parsed:0.000572):0.000001):0.000001,(((HA|HS1425_B.removed.parsed:0.000572,HA|HS1357_B.removed.parsed:0.000001):0.001727,HA|HS1465_B.removed.parsed:0.001146):0.000577,(HA|HS1294_B.removed.parsed:0.000572,(((HA|MH8589_A.removed.parsed:0.001731,(HA|MH8453_B.removed.parsed:0.001145,HA|HS1518_B.removed.parsed:0.001720):0.000001):0.000001,HA|HS1562_B.removed.parsed:0.000001):0.000001,HA|HS1457_B.removed.parsed:0.001146):0.001722):0.000001):0.000001):0.000001,(HA|HS1573_B.removed.parsed:0.001146,HA|HS1532_B.removed.parsed:0.001722):0.000001):0.000001,(HA|MH9045_A.removed.parsed:0.001743,HA|MH9547_B.removed.parsed:0.001151):0.002302):0.000572):0.000001,((HA|MH8922_A.removed.parsed:0.001719,(HA|HS1530_B.removed.parsed:0.000001,(HA|HS1377_B.removed.parsed:0.001145,HA|HS1484_B.removed.parsed:0.001145):0.000001):0.000572):0.000001,(HA|HS1427_B.removed.parsed:0.000001,HA|HS1455_B.removed.parsed:0.000572):0.001147):0.000001,HA|HS1402_B.removed.parsed:0.000001);'))
phydat_muscle <- read.phyDat(file="clustalw_segmented_compiled-HA.fasta", format="fasta")

tr_bestTree <- clean_tips(tr_bestTree)
perth_root <- root_tree(tr_bestTree)
plot(perth_root)
perth_parsimony  <- parsimony(perth_root, phydat_muscle, method = "fitch", site = "pscore")


## tree to nexus format
setwd(dir_s); getwd(); dir()