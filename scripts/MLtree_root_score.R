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
if(!require(phytools)){
  install.packages("phytools", dependencies = T)
  library(phytools)
}

#### Import data ####
dir_s <- paste("/Users/rieshunter/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/_Admin/Thesis/Coursework/Phylo/myProject/results")
dir_HK <- paste("/Users/rieshunter/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/_Admin/Thesis/Coursework/Phylo/myProject/results/HK")
dir_p <- paste("/Users/rieshunter/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/_Admin/Thesis/Coursework/Phylo/myProject/results/perth")
dir_c <- paste("/Users/rieshunter/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/_Admin/Thesis/Coursework/Phylo/myProject/results/cali09")

#### Functions ####
clean_cols <- function(x) {
  names(x) <- gsub(x = names(x), pattern = ".removed.parsed", replacement = "")  
  names(x) <- gsub(x = names(x), pattern = "HA\\|", replacement = "")  
  return(x)
}
clean_tips <- function(x) {
  x$tip.label <- gsub(".removed.parsed", "", x$tip.label)
  x$tip.label <- gsub("HA\\|", "", x$tip.label)
  is.rooted(x)
  return(x)
}
root_tree <- function (x) {
  #x <- root(x, outgroup = 1, resolve.root = TRUE)
  x <- midpoint.root(x)
  print(is.rooted(x))
  return(x)
}
likelihood <- function(x, y) {
  temp <- pml(x, phydat, site.rate = "gamma", model = "GTR")
  temp2 <- data.frame("logLik" = temp$logLik, "group" = y)
  return(temp2)
}

#### HK ####
## import
setwd(dir_HK); getwd(); dir()

#ClustalW through RAxML-NG
bestTree <- paste("HK_clustalw.raxml.bestTree")
tr_bestTree <- read.tree(bestTree)
phydat <- read.phyDat(file="clustalw_segmented_compiled-HA.fasta", format="fasta")
phydat <- clean_cols(phydat)
tr_bestTree <- clean_tips(tr_bestTree)
ClustalW <- root_tree(tr_bestTree)
plot(ClustalW)
df_likelihood <- likelihood(ClustalW, "HK_c_r")

#Muscle through RAxML-NG
bestTree <- paste("HK_muscle.raxml.bestTree")
tr_bestTree <- read.tree(bestTree)
phydat <- read.phyDat(file="muscle_segmented_compiled-HA.fasta", format="fasta")
phydat <- clean_cols(phydat)
tr_bestTree <- clean_tips(tr_bestTree)
Muscle <- root_tree(tr_bestTree)
plot(Muscle)
df_likelihood <- rbind(likelihood(Muscle, "HK_m_r"),df_likelihood)

## compare trees
plot_HK <- comparePhylo(ClustalW, Muscle, plot = T, use.edge.length = T)

#### cali09 ####
## import
setwd(dir_c); getwd(); dir()

#ClustalW through RAxML-NG
bestTree <- paste("cali09_clustalw.raxml.bestTree")
tr_bestTree <- read.tree(bestTree)
phydat <- read.phyDat(file="clustalw_segmented_compiled-HA.fasta", format="fasta")
phydat <- clean_cols(phydat)
tr_bestTree <- clean_tips(tr_bestTree)
ClustalW <- root_tree(tr_bestTree)
plot(ClustalW)
df_likelihood <- rbind(likelihood(ClustalW, "cali09_c_r"),df_likelihood)

#Muscle through RAxML-NG
bestTree <- paste("cali09_muscle.raxml.bestTree")
tr_bestTree <- read.tree(bestTree)
phydat <- read.phyDat(file="muscle_segmented_compiled-HA.fasta", format="fasta")
phydat <- clean_cols(phydat)
tr_bestTree <- clean_tips(tr_bestTree)
Muscle <- root_tree(tr_bestTree)
plot(Muscle)
df_likelihood <- rbind(likelihood(Muscle, "cali09_m_r"),df_likelihood)

## compare trees
comparePhylo(ClustalW, Muscle, plot = T, use.edge.length = T)

#### perth ####
## import
setwd(dir_p); getwd(); dir()

#ClustalW through RAxML-NG
bestTree <- paste("perth_clustalw.raxml.bestTree")
tr_bestTree <- read.tree(bestTree)
phydat <- read.phyDat(file="clustalw_segmented_compiled-HA.fasta", format="fasta")
phydat <- clean_cols(phydat)
tr_bestTree <- clean_tips(tr_bestTree)
ClustalW <- root_tree(tr_bestTree)
plot(ClustalW)
df_likelihood <- rbind(likelihood(ClustalW, "perth_c_r"),df_likelihood)

#Muscle through RAxML-NG
bestTree <- paste("perth_muscle.raxml.bestTree")
tr_bestTree <- read.tree(bestTree)
phydat <- read.phyDat(file="muscle_segmented_compiled-HA.fasta", format="fasta")
phydat <- clean_cols(phydat)
tr_bestTree <- clean_tips(tr_bestTree)
Muscle <- root_tree(tr_bestTree)
plot(Muscle)
df_likelihood <- rbind(likelihood(Muscle, "perth_m_r"),df_likelihood)

## compare trees
comparePhylo(ClustalW, Muscle, plot = T, use.edge.length = T)

#### Plot ####
setwd(dir_s); getwd(); dir()