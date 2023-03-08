## Clear Global Environment
rm(list = ls())

#### Session prep ####
## Install packages and load libraries as required
install.packages("tidyverse", dependencies = T)
library(tidyverse)
if(!require(ggplot2)){
  install.packages("ggplot2", dependencies = T)
  library(ggplot2)
}

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  library(BiocManager)
}
BiocManager::install("ggtree")
library(ggtree)

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
dna <- fasta2DNAbin(file="http://adegenet.r-forge.r-project.org/files/usflu.fasta")

#### Analyze data ####
D <- dist.dna(dna, model="TN93")
tre <- nj(D)
tre <- ladderize(tre)

#### Plot data ####
plot(tre, cex=.6)
title("A simple NJ tree")
