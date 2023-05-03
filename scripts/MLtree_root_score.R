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
if(!require(ggtree)){
  install.packages("ggtree", dependencies = T)
  library(ggtree)
}

#### Import data ####
dir_s <- paste("/Users/rieshunter/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/_Admin/Thesis/Coursework/Phylo/myProject/figures")
dir_HK <- paste("/Users/rieshunter/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/_Admin/Thesis/Coursework/Phylo/myProject/results/HK")
dir_p <- paste("/Users/rieshunter/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/_Admin/Thesis/Coursework/Phylo/myProject/results/perth")
dir_c <- paste("/Users/rieshunter/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/_Admin/Thesis/Coursework/Phylo/myProject/results/cali09")
dir <- paste("/Users/rieshunter/Library/CloudStorage/GoogleDrive-hries@wisc.edu/Shared drives/TCF lab/Current Lab Members/Hunter_Ries/_Admin/Thesis/Coursework/Phylo/myProject/results")
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
root_tree <- function (x, y) {
  x <- root(x, outgroup = y, 
            resolve.root = TRUE)
  #x <- midpoint.root(x)
  print(is.rooted(x))
  return(x)
}
likelihood <- function(x, y) {
  temp <- pml(x, phydat, site.rate = "gamma", model = "GTR")
  temp2 <- data.frame("logLik" = temp$logLik, "group" = y)
  return(temp2)
}
treelength <- function(x, y) {
  temp <- sum(node.depth.edgelength(x))
  temp2 <- data.frame("treelength" = temp, "group" = y)
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
ClustalW <- root_tree(tr_bestTree, "MH8156_A")
df_likelihood <- likelihood(ClustalW, "HK_c_r")
df_treelength <- treelength(ClustalW, "HK_c_r")

#Muscle through RAxML-NG
bestTree <- paste("HK_muscle.raxml.bestTree")
tr_bestTree <- read.tree(bestTree)
phydat <- read.phyDat(file="muscle_segmented_compiled-HA.fasta", format="fasta")
phydat <- clean_cols(phydat)
tr_bestTree <- clean_tips(tr_bestTree)
Muscle <- root_tree(tr_bestTree,"MH8156_A")
df_likelihood <- rbind(likelihood(Muscle, "HK_m_r"),df_likelihood)
df_treelength <- rbind(treelength(ClustalW, "HK_m_r"),df_treelength)

## compare trees
plot_HK <- comparePhylo(ClustalW, Muscle, plot = T, use.edge.length = T)

wasp.cophylo<-cophylo(ClustalW,Muscle)
plot(wasp.cophylo,link.type="curved",
     link.lwd=4,
     link.lty="solid",
     link.col=make.transparent("red", 0.25))
title(main="Hong Kong",font.main=3)

ggtree(Muscle) + 
  geom_tiplab(size = 2, offset = 0.00125, hjust = 1) + 
  theme_tree2() +
  labs(title = "Hong Kong\nMuscle")
ggtree(ClustalW) + 
  geom_tiplab(size = 2, offset = 0.00125, hjust = 1) + 
  theme_tree2() +
  labs(title = "Hong Kong\nClustalW")

#### cali09 ####
## import
setwd(dir_c); getwd(); dir()

#ClustalW through RAxML-NG
bestTree <- paste("cali09_clustalw.raxml.bestTree")
tr_bestTree <- read.tree(bestTree)
phydat <- read.phyDat(file="clustalw_segmented_compiled-HA.fasta", format="fasta")
phydat <- clean_cols(phydat)
tr_bestTree <- clean_tips(tr_bestTree)
ClustalW <- root_tree(tr_bestTree, "1336")
df_likelihood <- rbind(likelihood(ClustalW, "cali09_c_r"),df_likelihood)
df_treelength <- rbind(treelength(ClustalW, "cali09_c_r"),df_treelength)

#Muscle through RAxML-NG
bestTree <- paste("cali09_muscle.raxml.bestTree")
tr_bestTree <- read.tree(bestTree)
phydat <- read.phyDat(file="muscle_segmented_compiled-HA.fasta", format="fasta")
phydat <- clean_cols(phydat)
tr_bestTree <- clean_tips(tr_bestTree)
Muscle <- root_tree(tr_bestTree, "1336")
df_likelihood <- rbind(likelihood(Muscle, "cali09_m_r"),df_likelihood)
df_treelength <- rbind(treelength(ClustalW, "cali09_m_r"),df_treelength)

## compare trees
plot_cali09 <- comparePhylo(ClustalW, Muscle, plot = T, use.edge.length = T)

wasp.cophylo<-cophylo(ClustalW,Muscle)
plot(wasp.cophylo,link.type="curved",
     link.lwd=4,
     link.lty="solid",
     link.col=make.transparent("red", 0.25))
title(main="Cali09",font.main=3)

ggtree(Muscle) + 
  geom_tiplab(size = 2, offset = 0.00125, hjust = 1) + 
  theme_tree2() +
  labs(title = "Cali09\nMuscle")
ggtree(ClustalW) + 
  geom_tiplab(size = 2, offset = 0.00125, hjust = 1) + 
  theme_tree2() +
  labs(title = "Cali09\nClustalW")

#### perth ####
## import
setwd(dir_p); getwd(); dir()

#ClustalW through RAxML-NG
bestTree <- paste("perth_clustalw.raxml.bestTree")
tr_bestTree <- read.tree(bestTree)
phydat <- read.phyDat(file="clustalw_segmented_compiled-HA.fasta", format="fasta")
phydat <- clean_cols(phydat)
tr_bestTree <- clean_tips(tr_bestTree)
ClustalW <- root_tree(tr_bestTree, "Perth_pool")
df_likelihood <- rbind(likelihood(ClustalW, "perth_c_r"),df_likelihood)
df_treelength <- rbind(treelength(ClustalW, "perth_c_r"),df_treelength)

#Muscle through RAxML-NG
bestTree <- paste("perth_muscle.raxml.bestTree")
tr_bestTree <- read.tree(bestTree)
phydat <- read.phyDat(file="muscle_segmented_compiled-HA.fasta", format="fasta")
phydat <- clean_cols(phydat)
tr_bestTree <- clean_tips(tr_bestTree)
Muscle <- root_tree(tr_bestTree, "Perth_pool")
df_likelihood <- rbind(likelihood(Muscle, "perth_m_r"),df_likelihood)
df_treelength <- rbind(treelength(ClustalW, "perth_m_r"),df_treelength)

## compare trees
plot_perth <- comparePhylo(ClustalW, Muscle, plot = T, use.edge.length = T)

wasp.cophylo<-cophylo(ClustalW,Muscle)
plot(wasp.cophylo,link.type="curved",
     link.lwd=4,
     link.lty="solid",
     link.col=make.transparent("red", 0.25))
title(main="Perth",font.main=3)

ggtree(Muscle) + 
  geom_tiplab(size = 2, offset = 0.00125, hjust = 1) + 
  theme_tree2() +
  labs(title = "Perth\nMuscle")
ggtree(ClustalW) + 
  geom_tiplab(size = 2, offset = 0.00125, hjust = 1) + 
  theme_tree2() +
  labs(title = "Perth\nClustalW")

#### Stats ####
## df_run
setwd(dir); getwd(); dir()
df_run <- as.data.frame(read_tsv("run_data.tsv"))
df_run$Seqs_per_sec <- df_run$Seqs/df_run$Secs
df_run <- separate(df_run, "Group", c("dataset", "alignment"), sep = "_")
df_run$alignment <- gsub("C", "ClustalW", df_run$alignment)
df_run$alignment <- gsub("M", "MUSCLE", df_run$alignment)
df_run$da <- as.factor(paste(df_run$dataset, df_run$alignment, sep = "_"))
df_run$Seqs_per_sec <- round(df_run$Seqs_per_sec, 3)

## df_likelihood
df_likelihood <- separate(df_likelihood, "group", c("dataset", "alignment", "ml"), sep = "_")
df_likelihood$dataset <- as.factor(df_likelihood$dataset)
df_likelihood$alignment <- gsub("m", "MUSCLE", df_likelihood$alignment)
df_likelihood$alignment <- gsub("c", "ClustalW", df_likelihood$alignment)
df_likelihood$alignment <- as.factor(df_likelihood$alignment)
df_likelihood$da <- ""
df_likelihood$da <- paste(df_likelihood$dataset, df_likelihood$alignment, sep = "_")
df_likelihood$da <- as.factor(df_likelihood$da)
df_likelihood$logLik <- round(df_likelihood$logLik, 4)

## df_treelength
df_treelength <- separate(df_treelength, "group", c("dataset", "alignment", "ml"), sep = "_")
df_treelength$dataset <- as.factor(df_treelength$dataset)
df_treelength$alignment <- gsub("m", "MUSCLE", df_treelength$alignment)
df_treelength$alignment <- gsub("c", "ClustalW", df_treelength$alignment)
df_treelength$alignment <- as.factor(df_treelength$alignment)
df_treelength$da <- ""
df_treelength$da <- paste(df_treelength$dataset, df_treelength$alignment, sep = "_")
df_treelength$da <- as.factor(df_treelength$da)
df_treelength$treelength <- round(df_treelength$treelength, 4)

#### Plot ####
setwd(dir_s); getwd(); dir()
## plot_run_stats
plot_run_stats <- ggplot() + 
  geom_point(data = df_run, 
             aes(x = dataset, y = Seqs_per_sec,
                 color = alignment, group = alignment),
             position = position_dodge2(width = .2)) + 
  geom_text(data = df_run, 
            aes(x = dataset, y = Seqs_per_sec,
                color = alignment, group = alignment, label = Seqs_per_sec),
            position = position_dodge2(width = .2), vjust=-0.5, size = 2.5,
            show.legend = FALSE) + 
  scale_y_continuous(limits = c(0, 12)) + 
  labs(x = "", y = "Sequences per second") + 
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.title = element_blank())

## plot_likelihood
plot_likelihood <- ggplot() + 
  geom_point(data = df_likelihood, 
             aes(x = dataset, y = logLik, 
                 color = alignment, group = alignment),
             position = position_dodge2(width = .2)) + 
  geom_text(data = df_likelihood, 
            aes(x = dataset, y = logLik, 
                color = alignment, group = alignment, label = logLik),
            position=position_dodge2(0.9), vjust=-0.5, size = 2.5,
            show.legend = FALSE) + 
  scale_y_reverse(limits = c(-3400, -4000)) +
  labs(x = "", y = "Log10 Likelihood of Tree") + 
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.title = element_blank())

## plot_treelength
plot_treelength <- ggplot() + 
  geom_point(data = df_treelength, 
             aes(x = dataset, y = treelength, 
                 color = alignment, group = alignment),
             position = position_dodge2(width = .2)) + 
  geom_text(data = df_treelength, 
            aes(x = dataset, y = treelength, 
                color = alignment, group = alignment, label = treelength),
            position=position_dodge2(0.9), vjust=-0.5, size = 2.5,
            show.legend = FALSE) + 
  scale_y_continuous(limits = c(0, 2)) +
  labs(x = "", y = "Tree length") + 
  theme_bw() + 
  theme(legend.position = "bottom",
        legend.title = element_blank())
  

plots <- plot_grid(plot_run_stats, 
                   plot_likelihood, 
                   plot_treelength,
                   ncol = 1,
                   align = "v")
ggsave("Fig3.pdf", plots,
       width = 5, height = 10, 
       units = "in", dpi = 320)


