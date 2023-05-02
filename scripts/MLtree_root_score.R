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
plot_cali09 <- comparePhylo(ClustalW, Muscle, plot = T, use.edge.length = T)

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
plot_perth <- comparePhylo(ClustalW, Muscle, plot = T, use.edge.length = T)


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

setwd(dir_s); getwd(); dir()
ggsave("Fig_speed.pdf", plot_run_stats,
       width = 6.5, height = 3, 
       units = "in", dpi = 320)

ggsave("Fig_logLikelihood.pdf", plot_likelihood,
       width = 6.5, height = 3, 
       units = "in", dpi = 320)