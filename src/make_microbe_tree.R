# Quick visualizations of microbial phylogeny
## use ggtree to reconstruct qemistree and highlight features
library(ggplot2)
library(ggtree)
library(tidytree)
library(ape)
library(data.table)
library(ggstance)
library(ggtreeE)

# Read in tree ------------------------------------------------------


# read in tree file 
all_reads<-read.FASTA("data/raw/all_esv.fasta")

Ntip(tree_raw)
  # plot tree
ggtree(tree_raw, layout =  )


