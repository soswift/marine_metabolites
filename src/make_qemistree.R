## Generate a two-way heatmap of metabolites and microbes
## Microbes are arrange by phylogeny, metabolites are arranged by qemistree classifcation
## Match up all the tip names with heatmap data (correlations)
## Pass phylogenies to complexHeatmap as dendrograms
library(ggplot2)
library(ggtree)
library(tidytree)
library(ape)
library(data.table)
library(ggstance)
library(ggtreeE)
library(dendextend)
library(ComplexHeatmap)
library(phylogram)

detach("package:speedyseq")
unloadNamespace("speedyseq")
detach("package:phyloseq")
unloadNamespace("phyloseq")

# Read Qemistree Tables -------------------------------------------------

# tables of data related to metabolite qemistree 
class_tips <- fread("data/raw/qemistree/labels.tsv",
                    key = "id")
  
color_tips <- fread("data/raw/qemistree/colors.tsv", 
                    key = "id")

barplot_tips <- fread("data/raw/qemistree/barplots.tsv",
                      key = "id")

# table that provides key to matching qemistree tip labels to other metabolomics data(mmvec etc.)
tip_to_feature <- fread("data/raw/qemistree/Fingerprints to features.tsv",
                        key = "id")

# additional feature information (networks, classifications, etc.)
peak_data <- fread("data/processed/table_exports/all_marine_metabolite_tax_flat_table.csv")
setnames(peak_data,"V1", "id")

# put the qemistree tips together with information we have on the individual metabolites

tip_data <- merge(class_tips, color_tips)
tip_data <- merge(tip_data, barplot_tips)
tip_data <- merge(tip_to_feature, tip_data)

qem_tip_data <- tip_data[id %in% qemistree_raw$tip.label]

qem_id_map <- setNames(tip_data$featureID, tip_data$id)


# Read MMVEC tables -------------------------------

# mmvec data (ranks for metabolite microbe correlation)
mmvec_table <- fread("data/raw/mmvec/Ranks_result.tsv",
                    key = "featureid")

mmvec_table[ , featureid := sub("metabolite",
                               "",
                               featureid)]
setnames(mmvec_table,
         "featureid",
         "featureID")

# transform mmvec scores to matrix
mmvec_mat <- as(mmvec_table[ , .SD, .SDcols = !"featureID"], "matrix")

# set the mmvec rownames to ids (which match qemistree) instead of featureID
# TODO make sure all tips are in mmvec data!
row.names(mmvec_mat) <- mmvec_table$featureID


# Read and Organize Metabolite Qemistree Tree ---------------------------------------------------

# read in qemistree tree file 
qemistree_raw <- read.tree(file = "data/raw/qemistree/qemistree.tree")

# update tip labels to match the rest of our metabolite data
qemistree_raw$tip.label <- unname(qem_id_map[qemistree_raw$tip.label])


# get list of metabolite features that are in mmvec data
mmvec_feats <- row.names(mmvec_mat)

# check if there are more features in the mmvec table than in qemistree
length(mmvec_feats)
Ntip(qemistree_raw)

# get mmvec feats that are present in qemistree
qem_mmvec_feats <- mmvec_feats[mmvec_feats %in% qemistree_raw$tip.label]

# subset tree to tips that we have mmvec data for
qemistree_clean <- keep.tip(qemistree_raw, qem_mmvec_feats)
Ntip(qemistree_clean)

# Read and Organize Microbe FastTree Phylogenetic Tree --------------------------------------------

# read in fastree file
fastree_raw <- read.tree("data/raw/FastTree_100.tre")

# subset microbe tree to ASVs that were run with mmvec
mmvec_asvs <- colnames(mmvec_mat)
fastree_clean <- keep.tip(fastree_raw, mmvec_asvs)


# Heatmap of MMVEC Values-----------------------------------------------

# get microbe tree as dendrogram
micro_dendro <- as.dendrogram.phylo(fastree_clean)

# get metabolite qemistree as dendrogram
chem_dendro  <- as.dendrogram.phylo(qemistree_clean)

# subset mmvec data to only include ids that are in qemistree
mmvec_mat <- mmvec_mat[row.names(mmvec_mat) %in% qemistree_clean$tip.label , ]

# generate heatmap
ht <- Heatmap(
              mmvec_mat,
              cluster_rows = chem_dendro,
              cluster_columns = micro_dendro,
              # row parameters (metabolite)
              row_dend_width = unit(4, "in"),
              row_names_gp = gpar(fontsize = 6),
              show_row_names = F,
              # column parameters (microbe)
              column_dend_height = unit(2, "in"),
              column_names_gp = gpar(fontsize = 6),
              show_column_names = F)
        

# rownames of mmvec table match 'featureID' in tip_data, 
png(
  filename = "output/heatmap/metabolite_microbe_phylo_heatmap.png",
  width = 30,
  height = 30,
  res = 300,
  units = "in" )
print(ht)
dev.off()

# normal clustering heatmap
ht2 <- Heatmap( mmvec_mat)


png(
  filename = "output/heatmap/metabolite_microbe_clust_heatmap.png",
  width = 30,
  height = 25,
  res = 300,
  units = "in" )
print(ht2)
dev.off()
