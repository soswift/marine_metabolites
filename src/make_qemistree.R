## Generate a two-way heatmap of metabolites and microbes
## Microbes are arrange by phylogeny, metabolites are arranged by qemistree classifcation
## Match up all the tip names with heatmap data (correlations)
## Pass phylogenies to complexHeatmap as dendrograms
library(ape)
library(data.table)
library(dendextend)
library(ComplexHeatmap)
library(phylogram)
library(tidyr)

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

qem_id_map <- setNames(tip_data$featureID, tip_data$id)


# Read MMVEC and correlation tables -------------------------------

# mmvec data (ranks for metabolite microbe correlation)
mmvec_table <- fread("data/raw/mmvec/Ranks_result.tsv",
                    key = "featureid")

mmvec_table[ , featureid := sub("metabolite",
                               "",
                               featureid)]
setnames(mmvec_table,
         "featureid",
         "featureID")


# read Craig correlation table (long format showing pairwise comparisons)
correl_raw <-read.csv("data/processed/correl.csv")

# clean up and make into wide matrix
cors <- correl_raw[ ,c("Label", "weight")]

cors$id  <- gsub(".+d-(\\d+).+", "\\1", correl_raw$Label, perl = T)
cors$otu <- gsub(".*(Otu\\d+).+","\\1", correl_raw$Label, perl = T)
cors$otu <- gsub("Otu(\\d{4}$)","Otu0\\1", cors$otu)

cors_wide <- reshape(cors[ ,-1], idvar = "id", timevar = "otu", direction = "wide")

row.names(cors_wide) <- cors_wide$id
colnames(cors_wide) <- sub("weight.","",colnames(cors_wide))

cors_mat <- as(cors_wide[ , -1 ], "matrix")

# Read and Organize Metabolite Qemistree Tree ---------------------------------------------------

# read in qemistree tree file 
qemistree_raw <- read.tree(file = "data/raw/qemistree/qemistree.tree")

# update tip labels to match the rest of our metabolite data
qemistree_raw$tip.label <- unname(qem_id_map[qemistree_raw$tip.label])

qem_tip_data <- tip_data[featureID %in% qemistree_raw$tip.label]


# Read and Organize Microbe FastTree Phylogenetic Tree --------------------------------------------

# read in fastree file
fastree_raw <- read.tree("data/raw/FastTree_100.tre")

# Heatmap of MMVEC Values-----------------------------------------------
# subset mmvec data to only include ids that are in qemistree
mmvec_mat <- as(mmvec_table[ , .SD, .SDcols = !"featureID"], "matrix")
row.names(mmvec_mat) <- mmvec_table$featureID

mmvec_mat <- mmvec_mat[row.names(mmvec_mat) %in% qemistree_clean$tip.label , ]


# normalize  mmvec data
scale_dat <- function(mat_dat){
  mat_dat <- t(scale(t(mat_dat)))
  mat_dat[is.na(mat_dat)] <- 0
  return(mat_dat)
}

z_mmvec_mat <- scale_dat(mmvec_mat)
z_mmvec_mat[1:5, 1:5]

# generate_phymap() makes a complex Heatmap that 'clusters' using phylogenetic trees (or equivalent qemistree)
# takes metabolites as rows, otus as columns

generate_phymap <- function(micro_tree = fastree_raw, chem_tree = qemistree_raw, correlation_mat = cors_mat){
 
  ## microbe
  # subset microbe tree to ASVs that are present in correlation matrix and vice versa
  mat_asvs <- colnames(correlation_mat)
  
  micro_tree_asvs <-mat_asvs[ mat_asvs %in% micro_tree$tip.label]
  
  micro_tree_clean <- keep.tip(micro_tree, micro_tree_asvs)
  correlation_mat <- correlation_mat[ , micro_tree_asvs]
  
  # microbe tree as dendrogram
  micro_dendro <- as.dendrogram.phylo(micro_tree_clean)
  
  ## chem
  # subset metabolite tree to features in the correlation matrix and vice versa
  if(!is.null(chem_tree)){
  mat_feats <- row.names(correlation_mat)
  
  chem_tree_feats <-
    mat_feats[mat_feats %in% chem_tree$tip.label]
  
  chem_tree_clean <- keep.tip(chem_tree, chem_tree_feats)
  correlation_mat <- correlation_mat[chem_tree_feats , ]
  
  
  # metabolite tree as dendrogram
  chem_dendro  <- as.dendrogram.phylo(chem_tree_clean)
  } else {
    chem_dendro <- as.dendrogram( hclust ( dist( correlation_mat ) ) )
  }
  
  
  # generate heatmap
  ht <- Heatmap(
    correlation_mat,
    cluster_rows = chem_dendro,
    cluster_columns = micro_dendro,
    # row parameters (metabolite)
    row_dend_width = unit(4, "in"),
    row_names_gp = gpar(fontsize = 6),
    show_row_names = F,
    # column parameters (microbe)
    column_dend_height = unit(2, "in"),
    column_names_gp = gpar(fontsize = 6),
    show_column_names = F,
    # colors
    na_col = "white"
  )
  return(ht)
}

# generate heatmaps with and without qemistree

ht_cor <- generate_phymap(correlation_mat = cors_mat)

ht_cor_clust <- generate_phymap(correlation_mat = cors_mat,
                                chem_tree = NULL)

ht_mmvec <- generate_phymap(correlation_mat = z_mmvec_mat)
ht_mmvec_clust <- generate_phymap(correlation_mat = z_mmvec_mat, 
                                  chem_tree = NULL)

# plot phylogeny heatmap
png(
  filename = "output/heatmap/metabolite_microbe_phylo_heatmap.png",
  width = 30, height = 30, res = 300, units = "in" )

print(ht_cor)

dev.off()

# plot cluster heatmap

png(
  filename = "output/heatmap/metabolite_microbe_clust_heatmap.png",
  width = 30,
  height = 25,
  res = 300,
  units = "in" )
print(ht_cor_clust)
dev.off()
