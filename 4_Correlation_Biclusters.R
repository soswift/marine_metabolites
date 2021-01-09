## Generate a two-way heatmap of metabolites and microbes
## Microbes are arrange by phylogeny, metabolites are arranged by qemistree classifcation
## Match up all the tip names with heatmap data (correlations)
## Pass phylogenies to complexHeatmap as dendrograms
library(ape)
library(data.table)
library(dendextend)
library(ComplexHeatmap)
library(phylogram)
library(circlize)

# source functions
source("src/biclust_helper_functions.R")

# Read Qemistree Tables -------------------------------------------------
tip_data_file <- "data/raw/qemistree/Fingerprints to features.tsv"
peak_data_file <- "data/processed/table_exports/all_marine_metabolite_tax_flat_table.csv"
chem_abund_file <- "data/processed/table_exports/paired_marine_metabolite_abundance_flat_table.csv"

asv_tax_file <- "data/processed/table_exports/paired_marine_microbe_tax_flat_table.csv"
micro_abund_file <- "data/processed/table_exports/paired_marine_microbe_abundance_flat_table.csv"

micro_tree_file <- "data/raw/FastTree_100.tre"
chem_tree_file <- "data/raw/qemistree/qemistree.tree"

sample_data_file <- "data/processed/table_exports/paired_marine_microbe_sample_flat_table.csv"

mmvec_file <- "data/raw/mmvec/Ranks_result.tsv"

# key to matching qemistree tip labels to other metabolomics data(mmvec etc.)
tip_to_feature <- fread(tip_data_file, key = "id")

# additional feature information (networks, classifications, etc.)metabolites
peak_data <- fread(peak_data_file)
setnames(peak_data,"V1", "featureID")
peak_data[ , featureID:= gsub("id_","",featureID)]

# Merge the qemistree outputs with metabolite peak information
tip_data <- merge(tip_to_feature, peak_data, by = "featureID")
tip_data$componentindex <- sub(" +", "", tip_data$componentindex)
qem_id_map <- setNames(tip_data$featureID, tip_data$id)

# Read and Organize Metabolite Qemistree Tree ---------------------------------------------------
# read in qemistree tree file 
qemistree_raw <- read.tree(file = chem_tree_file)

# update tip labels to 'featureID' to match the rest of our metabolite data
qemistree_raw$tip.label <- unname(qem_id_map[qemistree_raw$tip.label])
qem_tip_data <- tip_data[featureID %in% qemistree_raw$tip.label]
qemistree_raw <- drop.tip(qemistree_raw, tip_data[class == "unclassified", featureID])

# Read and Organize Microbe FastTree Phylogenetic Tree --------------------------------------------
# read in fastree file
fastree_raw <- read.tree(micro_tree_file)

# Read in ASV and Metabolite Abundance Data-----------------------
# asv taxonomy data
asv_table <- fread(asv_tax_file,
                   header = T)
setnames(asv_table, "V1", "OTU_ID")

# asv relative abundance
micro_abund <- read.csv(micro_abund_file,
                        row.names = 1)
micro_abund <- as.matrix(micro_abund)
colnames(micro_abund) <- sub("X","",colnames(micro_abund))

# sample data
sample_dat <- read.csv(sample_data_file)

# metabolite relative abundance
chem_abund <- read.csv(chem_abund_file,
                       row.names = 1)
chem_abund  <- as.matrix(chem_abund)
colnames(chem_abund) <- sub("X","",colnames(chem_abund))
row.names(chem_abund) <- sub("id_","",row.names(chem_abund))

# match order of microbe samples and metabolite samples
micro_abund <- micro_abund[ , colnames(chem_abund)]

# check samples match
all(colnames(micro_abund) == colnames(chem_abund))

# get sample type sums for microbes and metabolites (sourced from biclust_helpers.R)
# use these for barplots on the heatmap and to identify realtionsips betwen sample types
micro_sums <- get_sums(micro_abund,
                        group_col = "sample_type",
                        id_col = "sample_barcode",
                        meta = sample_dat,
                        new_name = "OTU_ID")

chem_sums <- get_sums(chem_abund,
                        group_col = "sample_type",
                        id_col = "sample_barcode",
                        meta = sample_dat,
                       new_name = "featureID")

# Read MMVEC and correlation tables -------------------------------

# MMVEC data (matrix showing pairwise mmvec scores)
mmvec_table <- fread(mmvec_file,
                     key = "featureid")

mmvec_table[ , featureid := sub("metabolite",
                                "",
                                featureid)]
setnames(mmvec_table,
         "featureid",
         "featureID")

# Filter Matrix For Heatmap-----------------------------------------------
## Subset correlation matrices to match qemistree and filter out uninteresting metabolites

# transform MMVEC scores to matrix
mmvec_mat <- as( mmvec_table[ , .SD, .SDcols = !"featureID"], "matrix")
row.names(mmvec_mat) <- mmvec_table$featureID

# Subset mmvec data to only include ids that are in qemistree and ASV metadata
mmvec_mat <- mmvec_mat[row.names(mmvec_mat) %in% qemistree_raw$tip.label ,
                       colnames(mmvec_mat) %in% asv_table$OTU_ID]
mmvec_mat <- mmvec_mat[apply(mmvec_mat, 1, median) >= 2 , ]

# scale mmvec data
z_mmvec_mat <- scale_dat(mmvec_mat)

# Generate Categorical Matrix -------------------------------------------------------------------
# use get_top_type() to figure out the top type for each microbe and metabolite (based on relative abundance)
# type_match() takes this information and generates a categorical matrix
# each cell in the matrix gets a value based on two criteria:
#1) is the correlation score above a specified cutoff (e.g. 1 or 0.4)
#2) were the microbe and metabolite most abundant in the same sample type or different sample types?

# calculate sample type with highest summed relative abundance
micro_top_type <- get_top_type(micro_sums)
chem_top_type  <- get_top_type(chem_sums)

# find matches
z_cat_mat <- type_match(mi_type = micro_top_type,
                        me_type = chem_top_type,
                        cor_mat = z_mmvec_mat,
                        cutoff = 1)

# Generate Heatmaps ----------------------------------------------
# see generate_phymap() in biclust_helper_functions.R for more information
# a lot of information (e.g. trees, sample data, etc.) is passed to this function via default settings
ht_cat_mmvec <- generate_phymap(correlation_mat = z_cat_mat)
save_heatmap(ht_cat_mmvec, "sample_type_heatmap_zmmvec_qemistree")

# Matching Spearman to MMVE
ht_cat_mmvec <- generate_phymap(correlation_mat = cat_match_mat)
save_heatmap(ht_cat_mmvec, "sample_type_heatmap_matching_spearman_qemistree")

# Zoomed In Heatmap -------------------------------------

# single example
# select interesting groups of microbes and metabolites
cool_network <- "118" # remove as contaminant
cool_order <- c("BD2-11_terrestrial_group_or","Chitinophagales")
cool_class <- "Glycerophospholipids"

# Plot glycerophospholipds and prenol lipids

# subset by microbe orders that show algae/coral difference and metabolite class 
glyc_mat <- sub_cor_mat(cor_mat = z_cat_mat,
                        chem_col = "class",
                        micro_col = "order",
                        chem_val = cool_class,
                        micro_val = cool_order)

# generate heatmap
ht_glyc <- generate_phymap(correlation_mat = glyc_mat,
                           box_plot = T)
  
save_heatmap(ht_glyc, "gly-pho-lip_bicluster")


# loop through other interesting chemical classes

cool_classes <- c(
  "Fatty Acyls",
  "Carboxylic acids and derivatives",
  "Prenol lipids",
  "Steroids and steroid derivatives",
  "Benzene and substituted derivatives",
  "Organooxygen compounds"
)

# make it easy to plot a bunch of heatmaps for subsets of the metabolite data
ht_wrapper <-function(cl, chem_col){
  cl_mat <- sub_cor_mat(cor_mat = z_cat_mat,
                        chem_col = chem_col,
                        micro_col = NULL,
                        chem_val = cl,
                        micro_val = NULL)
  # generate heatmap
  ht_cl <- generate_phymap( correlation_mat = cl_mat, box_plot = T, )
  save_heatmap(ht_cl, paste0(cl, "_bicluster"))
}

lapply(cool_classes, ht_wrapper, chem_col = "class")

## cool networks 
# Based on Craig's stats, not necessarily specific to sample type
# all have library matches
cool_networks <- c(114, 5, 18, 40, 52, 69, 74, 107, 135, 144, 146, 147, 191, 201, 241, 242, 289, 380, 610, 617, 734, 801)

# Single sample type networks
cca_networks <- c(106, 127, 416, 501, 529, 54, 88 )
coral <- c(103)

# Pull out networks that are actually in the correlation matrix
z_mmvec_nets <- tip_data[tip_data$featureID %in% row.names(z_mmvec_mat),
                         .N, by =  componentindex][order(N, decreasing = T)]

nets_we_have <- z_mmvec_nets[N > 3 & componentindex != "-1", componentindex]

# which cool networks are in the mmvec selected networks?


lapply(nets_we_have[nets_we_have %in% cool_networks],
       ht_wrapper, chem_col = "componentindex")

# Identifying Noteworthy Metabolites ------------------------------------------

## merge all of our data for microbe metabolite relationships together
# 1) coerce all data to long format and standardize column names
# 2) merge all data together and filter pairwise associations by various criteria

# z-scored mmvec table
z_mmvec_table <- data.table(z_mmvec_mat, keep.rownames = T)
z_mmvec_table_long <- melt.data.table(z_mmvec_table, id.vars = "rn")
setnames(z_mmvec_table_long,
         c("rn","variable","value"),
         c("featureID","OTU_ID","z_mmvec"))

# raw mmvec table
mmvec_table_long <- melt.data.table(mmvec_table, id.vars = "featureID")
setnames(mmvec_table_long,
         c("variable","value"),
         c("OTU_ID", "raw_mmvec"))

# pearson correlations
cors_table <- as.data.table(all_cors_mat, keep.rownames = T)
setnames(cors_table, "rn", "featureID")
cors_table_long <- melt.data.table(cors_table, 
                                   id.vars = "featureID",
                                   variable.name = "OTU_ID",
                                   value.name = "spearman")

# put it all together so we can filter microbe/metabolite associtations by various criteria
pair_dat <- merge(z_mmvec_table_long, mmvec_table_long,
                  all.x = F,
                  all.y = F,
                  by = c("featureID","OTU_ID"))

pair_dat <- merge(pair_dat, cors_table_long,
                  all.x = F,
                  all.y = F,
                  by = c("featureID","OTU_ID"))

pair_dat <- merge(pair_dat, tip_data,
                  all.x = F,
                  all.y = F,
                  by = "featureID")

pair_dat <- merge(pair_dat, asv_table,
                  all.x = F,
                  all.y = F,
                  by = "OTU_ID")

## Filtering criteria:
# high raw mmvec probability
# high spearman
# library hit for the metabolite and well classified
names(pair_dat)

top_pairs <- cors_table_long[spearman > 0.9][1:20]

# write out plots for all selected pairs of metabolites and microbes
pdf("output/Correlations/top_pairs_linear_plots.pdf")
for (i in 1:nrow(top_pairs)) {
  graph_pair(top_pairs[i])
}
dev.off()


