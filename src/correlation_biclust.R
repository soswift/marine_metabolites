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
# Qemistree output files include plotting information, all organized by unique 'id'
# Read tables of data related to metabolite qemistree 
class_tips   <- fread("data/raw/qemistree/labels.tsv", key = "id")
color_tips   <- fread("data/raw/qemistree/colors.tsv",  key = "id")
barplot_tips <- fread("data/raw/qemistree/barplots.tsv", key = "id")

# table that provides key to matching qemistree tip labels to other metabolomics data(mmvec etc.)
tip_to_feature <- fread("data/raw/qemistree/Fingerprints to features.tsv", key = "id")

# additional feature information (networks, classifications, etc.)
peak_data <- fread("data/processed/table_exports/all_marine_metabolite_tax_flat_table.csv")
setnames(peak_data,"V1", "featureID")
peak_data[ , featureID:= gsub("id_","",featureID)]

# Merge the qemistree outputs together with information we have on the individual metabolites
tip_data <- merge(class_tips, color_tips)
tip_data <- merge(tip_data, barplot_tips)
tip_data <- merge(tip_to_feature, tip_data)
tip_data <- merge(tip_data, peak_data, by = "featureID")

qem_id_map <- setNames(tip_data$featureID, tip_data$id)

# Read and Organize Metabolite Qemistree Tree ---------------------------------------------------
# read in qemistree tree file 
qemistree_raw <- read.tree(file = "data/raw/qemistree/qemistree.tree")

# update tip labels to 'featureID' to match the rest of our metabolite data
qemistree_raw$tip.label <- unname(qem_id_map[qemistree_raw$tip.label])
qem_tip_data <- tip_data[featureID %in% qemistree_raw$tip.label]

qemistree_raw <- drop.tip(qemistree_raw, tip_data[class == "unclassified", featureID])

# Read and Organize Microbe FastTree Phylogenetic Tree --------------------------------------------
# read in fastree file
fastree_raw <- read.tree("data/raw/FastTree_100.tre")

# Read in ASV and Metabolite Abundance Data-----------------------
# asv taxonomy data
asv_table <- fread("data/processed/table_exports/paired_marine_microbe_tax_flat_table.csv",
                   header = T)
setnames(asv_table, "V1", "OTU_ID")

# asv relative abundance
micro_abund <- read.csv("data/processed/table_exports/paired_marine_microbe_abundance_flat_table.csv",
                        row.names = 1)
micro_abund <- as.matrix(micro_abund)
colnames(micro_abund) <- sub("X","",colnames(micro_abund))

# sample data
sample_dat <- read.csv("data/processed/table_exports/paired_marine_microbe_sample_flat_table.csv")

# metabolite relative abundance
chem_abund <- read.csv("data/processed/table_exports/paired_marine_metabolite_abundance_flat_table.csv",
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
mmvec_table <- fread("data/raw/mmvec/Ranks_result.tsv", key = "featureid")

mmvec_table[ , featureid := sub("metabolite",
                                "",
                                featureid)]
setnames(mmvec_table,
         "featureid",
         "featureID")

## Read Spearman correlations from 'parallel_cor.R'
# all_cors <- readRDS("data/processed/all_cors_cutoff.rds")

# Filter Matrix For Heatmap-----------------------------------------------
## Subset correlation matrices to match qemistree and filter out uninteresting metabolites

# transform MMVEC scores to matrix
mmvec_mat <- as( mmvec_table[ , .SD, .SDcols = !"featureID"], "matrix")
row.names(mmvec_mat) <- mmvec_table$featureID

# Subset mmvec data to only include ids that are in qemistree and ASV metadata
# Also remove metabolites with a median mmvec score of <2
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

# Figure 5B networks
c(558,
  306,
  152,
  101,
  45,
  243,
  218,
  68) # ,acyl carnitine

# Just CCA networks
c(47, 359, 141, 48, 50, 222, 211)

# plot glycerophospholipds and   prenol lipids


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

classes <- c("Fatty Acyls",
"Carboxylic acids and derivatives",
"Prenol lipids",
"Steroids and steroid derivatives",
"Benzene and substituted derivatives",
"Organooxygen compounds"
)

for(cl in classes){
  cl_mat <- sub_cor_mat(cor_mat = z_cat_mat,
                        chem_col = "class",
                        micro_col = "order",
                        chem_val = cl,
                        micro_val = NULL)
  # generate heatmap
  ht_cl <- generate_phymap(correlation_mat = cl_mat,
                             box_plot = T)
  
  save_heatmap(ht_cl, paste0(cl, "_bicluster"))
  
}



i# Identifying Noteworthy Metabolites ------------------------------------------

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


