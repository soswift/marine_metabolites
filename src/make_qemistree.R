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

detach("package:speedyseq")
unloadNamespace("speedyseq")
detach("package:phyloseq")
unloadNamespace("phyloseq")

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

# make new column that only shows top 5 classes for plot
top_classes <- tip_data[class != "unclassified" , .N, by = class][order(-N)][1:10, class]

tip_data$plot_class <- tip_data$class
tip_data$plot_class[!(tip_data$heatmap_class %in% top_classes)] <- "other"

# Read MMVEC and correlation tables -------------------------------

# MMVEC data (matrix showing pairwise mmvec scores)
mmvec_table <- fread("data/raw/mmvec/Ranks_result.tsv", key = "featureid")

mmvec_table[ , featureid := sub("metabolite",
                                "",
                                featureid)]
setnames(mmvec_table,
         "featureid",
         "featureID")

# Pearson correlation data (long format showing pairwise comparisons)
correl_raw <-read.csv("data/processed/2347nodes.780kedges.pvals.csv")

# clean
correl_raw$featureID  <- gsub( ".+d-(\\d+).*", "\\1", 
                               correl_raw$Label,
                               perl = T)
correl_raw$OTU_ID <- gsub( ".*(Otu\\d+).+","\\1", 
                           correl_raw$Label,
                           perl = T)
correl_raw$OTU_ID <- gsub( "Otu(\\d{4}$)","Otu0\\1", 
                           cors$OTU_ID)

# subset columns and coerce into wide matrix
cors <- correl_raw[ , c("featureID",
                        "OTU_ID", 
                        "pearson")]

cors_wide <- reshape(cors,
                     idvar = "featureID",
                     timevar = "OTU_ID",
                     direction = "wide")

row.names(cors_wide) <- cors_wide$featureID
colnames(cors_wide)  <- sub("pearson.","",colnames(cors_wide))


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

# metabolite relative abundance
chem_abund <- read.csv("data/processed/table_exports/paired_marine_metabolite_abundance_flat_table.csv",
                       row.names = 1)
chem_abund  <- as.matrix(chem_abund)
micro_abund <- micro_abund[ , colnames(chem_abund)]

# check samples match
colnames(micro_abund) == colnames(chem_abund)

# Filter Matrix For Heatmap-----------------------------------------------

# Correlation matrix
cors_mat <- as(cors_wide[ , -1 ], "matrix")

# MMVEC matrix
mmvec_mat <- as(mmvec_table[ , .SD, .SDcols = !"featureID"], "matrix")
row.names(mmvec_mat) <- mmvec_table$featureID

## Subset correlation matrices to match qemistree and filter out noise
# Subset mmvec data to only include ids that are in qemistree and ASV metadata
# Also remove metabolites with a median mmvec score of <2
mmvec_mat <- mmvec_mat[row.names(mmvec_mat) %in% qemistree_raw$tip.label , 
                       colnames(mmvec_mat) %in% asv_table$OTU_ID]
mmvec_mat <- mmvec_mat[apply(mmvec_mat, 1, median) >= 2 , ]

# For spearman data, subset by qemistree and ASVs
# Also remove metabolites with < 10 r values > 0.3
cors_mat <- cors_mat[row.names(cors_mat) %in% qemistree_raw$tip.label , 
                     colnames(cors_mat) %in% asv_table$OTU_ID]
cors_mat <- cors_mat[apply(cors_mat, 1, function(x) length(x > 0.3) > 10 ) , ]

## normalize mmvec data
scale_dat <- function(mat_dat){
  mat_dat <- t(scale(t(mat_dat)))
  mat_dat[is.na(mat_dat)] <- 0
  return(mat_dat)
}

z_mmvec_mat <- scale_dat(mmvec_mat)

# Generate Heatmaps ----------------------------------------------
# generate_phymap() makes a complex Heatmap that 'clusters' using phylogenetic trees (or equivalent qemistree)
# takes matrix where metabolites are rows, otus are columns, trees for both as dendrograms.
# If either tree is null, defaults to hclust.
# 'meta' files are used for annotation.
generate_phymap <- function(micro_tree = fastree_raw,
                            chem_tree = qemistree_raw,
                            correlation_mat = cors_mat,
                            micro_meta = asv_table,
                            chem_meta = tip_data,
                            chem_anno = "plot_class",
                            micro_anno = "class"){
  
  # color scheme
  col_fun <- circlize::colorRamp2( c(0, 2),
                       c("white","black"))
  
  ## Microbe Tree
  # if tree is provided
  # subset microbe tree to ASVs that are present in correlation matrix and vice versa
  if(!is.null(micro_tree)){
    mat_asvs         <- colnames(correlation_mat)
    micro_tree_asvs  <- mat_asvs[ mat_asvs %in% micro_tree$tip.label]
    micro_tree_clean <- keep.tip(micro_tree, micro_tree_asvs)
    correlation_mat  <- correlation_mat[ , micro_tree_asvs]
    
    # microbe tree as dendrogram
    micro_dendro <- as.dendrogram.phylo(micro_tree_clean)
    
  # if tree is null, then use hclust to generate dendrogram
  } else {
    dist_mat <- t(correlation_mat)
    dist_mat[is.na(dist_mat)] <- 0
    micro_dendro <- as.dendrogram( hclust ( dist( dist_mat ) ) )
  }
  
  ## Metabolite Tree
  # if tree is provided
  # subset metabolite tree to features in the correlation matrix and vice versa
  if(!is.null(chem_tree)){
    mat_feats       <- row.names(correlation_mat)
    chem_tree_feats <-
      mat_feats[mat_feats %in% chem_tree$tip.label]
    
    chem_tree_clean <- keep.tip(chem_tree, chem_tree_feats)
    correlation_mat <- correlation_mat[chem_tree_feats , ]
    
    # metabolite tree as dendrogram
    chem_dendro  <- as.dendrogram.phylo(chem_tree_clean)
    
    # if tree is null, use hclust to generate dendrogram
  } else {
    dist_mat <- correlation_mat
    dist_mat[is.na(dist_mat)] <- 0
    chem_dendro <- as.dendrogram( hclust ( dist( dist_mat ) ) )
  }
  
  # IMPORTANT: order correlation matrix to match dendrograms
  correlation_mat <- correlation_mat[ labels(chem_dendro),
                                      labels(micro_dendro)]
  
  # arrange metadata to match correlation matrix
  chem_meta <- chem_meta[ match( labels(chem_dendro),
                                 chem_meta$featureID) , ]
  
  micro_meta <- micro_meta[ match(labels(micro_dendro),
                                  micro_meta$OTU_ID) , ]
  
  
  # metabolite row annotation
  ha_row <- rowAnnotation(
    class = chem_meta[[chem_anno]])
  
  # microbe column annotation
  ha_col <- columnAnnotation(
    class = micro_meta[[micro_anno]])
  
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
    na_col = "white",
    col = col_fun,
    # annotations
    top_annotation =  ha_col,
    left_annotation = ha_row
  )
  
  print(paste0("Total tips in Qemistree:","  ",Ntip(chem_tree_clean)))
  print(paste0("Total tips in 16S phylogeny:", "  ", Ntip(micro_tree_clean)))
  
  return(ht)
}

# save_heatmap() saves a heatmap as png with specified filename
save_heatmap <- function(heatmap, filename, outdir = "output/heatmap/"){
  png(filename = paste0(outdir,filename),
      width = 30, height = 30, res = 300, units = "in" )
  print(heatmap)
  dev.off()
}

## Generate heatmaps with and without qemistree
# Pearson
ht_cor <- generate_phymap(correlation_mat = cors_mat)
save_heatmap( ht_cor, "micro_meta_heatmap_pearson_qemistree.png")

ht_cor_clust <- generate_phymap(correlation_mat = cors_mat,
                                chem_tree = NULL)
save_heatmap( ht_cor_clust, "micro_meta_heatmap_pearson_hclust.png")

# MMVEC
ht_mmvec <- generate_phymap(correlation_mat = z_mmvec_mat)
save_heatmap(ht_mmvec, "micro_meta_heatmap_mmvec_qemistree.png")

ht_mmvec_clust <- generate_phymap(correlation_mat = z_mmvec_mat, 
                                  chem_tree = NULL)
save_heatmap(ht_mmvec_clust, "micro_meta_heatmap_mmvec_hclust.png")

ht_raw_mmvec <- generate_phymap(correlation_mat = mmvec_mat)
save_heatmap(ht_raw_mmvec, "raw_mmvec_qemistree.png")

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
cors_table_long <- as.data.table(correl_raw)

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
# high pearson rm significant pearson p, and high n
# library hit for the metabolite and well classified
names(pair_dat)

top_pairs <- pair_dat[ 
  abs(raw_mmvec) > sd(raw_mmvec) + mean(raw_mmvec) 
  & abs(z_mmvec) > 0.5 
  & abs(pearson) > 0.7
  & pval < 0.05 
  & ms2_library_match != "missing"
  ]

# Graph of microbe abundance vs. metabolite abundance should look pretty
# graph_pair() takes a row from the paired data table and makes an x/y plot
graph_pair <-
  function(pairs_dat = top_pairs,
           chem_dat = chem_abund,
           micro_dat = micro_abund) {
    
    OTU = pairs_dat$OTU_ID
    Feature = pairs_dat$featureID
    
    print(plot(
      x = log(chem_dat[paste0("id_", Feature),]),
      y = log(micro_dat[OTU ,]),
      main = paste(OTU, "vs.", "Metabolite Feature", Feature),
      sub = paste("Taxonomic Family:", pairs_dat$family,
                  "Chemical Class:", pairs_dat$class.x),
      xlab = "Log Metabolite RA",
      ylab = "Log OTU RA"
    ))
  }

# write out plots for all selected pairs of metabolites and microbes
pdf("output/Correlations/top_pairs_linear_plots.pdf")
for (i in 1:nrow(top_pairs)) {
  
  graph_pair(top_pairs[i])
  
}
dev.off()




