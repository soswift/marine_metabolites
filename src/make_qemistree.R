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

# get_sums() aggregates a matrix and generates sums by a metadata category
get_sums <- function(abundance, group_col, id_col, meta, new_name){
  abund_dt <- as.data.table(t(abundance), keep.rownames = T)
  abund_dt$group_vec <- meta[ match(abund_dt$rn, meta[[id_col]]), group_col ]
  abund_dt[, rn:=NULL]
  abund_dt <- melt.data.table(abund_dt,
                              id.vars = "group_vec", 
                              variable.name = "ID",
                              value.name = "abund")
  abund_dt[ , sum(abund), by = .(group_vec, ID)]
  abund_dt <- dcast(abund_dt, ID ~ group_vec,
                    value.var = "abund",
                    fun.aggregate = sum)
  setnames(abund_dt, "ID", new_name)
  return(abund_dt)
}

# get sample type sums for microbes and metabolites
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

# top_levels() makes updates a vector so that it only shows the top n categories
# useful for limiting colors in a plot
top_levels <- function(Vec, N = 10, exclude = "unclassified"){
  top_levels <-names(sort(table(Vec), decreasing = T)[1:N])
  Vec[!(Vec %in% top_levels) | Vec %in% exclude] <- "other"
  return(Vec)
}


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
all_cors <- readRDS("data/processed/all_cors_cutoff.rds")

# Filter Matrix For Heatmap-----------------------------------------------
## Subset correlation matrices to match qemistree and filter out uninteresting metabolites

# transform MMVEC scores to matrix
mmvec_mat <- as(mmvec_table[ , .SD, .SDcols = !"featureID"], "matrix")
row.names(mmvec_mat) <- mmvec_table$featureID

# Subset mmvec data to only include ids that are in qemistree and ASV metadata
# Also remove metabolites with a median mmvec score of <2
mmvec_mat <- mmvec_mat[row.names(mmvec_mat) %in% qemistree_raw$tip.label ,
                       colnames(mmvec_mat) %in% asv_table$OTU_ID]
mmvec_mat <- mmvec_mat[apply(mmvec_mat, 1, median) >= 2 , ]

## normalize mmvec data
scale_dat <- function(mat_dat){
  mat_dat <- t(scale(t(mat_dat)))
  mat_dat[is.na(mat_dat)] <- 0
  return(mat_dat)
}

z_mmvec_mat <- scale_dat(mmvec_mat)

# Spearman correlations: transform to matrix
all_cors_mat <- t(do.call(rbind, all_cors))
all_cors_mat[is.na(all_cors_mat)] <- 0
all_cors_mat[all_cors_mat < 0] <- 0
row.names(all_cors_mat) <- sub("id_", "", row.names(all_cors_mat))

# subset by qemistree tips and ASV metadata table
# Also remove microbes/metabolites with < 10 r values > 0.33 and any empty columns ( 0 correlation microbes)

all_cors_mat <- all_cors_mat[row.names(all_cors_mat) %in% qemistree_raw$tip.label,
                             colnames(all_cors_mat) %in% asv_table$OTU_ID]

# To compare spearman and mmvec, get the spearman scores for a matrix matching the mmvec matrix
match_cors <- all_cors_mat[match(row.names(z_mmvec_mat), row.names(all_cors_mat)), 
                           match(colnames(z_mmvec_mat), colnames(all_cors_mat))]
any(is.na(c(colnames(match_cors),row.names(match_cors))))

# independent from mmvec, cull down the correlation scores by cutoff criteria
all_cors_mat <- all_cors_mat[ apply(all_cors_mat, 1, 
                                    function(x) length(x[x > 0.33]) > 10) ,
                              apply(all_cors_mat, 2, function(x) length(x[x >0.33]) > 10)]


write.csv(all_cors_mat, "data/processed/spearman_correlations.csv")


# Generate Categorical Matrix -------------------------------------------------------------------

# get_top_type() returns the name of the highest sum for each row
# assumes the first column is an ID name, subsequent columns are sums by type
get_top_type <- function(type_sums){
 type_names <- colnames(type_sums[ , -1])
 top_types <-apply( type_sums[ , -1], 1, function(x) type_names[which.max(x)])
 names(top_types) <- type_sums[[1]]
 return(top_types)
}

micro_top_type <- get_top_type(micro_sums)
chem_top_type  <- get_top_type(chem_sums)

# type_match() takes the assigned types of microbes/metabolites and a matrix showing correlations between the two
# if types are the same, returns type name, if different returns "NS"
# results are filtered based on a correlation cutoff point (e.g. r = 0.4)
type_match <- function(mi_type, me_type, cor_mat, cutoff = 0.4){

  # full matrices of microbe x metabolite indicating types for each
  mim <- t(matrix(rep(mi_type, times =  length(me_type)),
                ncol = length(me_type),
                dimnames = list(names(mi_type), names(me_type))))
  
  mem <- matrix(rep(me_type, length(mi_type)),
                ncol = length(mi_type),
                dimnames = list(names(me_type), names(mi_type)))
  # check equal
  all(colnames(mim) == colnames(mem))
  
  # if types in the matrices don't match, replace with "NS" for 'Not Same'
  mimem <- mem
  mimem[mimem != mim] <- "NS"
  
  # match to correlation matrix
  mimem <- mimem[row.names(cor_mat), colnames(cor_mat)]
  
  # replace value with "NA" if below cutoff
  mimem[cor_mat < cutoff] <- NA
  return(mimem)
}

type_cat_mat <- type_match(mi_type = micro_top_type,
                           me_type = chem_top_type,
                           cor_mat = all_cors_mat,
                           cutoff = 0.4)

z_cat_mat <- type_match(mi_type = micro_top_type,
                        me_type = chem_top_type,
                        cor_mat = z_mmvec_mat,
                        cutoff = 1)

cat_match_mat <- type_match(mi_type = micro_top_type,
                            me_type = chem_top_type,
                            cor_mat = match_cors,
                            cutoff = 0.4)


# Generate Heatmaps ----------------------------------------------
# generate_phymap() makes a complex Heatmap that 'clusters' using phylogenetic trees (or equivalent qemistree)
# takes matrix where metabolites are rows, otus are columns, trees for both as dendrograms.
# If either tree is null, defaults to hclust.
# 'meta' files are used for annotation.
generate_phymap <- function(micro_tree = fastree_raw,
                            chem_tree = qemistree_raw,
                            correlation_mat = all_cors_mat,
                            micro_meta = asv_table,
                            chem_meta = tip_data,
                            micro_b = micro_sums,
                            chem_b = chem_sums,
                            n_col = 15){
  
  # define sample type color scheme
  type_cols <- structure(c("#6469ed", "#e49c4c", "#7cc854","#808080"),
                         names = c("CCA","Coral","Limu","NS"))
  
  # color scheme depends on data type
  
  if(class(correlation_mat[1]) == "numeric"){
  col_fun <- circlize::colorRamp2( c(0, 1),
                       c("white","black"))
  }
  if(class(correlation_mat[1]) == "character"){
  col_fun <- type_cols
  }
  
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
  
  # arrange metadata and barplots to match the correlation matrix and dendrograms
  chem_meta <- chem_meta[ match( labels(chem_dendro),
                                 chem_meta$featureID) , ]
  chem_b <- as.matrix(chem_b[match(labels(chem_dendro),
                               chem_b$featureID) , .(CCA,Coral,Limu)])
  
  micro_meta <- micro_meta[ match(labels(micro_dendro),
                                  micro_meta$OTU_ID) , ]
  micro_b <- as.matrix( micro_b[ match(labels(micro_dendro),
                                  micro_b$OTU_ID) , .(CCA,Coral,Limu) ])
  
  # relativize bars
  chem_b  <- t(apply(chem_b, 1,
                     FUN =  function(x) x/sum(x)))
  micro_b <- t(apply(micro_b, 1,
                     FUN =  function(x) x/sum(x)))
  
   # metabolite row annotation
  ha_row <- rowAnnotation(
   Class = top_levels(chem_meta$class),
   Network = top_levels(chem_meta$componentindex,
                        exclude = "  -1",
                        N = 30),
   # stacked barplots
   Sample_Type = anno_barplot(
     apply(chem_b, 2, as.numeric),
     gp = gpar(
       fill = type_cols,
       col = type_cols,
       option = "A"),
     border = F,
     bar_width = 0.7,
     width = unit(1, "in")
    )
   )
  
  # microbe column annotation
  ha_col <- columnAnnotation(
    Class = top_levels(micro_meta$class, N = n_col),
    Order = top_levels(micro_meta$order, N = n_col),
    # stacked barplots
    Sample_type = anno_barplot(
      micro_b,
      gp = gpar(
        fill = type_cols,
        col = type_cols,
        option = "A"),
      border = F,
      bar_width = 0.7,
      height = unit(1, "in")
      )
    )
  
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
    column_dend_height = unit(4, "in"),
    column_names_gp = gpar(fontsize = 6),
    show_column_names = F,
    
    # colors
    na_col = "white",
    col = col_fun,
    # annotations
    top_annotation =  ha_col,
    right_annotation = ha_row
  )
  
  print(paste0("Total tips in Qemistree:","  ",length(labels(chem_dendro))))
  print(paste0("Total tips in 16S phylogeny:", "  ", length(labels(micro_dendro))))
  
  return(list(heatmap = ht, matrix = correlation_mat))
}



# save_heatmap() saves a heatmap as png with specified filename
save_heatmap <- function(heatmap,
                         filename,
                         outdir = "output/heatmap/",
                         csv = F){
  # save image
  png(filename = paste0(outdir,filename,".png"),
      width = 30, height = 25, res = 300, units = "in" )
  print(heatmap[[1]])
  dev.off()
  if(isTRUE(csv)){
  # write data
  write.csv(heatmap[2], paste0(outdir, filename, ".csv"))
  }
}


## Generate heatmaps with and without qemistree
# Spearman
# ht_cor <- generate_phymap(correlation_mat = all_cors_mat)
# save_heatmap( ht_cor, "micro_meta_heatmap_spearman_qemistree")

# Categorical Spearman (i.e. top types)
ht_cat_cor <- generate_phymap(correlation_mat = type_cat_mat)
save_heatmap( ht_cat_cor, "sample_type_heatmap_spearman_qemistree")


# MMVEC
# ht_mmvec <- generate_phymap(correlation_mat = z_mmvec_mat)
# save_heatmap(ht_mmvec, "micro_meta_heatmap_mmvec_qemistree.png")

ht_cat_mmvec <- generate_phymap(correlation_mat = z_cat_mat)
save_heatmap(ht_cat_mmvec, "sample_type_heatmap_zmmvec_qemistree")

# Matching Spearman to MMVE
ht_cat_mmvec <- generate_phymap(correlation_mat = cat_match_mat)
save_heatmap(ht_cat_mmvec, "sample_type_heatmap_matching_spearman_qemistree")


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

# Graph of microbe abundance vs. metabolite abundance should look pretty
# graph_pair() takes a row from the paired data table and makes an x/y plot
graph_pair <-
  function(pairs_dat = top_pairs,
           chem_dat = chem_abund,
           micro_dat = micro_abund) {
    
    OTU = pairs_dat$OTU_ID[1]
    Feature = pairs_dat$featureID[1]
    
    x = chem_dat[row.names(chem_dat) == Feature ,]
    y = micro_dat[row.names(micro_dat) == OTU ,]
    
    print(plot(
      x = x,
      y = y,
      main = paste(OTU, "vs.", "Metabolite Feature", Feature),
      sub = paste("Taxonomic Family:", pairs_dat$family,
                  "Chemical Class:", pairs_dat$class.x, "Spearman:",
                  cor(x,y, method = "spearman")),
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




