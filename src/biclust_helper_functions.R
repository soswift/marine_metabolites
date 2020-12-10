# functions that help lot the correlation bicluster
library(data.table)
library(ComplexHeatmap)

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

# top_levels() makes updates a vector so that it only shows the top n categories
# useful for limiting colors in a plot
top_levels <- function(Vec, N = 10, exclude = "unclassified"){
  top_levels <-names(sort(table(Vec), decreasing = T)[1:N])
  Vec[!(Vec %in% top_levels) | Vec %in% exclude] <- "other"
  return(Vec)
}

# scale_dat() normalizes mmvec data by centering to mean for each metabolite
scale_dat <- function(mat_dat){
  mat_dat <- t(scale(t(mat_dat)))
  mat_dat[is.na(mat_dat)] <- 0
  return(mat_dat)
}

# get_top_type() returns the name of the highest sum for each row
# assumes the first column is an ID name, subsequent columns are sums by type
get_top_type <- function(type_sums){
  type_names <- colnames(type_sums[ , -1])
  top_types <-apply( type_sums[ , -1], 1, function(x) type_names[which.max(x)])
  names(top_types) <- type_sums[[1]]
  return(top_types)
}

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

# generate_phymap() makes a complex Heatmap that 'clusters' using phylogenetic trees (or equivalent qemistree)
# takes matrix where metabolites are rows, otus are columns, trees for both as dendrograms.
# If either tree is null, defaults to hclust.
# 'meta' files are used for annotation.
# returns a list where the first item is the heatmap plot, the second is the underlying matrix arranged by the dendrograms
generate_phymap <- function(
  # dendrograms
  micro_tree = fastree_raw,
  chem_tree = qemistree_raw,
  # correlation matrix
  correlation_mat = all_cors_mat,
  # annotation data
  micro_meta = asv_table,
  chem_meta = tip_data,
  # number of annotation values to display
  n_col = 20,
  # barchart sums
  micro_b = micro_sums,
  chem_b = chem_sums,
  # add boxplots
  box_plot = F
){
  
  # define sample type color scheme
  type_cols <- structure(c("#6469ed", "#e49c4c", "#7cc854","#808080"),
                         names = c("CCA","Coral","Limu","NS"))
  
  # the color scheme depends on the data type (categorical or numeric)
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
  
  # IMPORTANT: manually order correlation matrix to match dendrograms
  # when complexHeatmap does the clustering, it does this for you
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
  
  # relativize values for stacked barplots
  chem_b  <- t(apply(chem_b, 1,
                     FUN =  function(x) x/sum(x)))
  micro_b <- t(apply(micro_b, 1,
                     FUN =  function(x) x/sum(x)))
  
  # generate boxplot values
  if(isTRUE(box_plot)){
    # use abund_list() to generate boxplot data that matches heatmap matrix
    box_abunds <- abund_by_group(cor_mat = correlation_mat)
    
    # identify separate wrapper functions of box_plots() for metabolites and microbes
    # these will be called in columnAnnotation/rowAnnotation
    chem_boxplot <- function(index){
      box_plots(index = index, 
                type = "chem",
                abund_list = box_abunds,
                box_direction = "horizontal")
    }
    micro_boxplot <- function(index){
      box_plots(index = rev(index),
                type = "micro",
                abund_list = box_abunds,
                box_direction = "vertical")
    }
  }else{
    # if box_plot == F, don't generate box plots as row/column annotations, instead set to NULL
      micro_boxplot = NULL
      chem_boxplot = NULL
    }
  
  # metabolite row annotation
  ha_row <- rowAnnotation(
    Class = top_levels(chem_meta$class),
    Network = top_levels(chem_meta$componentindex,
                         exclude = "  -1", # drop single networks
                         N = 30),
    # stacked barplots
    Sample_Type = anno_barplot(
      apply(chem_b, 2, as.numeric),
      gp = gpar(fill = type_cols,
                col = type_cols,
                option = "A"),
      border = F,
      bar_width = 0.7
    ),
    # boxplot
    boxplot = chem_boxplot,
    annotation_width = c(1,1,5,20),
    width = unit(3, "in")
  )
  
  # microbe column annotation
  ha_col <- columnAnnotation(
    # boxplots
    boxplot = micro_boxplot,
    # stacked barplots
    Sample_type = anno_barplot(
      micro_b,
      gp = gpar(
        fill = type_cols,
        col = type_cols,
        option = "A"),
      border = F,
      bar_width = 0.7
    ),
    Order = top_levels(micro_meta$order, N = n_col),
    Class = top_levels(micro_meta$class, N = n_col),
    annotation_height = c(20,5,1,1),
    height = unit(3, "in")
  )
  
  # generate heatmap
  ht <- Heatmap(
    correlation_mat,
    cluster_rows = chem_dendro,
    cluster_columns = micro_dendro,
    
    # row parameters (metabolite)
    row_dend_width = unit(1, "in"),
    row_names_gp = gpar(fontsize = 6),
    show_row_names = F,
    
    # column parameters (microbe)
    column_dend_height = unit(1, "in"),
    column_dend_side = "bottom",
    column_names_gp = gpar(fontsize = 6),
    show_column_names = F,
    
    # colors
    na_col = "white",
    col = col_fun,
    rect_gp = gpar(col = "white"),
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

# sub_cor_mat() subsets corrlation matrix based on metadata groups
sub_cor_mat <- function(cor_mat = z_cat_mat,
                        chem_col = "class",
                        micro_col = "order",
                        chem_val = NULL,
                        micro_val = NULL) {
  
  if(!is.null(chem_val)){
  # subset by specified metabolite group
  cor_mat <- cor_mat[row.names(cor_mat) %in%
                     tip_data$featureID[tip_data[[chem_col]] %in% chem_val], ]
  }
  
  if(!is.null(micro_val)){
  # subset by specified microbe group
  cor_mat <- cor_mat[ , colnames(cor_mat) %in%
                        asv_table$OTU_ID[ asv_table[[micro_col]] %in% micro_val] ]
  }
  # drop empty rows/columns
  cor_mat <- cor_mat[apply(cor_mat, 1, function(x) !(all(is.na(x)))),
                      apply(cor_mat, 2, function(x) !(all(is.na(x))))]
  
  if (dim(cor_mat)[1] == 0) {
    stop("no rows, something is up with your metabolite subset")
  }
  if (dim(cor_mat)[2] == 0) {
    stop("no columns, something is up with your microbe subset")
  }
  return(cor_mat)
}

# abund_by_group() pulls out abundance matrices for metabolites/microbes for each sample group
# takes:
# 1) correlation matrix where microbes are columns and metabolites are rows
# 2) abundance tables for microbes/metabolites where samples are columns
# 3) metadata table with informaiton on sample groups
# id col matches 
abund_by_group <- function(cor_mat = z_cat_mat,
                           m_abund = micro_abund,
                           c_abund = chem_abund,
                           meta = sample_dat,
                           group = "sample_type",
                           sample_col = "sample_barcode"
) {
  # match abundance matrices to the correlation matrix
  m_abund <- m_abund[match(colnames(cor_mat), row.names(m_abund)), ]
  c_abund <- c_abund[match(row.names(cor_mat), row.names(c_abund)), ]
  
  # transform by log10 
  m_abund <- log10(1e-06 + m_abund)
  c_abund <- log10(1e-06 + c_abund)
  # vector of grouping ids
  group_vec <- meta[[group]]
  # for each group, get sample names
  group_samples <- lapply(unique(group_vec), function(x) {
    meta[[sample_col]][group_vec == x]
  })
  names(group_samples) <- unique(group_vec)
  # pull out microbe and metabolite abundances for each group of samples
  micro_list  <- lapply(group_samples, function(grp) {
      micro = m_abund[ , colnames(m_abund) %in% grp]
  }
  )
  chem_list <-  lapply(group_samples, function(grp){  
    chem = c_abund[ , colnames(c_abund) %in% grp]
  }
  )
  abund_list <- list(chem = chem_list,
                     micro = micro_list)
  
  # return list of microbe and metabolite abundance matrices
  return(abund_list)
}

# box_plots() generates metabolite and microbe box plots for annotating heatmap
# gets index passed from complexHeatmap (i.e. which row/column is being plotted)
# specify range using 'rg' so plots are legible
# specify type ("chem" or "micro")
# box_direction can be "vertical" (for columns) or "horizontal (for rows)
box_plots = function(index, type, abund_list = abund_list, box_direction = "horizontal") {
  # determine number of rows/columns
  nr = length(index)

  # subset data to type (chem or micro) and determine range of values
  type_abund <- abund_list[[type]]
  rg <<- range(abund_list)
  
  # define plot area for boxplots
  # note: if plotting for columns, the index order is reversed
  if(box_direction == "horizontal"){
  pushViewport(viewport(xscale = rg,
                        yscale = c(0.5, nr + 0.5)))
    index_vec <- seq_along(index)
    grid.xaxis()
  }
  if(box_direction == "vertical"){
  pushViewport(viewport(xscale = c(0.5, nr + 0.5),
                        yscale = rg))
    index_vec <- seq_along(index)
    grid.yaxis()
  }
  
  # plot the three boxplots, then shift to next row/colunn
  for(i in index_vec) {
    # CCA
    grid.boxplot(type_abund$CCA[index[i], ],
                 pos = nr-i+1 - 0.2,
                 box_width = 0.2, 
                 gp = gpar(fill = "blue"),
                 direction = box_direction,
                 pch = NA,
                 outline = F)
    # Coral
    grid.boxplot(type_abund$Coral[index[i], ],
                 pos = nr-i+1 + 0,
                 box_width = 0.2, 
                 gp = gpar(fill = "orange"),
                 direction = box_direction,
                 pch = NA,
                 outline = F)
    # Limu
    grid.boxplot(type_abund$Limu[index[i], ],
                 pos = nr-i+1 + 0.2,
                 box_width = 0.2, 
                 gp = gpar(fill = "green"),
                 direction = box_direction,
                 pch = NA,
                 outline = F)
  }
  popViewport()
}


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
