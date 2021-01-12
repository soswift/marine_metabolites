## Load libraries --------------------------
# General utility
library(phyloseq)
library(data.table)
library(tidyr)
library(biomformat)
library(parallel)
# Plotting and analysis
library(vegan)
library(pairwiseAdonis)

# source lots of helper functions for plotting and subsetting
source("src/helper_functions.R")


# function for reading unifrac dists from flat tables
source("src/read_unifrac.R")


# set color scheme for sample types
sample_type_cols <- c(CCA = "#6469ed",
                      Coral = "#e49c4c",
                      Limu = "#7cc854" )
# set cores for parallel 
mc.cores = 3
# Load Cleaned Microbe and Metabolite Data----------------------------------------------------------------
final_marine_phy <- readRDS("data/processed/final_marine_phy.rds")
micro_phy <- final_marine_phy

final_unifrac <- readRDS("data/processed/final_unifrac.rds")

chem_phy <- readRDS("data/processed/chem_phy.rds")


# 13. Statistics on Metabolite Networks by Sample Type -----------------------------------------------
# For each metabolite network, run anova and look for fold changes in abundance between sample types

peak_data <- as.data.frame(
                      as(tax_table(chem_phy),
                         "matrix") )

# sum relative abundance of features in each network
# phyloseq treats taxonomy as hierarchical, so we need to drop most columns
# we only need 'componentindex', which is the network and 'cluster index' which is a unique id

tax_mat <- as(tax_table(chem_phy),
                  "matrix")
tax_mat <- tax_mat[ , c("componentindex","cluster index")]

net_tax_table <- tax_table(tax_mat)

net_phy <-phyloseq(sample_data(chem_phy),
                   otu_table(chem_phy),
                   net_tax_table)

network_merge <- tax_glom(net_phy,
                            "componentindex")

network_merge <- subset_taxa(network_merge,
                             componentindex != "  -1")

# the column componentindex identifies the network
unique_networks <- unique(as(tax_table(network_merge),
                             "matrix")[,"componentindex"]) 

# for each network, subset data and run anova 
net_data <- list()
for(a_network in unique_networks) {
  # subset compounds by network
  network_phy <- subset_taxa(network_merge,
                             componentindex == a_network)
  
  # create a map identifying samples by sample types
  sample_map <- setNames(sample_data(network_phy)[["sample_type"]],
                         sample_names(network_phy))
  sample_types <-c("Limu",
                   "CCA",
                   "Coral")
  
  # get mean relative abundance for each sample type present in the network
  mean_RAs <- list() 
  for(a_sample_type in sample_types){
  mean_RAs[[a_sample_type]] <- mean(
                                 otu_table(
                                 subset_samples(network_phy,
                                                sample_type == a_sample_type)))
  }
  
  # prior to modeling, normalize data by transform relative abundances
  arcsin_phy <- transform_sample_counts(network_phy,
                                          fun = function(x)
                                          asin(sqrt(x)) )
  # pull out transformed counts of sample_types, clean up
  sample_abunds <- data.frame(t( as(otu_table(arcsin_phy),
                        "matrix")))
  sample_abunds[is.na(sample_abunds)] <- 0
  
  # drop samples that don't contain this network
  #sample_abunds <-sample_abunds[rowSums(sample_abunds) > 0 , ]
  
  # if the network was only found in one sample, skip it
  #if(length(sample_abunds) == 1) next
  
  # make table with columns for sample type and abundance
  sample_table <- data.frame(sample_type = sample_map[row.names(sample_abunds)],
                             abundance = sample_abunds[[1]]
                             )
  # don't run anova if the network only occurs in one sample type
  # if(length( unique( sample_table$sample_type ) ) < 2){
  #   p_val = NA
  # } else {
  # run anova on sample type and get relevant outputs
  aov_out <-  aov(abundance ~ sample_type, data = sample_table)
  aov_sum <- summary(aov_out)
  
  p_val  <- aov_sum[[1]][["Pr(>F)"]][[1]]
  sum_sq <- aov_sum[[1]]["Sum Sq"][[1,1]]
  f_val <- aov_sum[[1]]["F value"][1,1]
  
  # run tukey HSD to test anova resulst
  tuk_out <-TukeyHSD(aov_out)
  
  # function to calculate fold change while handling zeroes
  fc <- function(x,y){
    if(y==0 && x == 0){
      return(0)
    }else if(y == 0){
      return(Inf)
    }else{
      return(log2(x/y))
    }
  }
  
  # assemble the output as a list
  # it shows for each network, mean RA by sample type, log2 fold changes, anova + tukey results
  net_data[[a_network]] <- list(
           network = a_network,
           limu_RA = mean_RAs$Limu,
           cca_RA = mean_RAs$CCA,
           coral_RA = mean_RAs$Coral,
           FC_LimuVCoral = fc(mean_RAs$Limu, mean_RAs$Coral),
           FC_LimuVCCA = fc(mean_RAs$Limu, mean_RAs$CCA),
           FC_CoralVCCA = fc(mean_RAs$Coral, mean_RAs$CCA),
           tuk_Coral_CCA_diff  = tuk_out[[1]]["Coral-CCA","diff"],
           tuk_Coral_CCA_p     = tuk_out[[1]]["Coral-CCA","p adj"],
           tuk_Limu_CCA_diff   = tuk_out[[1]]["Limu-CCA","diff"],
           tuk_Limu_CCA_p      = tuk_out[[1]]["Limu-CCA","p adj"],
           tuk_Limu_Coral_diff = tuk_out[[1]]["Limu-Coral","diff"],
           tuk_Limu_Coral_p    = tuk_out[[1]]["Limu-Coral","p adj"],
           f_stat = f_val,
           sum_of_sq = sum_sq,
           p_val = p_val
           )
  
  
# }
}

# put the fold changes for all networks into a nice data.frame
# note: this does not include networks that only showed up in one sample

network_fold_changes <- do.call(rbind, net_data)

network_fold_changes <-as.data.frame(apply(network_fold_changes, 2, as.numeric))
# adjust p values using 
network_fold_changes$adj_p_val <- p.adjust(network_fold_changes$p_val, method = "BH")

write.csv(network_fold_changes,
          "data/processed/network_anova_and_fold_changes.csv",
          row.names = F)


# 14. Variable Selection Using Random Forest ----------------------------
# Identify metabolite features that are highly correlated with sample types

# This section of code was run on the UH-HPC to speed things up

library(VSURF)

# set seed for parallel
# set.seed(2020, "L'Ecuyer-CMRG")
# 
# # read in metabolite abundance (relative abundance of peak areas) and sample data (information on samples)
# # abundance of metabolites will be used to predict sample type (Limu, Coral, CCA)


metabolite_abundance_file <- "data/processed/table_exports/all_marine_metabolite_abundance_flat_table.csv"
sample_data_file          <- "data/processed/table_exports/all_marine_metabolite_sample_flat_table.csv"
metabolite_metadata_file  <- "data/processed/table_exports//all_marine_metabolite_tax_flat_table.csv"

abund_raw <- read.csv(metabolite_abundance_file,
                      header = T,
                      row.names = 1)
sam_dat   <- read.csv(sample_data_file,
                    header = T,
                    row.names = 1)
sam_dat$sample_barcode <- paste0("X", sam_dat$sample_barcode)
chem_dat  <- read.csv(metabolite_metadata_file,
                     header = T,
                     row.names = 1)
chem_dat$featureID <- row.names(chem_dat)

# clean and arrange abundance data for VSURF random forest
abund <- as.data.frame(t(abund_raw))
row.names(sam_dat) <- paste0("X",row.names(sam_dat))
abund_clean <- abund[ row.names(sam_dat), ]
all(row.names(abund) == row.names(sam_dat))

# classify each metabolite to a sample type based on summed relative abundance

met_sums <- get_sums(t(abund_clean),
                        group_col = "sample_type",
                        id_col = "sample_barcode",
                        meta = sam_dat,
                        new_name = "featureID")
met_sample_class<- get_top_type(met_sums)

chem_dat$sample_class <- met_sample_class[chem_dat$featureID]

# file 'run_VSURF.R' run on computing cluster
# runs VSURF and random.Forest on cleaned metabolite data from above


# 15. LM on Random Forest Filtered Metabolits --------------------------------
# Run linear models on each metabolite that has high RF score

# Read in vsurf results (Variable Importance scores)
sample_type_vsurf <-readRDS("data/processed/sample_type_vsurf.rds")

sample_type_RF <- readRDS("data/processed/sample_type_rf.rds")



# pull out scores for each metabolite
all_metabolite_scores <-
  data.table(featureID = colnames(abund[, sample_type_vsurf$imp.mean.dec.ind]),
             VSURF_score = sample_type_vsurf$imp.mean.dec)
all_metabolite_scores[ , MeanDecreaseAccuracy := sample_type_RF$importance[featureID, "MeanDecreaseAccuracy" ]]

all_metabolite_scores <-  merge(all_metabolite_scores,
                                chem_dat, by = "featureID")

# identify scores higher than 1 Standard Deviation from mean
high_metabolite_scores <- 
              all_metabolite_scores[ MeanDecreaseAccuracy > mean(MeanDecreaseAccuracy) + 2*sd(MeanDecreaseAccuracy), ]

HS_feats <- high_metabolite_scores$featureID

# write out
write.csv(sample_type_RF$importance, "data/processed/all_metabolite_RF_importance.csv")
write.csv(sample_type_RF$localImportance, "data/processed/all_metabolite_RF_local_importance.csv")
fwrite(all_metabolite_scores, "data/processed/all_metabolite_vsurf_scores.csv")
fwrite(high_metabolite_scores, "data/processed/high_metabolite_vsurf_scores.csv")

# transform to normal-ish distribution
norm_abund <- as.data.frame(apply(MARGIN = 2, X= abund_clean, FUN =  function(x) asin(sqrt(x))))

norm_abund$sample_type <- as.character(sam_dat$sample_type)
norm_abund$site_name   <- as.character(sam_dat$site_name)
# get abundance data and names for high scoring metabolites
HS_abund <- norm_abund[ , colnames(norm_abund) %in% c(HS_feats, "sample_type", "site_name")]

# run linear models on each feature with sample type as independent variable
run_lm <-function(feat_name, lm_abund){
       formula_call <- as.formula(paste0(feat_name, " ~ sample_type")) 
    
        lm_out <- lm(formula_call,
                         data = lm_abund)
        return(lm_out)
}

# TODO add lines to calculate per sample type mean abundance

lm.HS <- lapply(
                HS_feats,
                run_lm,
                lm_abund = HS_abund)
names(lm.HS) <- HS_feats

# run linear models on all features
all_feats <- colnames(abund_clean)
lm.all <- lapply(
                all_feats,
                run_lm,
                lm_abund = norm_abund)
names(lm.all) <- all_feats

# summarize_lm() takes a feautre name and the output from run_lm() and generates summary stats for that feature
# when run with lapply, generates a list of stats for the provided features

summarize_lm <- function(i, lm.out){

  # summarize model output
  lm_sum <- summary(aov(lm.out[[i]]))

  # run tukey HSD
  feat_tuk_out <- TukeyHSD(aov(lm.out[[i]]))
  
  summary_out <-
  data.frame(
    featureID = i,
    p_val  = lm_sum[[1]][["Pr(>F)"]][[1]],
    sum_sq = lm_sum[[1]]["Sum Sq"][[1,1]],
    f_val  = lm_sum[[1]]["F value"][1,1],
    tuk_Coral_CCA_diff  = feat_tuk_out[[1]]["Coral-CCA","diff"],
    tuk_Coral_CCA_p     = feat_tuk_out[[1]]["Coral-CCA","p adj"],
    tuk_Limu_CCA_diff   = feat_tuk_out[[1]]["Limu-CCA","diff"],
    tuk_Limu_CCA_p      = feat_tuk_out[[1]]["Limu-CCA","p adj"],
    tuk_Limu_Coral_diff = feat_tuk_out[[1]]["Limu-Coral","diff"],
    tuk_Limu_Coral_p    = feat_tuk_out[[1]]["Limu-Coral","p adj"]
  )
  return(summary_out)
}


# get results for features that had a high random forest score
# organize results in a data.frame
HS_results <- lapply(HS_feats, summarize_lm, lm.out = lm.HS)
HS_results <- do.call("rbind", HS_results)

HS_results <- merge(HS_results,
                   all_metabolite_scores,
                   all.x = T, by = "featureID")


# get results for all features
all_results <- lapply(all_feats, summarize_lm, lm.out = lm.all)
all_results <- do.call("rbind", all_results)

all_results <- merge(all_results,
                    all_metabolite_scores,
                    all.x = T, by = "featureID")

# adjust p values when considering all features
all_results$adj_p_val <- p.adjust(all_results$p_val, method = "BH")

pdf(file = "output/randomforest/pvalue_vs_RF.pdf")
plot(log10(all_results$adj_p_val), all_results$VSURF_score,
     xlab =  "LM log10 P-value",
     ylab = "VSURF Importance Score")

plot(log10(all_results$adj_p_val), all_results$MDA,
     xlab =  "LM log10 P-value",
     ylab = "RandomForest MDA Score")
dev.off()

# write out results
write.csv(HS_results, "data/processed/RF_filt_metabolites_sample_type_anovas.csv", row.names = F)
write.csv(all_results, "data/processed/all_metabolites_sample_type_anovas.csv", row.names = F)


## Compare randomForest selected metabolites with MMVEC ordination --------------------------------------
# RF identified metabolites that are important for classifying sample types
# Do these metabolites also significantly separate out in the MMVEC ordination?
source("src/biclust_helper_functions.R")
# read in
# RF importance scores and linear model summary informations
HS_results  <- fread("data/processed/RF_filt_metabolites_sample_type_anovas.csv")
all_results <- fread("data/processed/all_metabolites_sample_type_anovas.csv")
all_results

# MMVEC ordination eigenvalues for microbes and metabolites
mic_mmvec_ord <- fread("data/raw/mmvec/microbe_ordination.tsv") 
met_mmvec_ord <- fread("data/raw/mmvec/metabolite_ordination.tsv") 

# clean up microbe and feature names
met_mmvec_ord[ ,featureID:= sub("metabolite", "id_", featureID) ]

# check that features selected by RF are in the ordination
all(HS_results$featureID %in% met_mmvec_ord$featureID)

# merge ordination scores with RF/linear model information
met_rf_ord <- merge(met_mmvec_ord, HS_results, by = "featureID", all.x = T)
met_lm_ord <- merge(met_mmvec_ord, all_results, by = "featureID", all.x = T)

# add columns that distinguish alpha and color
met_rf_ord$RF_selection <- ifelse(is.na(met_rf_ord$p_val), "Not selected", "Selected by RF")
met_rf_ord$transp <- ifelse(is.na(met_rf_ord$p_val), 0.2, 0.8)

met_lm_ord$LM_selection <- ifelse(met_lm_ord$adj_p_val < 0.05 & !is.na(met_lm_ord$adj_p_val), "Not significant", "Significant")
met_lm_ord$transp <- ifelse(met_lm_ord$LM_selection == "Not significant", 0.2 , 0.8)

# plot the ordination with high RF scored metabolites
p<- ggplot(data = met_rf_ord,
           aes(x = X, y = Y, shape = RF_selection, col = sample_class))+
          geom_point(alpha = met_rf_ord$transp) +
          scale_shape_manual(values = c(1,19))+
          theme_minimal()
    
p
ggsave("output/randomforest/RF_variables_mmvec_biplot.pdf", 
       plot = p, width = 9, height = 8)

# plot the ordination with low FDR pvals
p2 <- ggplot(data = met_lm_ord,
           aes(x = X, y = Y, shape = LM_selection, col = sample_class))+
          geom_point(alpha = met_lm_ord$transp) +
          scale_shape_manual(values = c(1,19))+
          theme_minimal()
    
p2
ggsave("output/randomforest/LM_variables_mmvec_biplot.pdf", 
       plot = p2, width = 9, height = 8)
