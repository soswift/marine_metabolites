# This R script will generate the figures and analyses associated with the Waimea Valley marine microbial communities
# Began 08/29/2020
# Setup environment ----------------------------------------

# Load libraries

library(phyloseq)
library(data.table)
library(eulerr)
library(vegan)
library(pairwiseAdonis)
library(tidyr)
library(biomformat)

#library(speedyseq)

## identify data files

# microbial data files
abundance_file   <- "data/raw/abundance_table_100.shared"

sample_data_file <- "data/raw/new_marine_meta.csv"

taxonomy_file    <- "data/raw/annotations_100.taxonomy"

unifrac_file     <- "data/raw/unifrac_unweighted_100.csv"

tree_file        <- "data/raw/FastTree_100.tre"

# source lots of helper functions for plotting and subsetting
source("src/helper_functions.R")

# source function for reading and cleaning abundance table
source("src/clean_and_query_16S.R")

# source helper function for making phyloseq objects
source("src/make_phyloseq.R")

# source funciton for writing out phyloseq objects
source("src/physeq_csv_out.R")

# function for reading unifrac dists from flat tables
source("src/read_unifrac.R")

# function for formatting microbial data for mmvec
source("src/format_mmvec.R")


# set color scheme for sample types
sample_type_cols <- c(CCA = "#6469ed",
                      Coral = "#e49c4c",
                      Limu = "#7cc854" )

#setwd("/home/sean/Documents/Bioinformatics/Projects/SuperTransect/Analysis/marine_manuscript")
# 1. Load and clean microbial data ----------------------------
# Load and clean data from the metaflowimcs pipeline.
# Expects standards pipeline output format.
# Data is subset to only include desired marine samples.
# Read counts are subsampled to even depth and relative abundance transformed.

# read in data and sample data, 

marine_data <- clean_16S_tables(abundance_file = abundance_file,
                                 taxonomy_file = taxonomy_file,
                                 metadata_file = sample_data_file,
                                 description = "ST_marine",
                                 output_dir = "data/processed",
                                 id_column = "sequencing_id",
                                 cull = list(min.num = 3,
                                             min.abund = 0.00001,
                                             min.single.abund = 0.001))


marine_phy <- make_phyloseq(abund_file = "data/processed/ST_marine_abundance_table.csv",
                            meta_file  = "data/processed/ST_marine_metadata_table.csv",
                            tax_file   = "data/processed/ST_marine_taxonomy_table.csv",
                            id_column  = "sequencing_id")

rm(marine_data)

# check variation in sequencing depths and subsample to standard value
nrow(otu_table(marine_phy))
length(sample_names(marine_phy))
hist(sample_sums(marine_phy), breaks = 30)


marine_phy_rar <- rarefy_even_depth(marine_phy,
                                    rngseed = 1,
                                    sample.size = 15000,
                                    replace = F)
# check total OTUs
nrow(otu_table(marine_phy_rar))

# check which samples were dropped due to low read counts
marine_data[["metadata"]][!(sequencing_id %in% sample_names(marine_phy_rar))]
length(sample_names(marine_phy_rar))


# convert the OTU counts to relative abundance
final_marine_phy <- transform_sample_counts(marine_phy_rar,
                                            function(x) x / sum(x))

# drop sediment and water samples
final_marine_phy <- subset_samples(final_marine_phy,
                                   sample_type %in% c("Limu","CCA","Coral"))

# Re-classify genera based on our target genera for Limu and Coral
keep_genera <- c("Jania", "Halimeda","Galaxaura", "Porites", "Monitipora", "Pocillopora")

genus_data <- sample_data(final_marine_phy)$genus

sample_data(final_marine_phy)$genus <- ifelse(genus_data %in% keep_genera,
                                              genus_data,
                                              "Other" )

final_marine_phy <- prune_taxa(x = final_marine_phy,
                               taxa = taxa_sums(final_marine_phy) > 0)


saveRDS(final_marine_phy, file = "data/processed/final_marine_phy.rds")

rm(marine_phy_rar)

# write marine microbial sample tables out as flat files
physeq_csv_out(final_marine_phy,
               description = "all_marine_microbe",
               outdir = "data/processed/table_exports/")

# finally, get unifrac distance object for the selected samples
final_unifrac <- read_unifrac(unifrac_file = "data/raw/unifrac_weighted_100.csv", 
                              phyloseq_obj = final_marine_phy)

# write out unifrac subset to marine samples
saveRDS(final_unifrac, "data/processed/final_unifrac.rds")
write.csv(as.matrix(final_unifrac),
          "data/processed/table_exports/all_marine_microbe_unifrac_dist.csv")

# Final_unifrac contains all the relevant marine microbial samples.
# Some of these samples may not have corresponding metabolomics data.
# Incomparable samples will be dropped for joint analyses.


# 2. export data for MMVEC -----------------------------------------

# to match metabolomics data, we need to replace the id used for sequencing with sample names
# from the metadata, pull out sample barcodes for each sequencing id, matching the order in the abundance table

new_ids  <-  marine_data[["metadata"]][ match( colnames( get_taxa( final_marine_phy ) ), sequencing_id ),
             sample_barcode]

format_mmvec(phyloseq_obj = final_marine_phy,
             id_column = "Group",
             new_ids = new_ids,
             output_dir = "data/processed",
             description = "marine")

# 3. Load and Cleaning metabolomics data ------------------------------------------------

# Metabolomics data provided in one big table with peak identity and peak areas in the same file.
# The sample data file matches the microbial sample data, but has additional information about metabolomics sample ids(?)


# metabolomics data files
chem_raw_file    <- "data/raw/Helena_CCA_Coral_Limu_FeatureTable.txt"

chem_sample_file <- "data/raw/new_marine_meta.csv"


# this table combines abundance and metabolite feature metadata
chem_raw <- fread(chem_raw_file)

# read in metabolomics sample data
chem_meta <- fread(chem_sample_file)

# clean up sample genus info
keep_genera <- c("Jania", "Halimeda","Galaxaura", "Porites", "Monitipora", "Pocillopora")

chem_meta[ , genus := ifelse(genus %in% keep_genera, genus, "Other" )]


# format metabolite data using custom cleaning functions

## drop bad samples
bad_samples <- c(2638,
                 2839,
                 2684,
                 2750,
                 2815,
                 2650,
                 2795,
                 2862,
                 2905)

# call function to pull out abundance

chem_abund_w_blanks <- get_chem_abundance(chem_data = chem_raw,
                                          bad_samples = bad_samples)
# call function to get chem peak data (i.e. classifications)
chem_taxa_w_blanks <- get_chem_peak_data(chem_data = chem_raw)

# arrange metadata
chem_meta_w_blanks <- chem_meta[match(colnames(chem_abund_w_blanks), chem_meta$sample_barcode)]
row.names(chem_meta_w_blanks) <- chem_meta_w_blanks[ , sample_barcode]

chem_phy_w_blanks <- make_chem_phyloseq(chem_abund = chem_abund_w_blanks,
                                     chem_meta = chem_meta_w_blanks,
                                     chem_tax = chem_taxa_w_blanks,
                                     id_column = "sample_barcode")

# drop blanks and remove chemical features no longer present
not_blanks <- !grepl("blank",
               colnames(chem_abund_w_blanks))

chem_abund <- chem_abund_w_blanks[ , not_blanks]
chem_abund <- chem_abund[rowSums(chem_abund) > 0, ]
write.csv(chem_abund, "data/processed/table_exports/all_marine_metabolite_raw_abundance.csv")

# subset metadata to match cleaned abundance
chem_meta <- chem_meta[match(colnames(chem_abund), chem_meta$sample_barcode)]


# call function to make phyloseq of metabolite data
chem_phy <- make_chem_phyloseq(chem_abund = chem_abund,
                               chem_meta = chem_meta, 
                               chem_tax = chem_taxa_w_blanks,
                               id_column = "sample_barcode")

# convert to relative abundance
chem_phy <- transform_sample_counts(chem_phy,
                                    function(x) x / sum(x))


# write out as flat tables
physeq_csv_out(chem_phy,
               description = "all_marine_metabolite",
               outdir = "data/processed/table_exports")

physeq_csv_out(chem_phy_w_blanks,
               description = "all_marine_metabolite_w_blanks",
               outdir = "data/processed/table_exports/")

saveRDS(chem_phy, "data/processed/chem_phy.rds")

# clean up
rm(chem_abund_w_blanks, 
        chem_taxa_w_blanks,
        chem_meta_w_blanks,
        chem_abund,
        chem_meta)

# 4. All metabolite NMDS and calculate permanovas -----------------------------------------------------

# initalize chem distances list
chem_dist <- list()

# initialize plots list
chem_p <- list()

# NMDS for all samples, by type
chem_dist$canberra <- vegdist(veganifyOTU(chem_phy),
                              method = "canberra")
chem_dist$bray <- vegdist(veganifyOTU(chem_phy),
                          method = "bray")
chem_dist$euclidean <- vegdist(veganifyOTU(chem_phy),
                               method = "euclidean")

# run NMDS on all three distance objects

for(a in seq_along(chem_dist) ){
  
  dist_name <- names(chem_dist[a])
  
  chem_ord  <- ordinate(chem_phy,
                        method = "NMDS",
                        distance = chem_dist[[a]])
  
  chem_p[[a]] <- plot_ordination(
    chem_phy,
    chem_ord,
    color = "sample_type",
    shape = "site_name",
    title = paste0("All Metabolite NMDS, ", dist_name)
  ) + scale_color_manual(values = sample_type_cols)
  
}

g <-arrangeGrob(chem_p[[1]], chem_p[[2]], chem_p[[3]],
                nrow = 1, ncol = 3)

ggsave(paste0("output/NMDS/all_metabolites_NMDS_dists.pdf"),
       plot = g,
       width = 15,
       height = 5)

# run NMDS on metabolites with blanks

dist_blanks <- vegdist(veganifyOTU(chem_phy_w_blanks),
                       method = "bray")

chem_ord_w_blanks  <- ordinate(chem_phy_w_blanks,
                      method = "NMDS",
                      distance = dist_blanks)

blank_p <- plot_ordination(
  chem_phy_w_blanks,
  chem_ord_w_blanks,
  color = "sample_type",
  shape = "site_name",
  title = "All Metabolite NMDS, canberra"
) + scale_color_manual(values = c("darkgray",sample_type_cols))

ggsave("output/NMDS/all_metabolites_NMDS_w_blanks.pdf",
       plot = blank_p)

# chem permanovas 
# do permanova for all samples, by type and site

chem_anova_all <-
  do_permanova(
    chem_phy,
    var_name = "sample_type*site_name",
    dist_obj = chem_dist$bray,
    description = "All microbe samples"
  )
write.csv(chem_anova_all,
          "output/permanova/all_metabolite_permanova_by_site_type.csv")


# do pairwise permanova for all samples by type
chem_anova_pairwise <- pairwise.adonis2(chem_dist$bray ~ sample_type,
                                        data = as(sample_data(chem_phy),
                                                  "data.frame"))
write.csv(chem_anova_pairwise,
          "output/permanova/all_metabolite_pairwise_permanova_by_type.csv")

# calculate beta dispersion by group and test significance
chem_sam_dat <- as(sample_data(chem_phy), "data.frame")

chem_groups <- chem_sam_dat[match(labels(chem_dist$bray),
                                  chem_sam_dat$sample_barcode), 
                            "sample_type"]

chem_betadisp <- betadisper(chem_dist$bray, chem_groups)

# write out beta dispersion summary and tukey test resulst
sink("output/permanova/all_metabolite_betadispersion.txt" )
print(chem_betadisp)
print(TukeyHSD(chem_betadisp))
sink()

# write out distance matrices
write.csv(as.matrix(chem_dist$bray),
          "data/processed/table_exports/all_marine_metabolite_bray_dist.csv")
write.csv(as.matrix(dist_blanks),
          "data/processed/table_exports/all_marine_metabolite_w_blanks_bray_dist.csv")

# 5. Metabolite by sample type NMDS and permanova --------------------------------------------
sample_types <- unique(chem_phy@sam_data$sample_type)

chem_type_p <- list()

chem_anovas <- list()

for(i in sample_types){
  # subset
  a_phy  <- subset_samples(chem_phy, sample_type == i)
  a_phy  <- prune_taxa(taxa_sums(a_phy) > 0, a_phy)
  
  a_dist <- subset_dist(chem_dist[[1]], a_phy)
  
  # plot by site
  chem_type_p[[paste0(i,"_site")]] <- plot_NMDS(
    a_phy,
    "chem",
    color_var = "site_name",
    dist_method = "bray",
    dist_obj = a_dist
  )
  
  
  # plot by genus
  chem_type_p[[paste0(i,"_genus")]] <- plot_NMDS(
    a_phy,
    "chem",
    color_var = "genus",
    dist_method = "bray",
    dist_obj = a_dist
  )
  
  # calcluate permanovas for site
  chem_anovas[[i]]<- do_permanova(a_phy, 
                                  var_name = "site_name",
                                  dist_obj = a_dist, 
                                  description =  paste0(i, " metabolite samples"))
  
}

write.csv(bind_rows(chem_anovas),
          "output/permanova/sample_type_metabolite_permanova_results_by_site.csv")


g <-arrangeGrob(chem_type_p[[1]], chem_type_p[[3]], chem_type_p[[5]],
                chem_type_p[[2]], chem_type_p[[4]], chem_type_p[[6]],
                nrow = 2, ncol = 3)


ggsave(paste0("output/NMDS/sample_type_metabolites_NMDS.pdf"),
       plot = g, 
       width = 15,
       height = 5)



# Save Point: Load Cleaned Microbe and Metabolite Data----------------------------------------------------------------
final_marine_phy <- readRDS("data/processed/final_marine_phy.rds")
micro_phy <- final_marine_phy

final_unifrac <- readRDS("data/processed/final_unifrac.rds")

chem_phy <- readRDS("data/processed/chem_phy.rds")


# 6. All microbe NMDS and permanova ---------------------------------------------------------

# clean up genera names
sample_data(micro_phy)$genus <- ifelse(sample_data(micro_phy)$genus %in% keep_genera,
                                       sample_data(micro_phy)$genus, "Other")

# unifrac distances
micro_dist <- final_unifrac

# NMDS for all samples, by type
micro_ord  <- ordinate(micro_phy,
                       method = "NMDS",
                       distance = micro_dist)

plot_ordination(micro_phy, 
                micro_ord, 
                color = "sample_type",
                shape = "site_name",
                title = "All Microbe NMDS, Unifrac") +
  scale_color_manual(values = sample_type_cols)

ggsave("output/NMDS/all_microbes_NMDS.pdf",width = 7, height = 5)


# do permanova on all microbe sample type separation

micro_anova_all <- do_permanova(micro_phy,
                                  var_name = "sample_type*site_name",
                                  micro_dist,
                                  description = "all microbe samples")
write.csv(micro_anova_all, "output/permanova/all_microbe_permanova_by_site_type.csv")

micro_anova_pairwise <- pairwise.adonis2(micro_dist ~ sample_type,
                                         data = as(sample_data(micro_phy),
                                                   "data.frame"))
write.csv(micro_anova_pairwise, "output/permanova/all_microbe_pairwise_permanova_by_type")

# do betadispersion by group on all microbe samples
micro_sam_dat <- as(sample_data(micro_phy), "data.frame")

micro_groups <- micro_sam_dat[match(labels(micro_dist),
                                  micro_sam_dat$sequencing_id), 
                            "sample_type"]

micro_betadisp <- betadisper(micro_dist, micro_groups)

# write out beta dispersion summary and tukey test resulst
sink("output/permanova/all_microbe_betadispersion.txt" )
print(micro_betadisp)
print(TukeyHSD(micro_betadisp))
sink()

# 7. Microbe by sample type NMDS and permanovas ---------------------------

# initialize list of plots for microbial samples types
micro_type_p <- list()

micro_anovas <- list()

for(i in sample_types){
  
  # subset
  a_phy  <- subset_samples(micro_phy, sample_type == i)
  a_phy  <- prune_taxa(taxa_sums(a_phy) > 0, a_phy)
  
  # get distances
  a_dist <- subset_dist(micro_dist, a_phy)
  # plot by site
  micro_type_p[[paste0(i,"_site")]] <- plot_NMDS(a_phy, "microbe",
                        color_var = "site_name",
                        dist_method = "unifrac",
                        dist_obj = a_dist)
  
  # plot by genus
  micro_type_p[[paste0(i, "_genus")]] <- plot_NMDS(a_phy, "microbe",
                        color_var = "genus",
                        dist_method = "unifrac",
                        dist_obj = a_dist)
              
  # calcluate permanovas for site
  micro_anovas[[i]]<- do_permanova(a_phy, 
                                   var_name = "site_name",
                                   dist_obj = a_dist, 
                                   description =  paste0(i, " microbe samples"))
  
}

g <-arrangeGrob(micro_type_p[["Limu_site"]], micro_type_p[["Coral_site"]], micro_type_p[["CCA_site"]],
                micro_type_p[["Limu_genus"]], micro_type_p[["Coral_genus"]], micro_type_p[["CCA_genus"]],
                nrow = 2, ncol = 3)


ggsave("output/NMDS/sample_type_microbes_NMDS.pdf",
       plot = g, width = 15, height = 5)


write.csv(bind_rows(micro_anovas),
          "output/permanova/sample_type_microbe_permanova_by_site.csv")


# 8. Pair up microbes and metabolites ------------------------------------------

# find mismatches
micro_sample_bcodes <- micro_phy@sam_data$sample_barcode
print( paste("Total micro samples = ", length(micro_sample_bcodes)))

chem_sample_bcodes  <- chem_phy@sam_data$sample_barcode
print( paste("Total chem samples = ", length(chem_sample_bcodes)))

chem_no_micro <- chem_phy@sam_data[ ! (chem_phy@sam_data$sample_barcode %in% micro_sample_bcodes) ] %>%
  as.matrix() %>%
  as.data.frame()

print( paste("In chem data, but not in micro = ", nrow(chem_no_micro) ))

micro_no_chem <- micro_phy@sam_data[ ! (micro_phy@sam_data$sample_barcode %in% chem_sample_bcodes)]%>%
  as.matrix() %>%
  as.data.frame()

print( paste("In micro data, but not in chem = ", nrow(micro_no_chem)))

# subset both datasets
pair_chem_phy <- subset_samples(chem_phy,
                                sample_barcode %in% micro_phy@sam_data$sample_barcode)

pair_micro_phy  <- subset_samples(micro_phy,
                                  sample_barcode %in% chem_phy@sam_data$sample_barcode)

# check number samples in both data sets
nsamples(pair_chem_phy)
nsamples(pair_micro_phy)

# switch micro samples over to chem barcodes so everything is unified
## pull out otu table and sample data
otus    <- pair_micro_phy@otu_table
samples <- pair_micro_phy@sam_data
taxa    <- pair_micro_phy@tax_table

## change sample ids to ms barcodes
samples$sequence_id <- row.names(samples)
row.names(samples)  <- samples$sample_barcode
colnames(otus)      <- samples$sample_barcode

pair_micro_phy  <- phyloseq(samples, otus, taxa)

# change unifrac names
pair_micro_dist <- as.matrix(micro_dist)[samples$sequence_id, samples$sequence_id]
colnames(pair_micro_dist) <- samples$sample_barcode
row.names(pair_micro_dist) <- samples$sample_barcode

pair_micro_dist <- as.dist(pair_micro_dist)

rm(samples, otus, taxa)

# generate paired chem dist
pair_chem_dist <- vegdist(veganifyOTU(pair_chem_phy), method = "bray")

## check that names match
if (all( sample_names(pair_micro_phy) %in% sample_names(pair_chem_phy) ) ){
  message("chem and micro sample names match")
} else (
  warning("chem and micro sample names are different")
)

# write out paired tables
# microbes
physeq_csv_out(pair_micro_phy, 
               description = "paired_marine_microbe",
               outdir = "data/processed/table_exports")

write.csv(as.matrix(pair_micro_dist), 
          "data/processed/table_exports/paired_marine_microbe_unifrac_dist.csv")

saveRDS(pair_micro_phy, "data/processed/table_exports/paired_chem_phyloseq.rds")
write_biom(make_biom(data = otu_table(pair_micro_phy)), "data/processed/paired_microbe.biom")

# metabolites
physeq_csv_out(pair_chem_phy,
               description = "paired_marine_metabolite",
               outdir = "data/processed/table_exports")

write.csv(as.matrix(pair_chem_dist),
          "data/processed/table_exports/paired_marine_metabolite_bray_dist.csv")

saveRDS(pair_chem_phy, "data/processed/table_exports/paired_chem_phyloseq.rds")
write_biom(make_biom(data = otu_table(pair_chem_phy)), "data/processed/paired_metabolite.biom")


# 9. NMDS and Mantel of microbes and metabolites ----------------------------------------------

# Initialize plot list
p <- list()

# All samples
p$all <- paired_ordination(
  microbe_phy = pair_micro_phy,
  chem_phy = pair_chem_phy,
  description = "All",
  shape = "site_name",
  color = "sample_type"
)

sample_types <- unique(pair_chem_phy@sam_data$sample_type)

for(i in sample_types){
  # subset
  c_phy  <- subset_samples(pair_chem_phy, sample_type == i)
  c_phy  <- prune_taxa(taxa_sums(c_phy) > 0, c_phy)
  
  m_phy  <- subset_samples(pair_micro_phy, sample_type == i)
  m_phy   <- prune_taxa(taxa_sums(m_phy) > 0, m_phy)
  
  p[[i]] <- paired_ordination(
    microbe_phy = m_phy,
    chem_phy = c_phy,
    description = i,
    shape = "site_name",
    color = "genus"
  )
}


grid.arrange(p$all$micro      , p$all$chem,     p$all$proc,
             p$Limu$micro     , p$Limu$chem,    p$Limu$proc,
             p$Coral$micro    , p$Coral$chem,   p$Coral$proc,
             p$CCA$micro      , p$CCA$chem,     p$CCA$proc,
             nrow = 4, ncol = 3
)

g <-arrangeGrob(p$all$micro      , p$all$chem,     p$all$proc,
                p$Limu$micro     , p$Limu$chem,    p$Limu$proc,
                p$Coral$micro    , p$Coral$chem,   p$Coral$proc,
                p$CCA$micro      , p$CCA$chem,     p$CCA$proc,
                nrow = 4, ncol = 3
)

ggsave("output/NMDS/combined_microbe_metabolite_procrust.png", g, width = 15, height = 20)


# 10. Identify outlier CCA -----------------------------------------------------

cca_p <- list()

# Plot metabolite ordination with CCA labelled

chem_ord  <- ordinate(chem_phy,
                      method = "NMDS",
                      distance = chem_dist[[2]])

cca_p[["chem"]] <- plot_ordination(
  chem_phy,
  chem_ord,
  color = "sample_type",
  shape = "site_name",
  title = paste0("All Metabolite NMDS, ", "bray")
) + 
  scale_color_manual(values = sample_type_cols) + 
  geom_text(aes(label = ifelse(sample_type == "CCA",sample_barcode,"")))


# Plot microbe ordination with CCA labelled


micro_ord  <- ordinate(micro_phy,
                       method = "NMDS",
                       distance = micro_dist)

cca_p[["micro"]] <- plot_ordination(
  micro_phy,
  micro_ord,
  color = "sample_type",
  shape = "site_name",
  title = "All Microbe NMDS, Unifrac"
) +
  scale_color_manual(values = sample_type_cols) +
  geom_text(aes(label = ifelse(sample_type == "CCA", sample_barcode, "")))


g <-arrangeGrob(grobs = cca_p)

ggsave("output/NMDS/CCA_outlier_comparison.pdf", plot = g)

# plot paired ordination with CCA removed

no_cca_micro <- subset_samples(pair_micro_phy, sample_type != "CCA")
no_cca_chem <- subset_samples(pair_chem_phy, sample_type != "CCA")

p <- paired_ordination(
  microbe_phy = no_cca_micro,
  chem_phy = no_cca_chem,
  description = "All",
  shape = "site_name",
  color = "sample_type",
  unifrac_dist = pair_micro_dist
)

g <- arrangeGrob(p$micro, p$chem, p$proc, ncol = 3, nrow = 1)

ggsave("output/NMDS/no_cca_procrustes.pdf", plot = g, width = 15)
# 11. All microbe and metabolite heatmaps ------------------------------------------------------

plot_heatmap(phyloseq_obj = micro_phy, description = "microbe")

plot_heatmap(phyloseq_obj = chem_phy, description = "metabolite",
             dist_method = "canberra")

# 12. Metabolite and Microbe Eulerr diagrams -----------------------------------------------
# For corals and algae, generate euler plots by Genus

# merge all samples by sample type (coral, limu, cca)
all_chem_eul  <- euler_subset(chem_phy, group_by = "sample_type")

all_micro_eul <- euler_subset(micro_phy, group_by = "sample_type")

pdf("output/Euler/Metabolite_Sample_Type_Euler.pdf")
print(plot(all_chem_eul, fills = sample_type_cols, quantities = T))
dev.off()

pdf("output/Euler/Microbe_Sample_Type_Euler.pdf")
print(plot(all_micro_eul, fills = sample_type_cols, quantities = T))
dev.off()



# make a list of euler plots by host genus
# for each sample type (coral/limu) and data type (metabolite, microbe) 
eul_plots <- list()

# metabolite
eul_plots$chem_coral <- euler_subset(chem_phy, group_by = "genus",
                                     sample_type == "Coral")
eul_plots$chem_limu  <- euler_subset(chem_phy, group_by = "genus",
                                     sample_type == "Limu" & genus != "Other")
# microbe
eul_plots$micro_coral <- euler_subset(micro_phy, group_by = "genus",
                                        sample_type == "Coral")
eul_plots$micro_limu <- euler_subset(micro_phy, group_by = "genus",
                                        sample_type == "Limu" & genus != "Other")
# write out
for(p in names(eul_plots)){
  pdf(paste0("output/Euler/", p, "_genus_euler.pdf"))
  print(plot(eul_plots[[p]], quantities = T))
  dev.off()
}

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
# 
metabolite_abundance_file <- "data/processed/table_exports/all_marine_metabolite_abundance_flat_table.csv"
sample_data_file <- "data/processed/table_exports/all_marine_metabolite_sample_flat_table.csv"
metabolite_metadata_file <- "data/processed/table_exports//all_marine_metabolite_tax_flat_table.csv"

abund_raw <- read.csv(metabolite_abundance_file,
                      header = T,
                      row.names = 1)
sam_dat <- read.csv(sample_data_file,
                    header = T,
                    row.names = 1)

chem_dat <- read.csv(metabolite_metadata_file,
                     header = T,
                     row.names = 1)
chem_dat$featureID <- row.names(chem_dat)

# clean and arrange abundance data for VSURF random forest
abund<- as.data.frame(t(abund_raw))
row.names(sam_dat) <- paste0("X",row.names(sam_dat))
abund_clean <- abund[ row.names(sam_dat), ]
all(row.names(abund) == row.names(sam_dat))

# run VSURF on computing cluster
# see: https://journal.r-project.org/archive/2015/RJ-2015-018/RJ-2015-018.pdf

# vsurf.out <- VSURF(x = abund, y = sam_dat$sample_type,
#                    ncores = 8, parallel = T, clusterType = "FORK")

#saveRDS(vsurf.out, "out/sample_type_vsurf.rds")


# 15. LM on Random Forest Filtered Metabolits --------------------------------
# Run linear models on each metabolite that has high RF score

# Read in vsurf results (Variable Importance scores)
sample_type_vsurf <-readRDS("data/processed/sample_type_vsurf.rds")

sample_type_RF <- readRDS("data/processed/sample_type_rf.rds")

# pull out scores for each metabolite
all_metabolite_scores <-
  data.table(featureID = colnames(abund[, sample_type_vsurf$imp.mean.dec.ind]),
             RF_score = sample_type_vsurf$imp.mean.dec)

all_metabolite_scores <-  merge(all_metabolite_scores,
                                chem_dat, by = "featureID")

# identify scores higher than 1 Standard Deviation from mean
high_metabolite_scores <- 
              all_metabolite_scores[ RF_score > 2*sd(RF_score), ]

HS_feats <- high_metabolite_scores$featureID

# write out
write.csv(sample_type_RF$importance, "data/processed/all_metabolite_RF_importance.csv")
write.csv(sample_type_RF$localImportance, "data/processed/all_metabolite_RF_local_importance.csv")
fwrite(all_metabolite_scores, "data/processed/all_metabolite_vsurf_scores.csv")
fwrite(high_metabolite_scores, "data/processed/high_metabolite_vsurf_scores.csv")

# get abundance data and names for high scoring metabolites
HS_abund <- abund_clean[ , HS_feats]
# transform to normal-ish distribution
HS_abund <- as.data.frame(apply(MARGIN = 2, X= HS_abund, FUN =  function(x) asin(sqrt(x))))

HS_abund$sample_type <- as.character(sam_dat$sample_type)
HS_abund$site_name   <- as.character(sam_dat$site_name)

# run linear models on each feature with sample type as independent variable
run_lm <-function(feat_name, lm_abund){
       formula_call <- as.formula(paste0(feat_name, " ~ sample_type")) 
    
        lm_out <- lm(formula_call,
                         data = lm_abund)
        return(lm_out)
}

lm.out <- lapply(
                HS_feats,
                run_lm,
                lm_abund = HS_abund)

# summarize output of lm and anova for each of the selected metabolites

lm_results <- list()
for( i in 1:length(lm.out)){

  # summarize model output
  lm_sum <- summary(aov(lm.out[[i]]))

  # run tukey HSD
  feat_tuk_out <- TukeyHSD(aov(lm.out[[i]]))
  
  lm_results[[HS_feats[[i]]]]<-
  list(
    featureID = HS_feats[[i]],
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
}

# adjust p values
lm_results_df <-as.data.frame(do.call(rbind, lm_results))
lm_results_df$adj_p_val <- p.adjust(lm_results_df$p_val, method = "BH")
lm_results_df<- apply(lm_results_df, 2, unlist)

# write out
write.csv(lm_results_df, "data/processed/RF_filt_metabolites_sample_type_anovas.csv")
