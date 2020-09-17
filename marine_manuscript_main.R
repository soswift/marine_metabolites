# This R script will generate the figures and analyses associated with the Waimea Valley marine microbial communities
# Began 08/29/2020


# Load libraries

library(phyloseq)
library(data.table)
library(dplyr)
library(eulerr)


setwd("/home/sean/Documents/Bioinformatics/Projects/SuperTransect/Analysis/marine_manuscript")
# 1. Loading and clean microbial data ----------------------------

## identify data files

# microbial data files
abundance_file   <- "data/raw/abundance_table_100.shared"

sample_data_file <- "data/raw/new_marine_meta.csv"

taxonomy_file    <- "data/raw/annotations_100.taxonomy"

unifrac_file     <- "data/raw/unifrac_unweighted_100.csv"

tree_file        <- "data/raw/FastTree_100.tre"




# source function for reading and cleaning abundance table
source("src/clean_and_query_16S.R")

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
source("src/make_phyloseq.R")

marine_phy <- make_phyloseq(abund_file = "data/processed/ST_marine_abundance_table.csv",
                            meta_file  = "data/processed/ST_marine_metadata_table.csv",
                            tax_file   = "data/processed/ST_marine_taxonomy_table.csv",
                            id_column  = "sequencing_id")

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
final_marine_phy <- transform_sample_counts(marine_phy_rar, function(x) x / sum(x))

# drop sediment and water samples
final_marine_phy <- subset_samples(final_marine_phy, sample_type %in% c("Limu","CCA","Coral"))
final_marine_phy <- prune_taxa(taxa_sums(final_marine_phy) > 0, final_marine_phy)


saveRDS(final_marine_phy, file = "data/processed/final_marine_phy.rds")



# read in unifrac distances and match to samples in processed phyloseq object

read_unifrac <- function(unifrac_file, phyloseq_obj){

unifrac <- read.csv(unifrac_file, header = T, row.names = 1)

colnames(unifrac) <- sub("X","",colnames(unifrac))

unifrac <- unifrac[sample_names(phyloseq_obj), sample_names(phyloseq_obj)]

unifrac <- as.dist(as.matrix(unifrac))

print(unifrac)


}

final_unifrac <- read_unifrac(unifrac_file = "data/raw/unifrac_weighted_100.csv", 
                              phyloseq_obj = final_marine_phy)

saveRDS(final_unifrac, "data/processed/final_unifrac.rds")

# 2. export data for MMVEC -----------------------------------------
source("src/format_mmvec.R")

# to match metabolomics data, we need to replace the id used for sequencing with sample names
# from the metadata, pull out sample barcodes for each sequencing id, matching the order in the abundance table

new_ids  <-  marine_data[["metadata"]][ match( colnames( get_taxa( final_marine_phy ) ), sequencing_id ),
             sample_barcode]

format_mmvec(phyloseq_obj = final_marine_phy,
             id_column = "Group",
             new_ids = new_ids,
             output_dir = "data/processed",
             description = "marine")

# 3. Loading and Cleaning metabolomics data ------------------------------------------------


# metabolomics data files
chem_raw_file <- "data/raw/Helena_CCA_Coral_Limu_FeatureTable.txt"

chem_sample_file <- "data/raw/new_marine_meta.csv"

source("src/make_pcoa.R")

# this table combines abundance and metabolite feature metadata
chem_raw <- fread(chem_raw_file)

# read in metabolomics sample data
chem_meta <- fread(chem_sample_file)

# clean up sample genus info
keep_genera <- c("Jania", "Halimeda","Dichotoma", "Porites", "Monitipora", "Pocillopora")

chem_meta[ , genus := ifelse(genus %in% keep_genera, genus, "Other" )]


# format metabolite data using tailored cleaning functions

# call function to pull out abundance
chem_abund_w_blanks <- get_chem_abundance(chem_data = chem_raw)

chem_meta_w_blanks <- chem_meta[match(colnames(chem_abund_w_blanks), chem_meta$sample_barcode)]

chem_phy_w_blanks <- make_chem_phyloseq(chem_abund = chem_abund_w_blanks,
                   chem_meta = chem_meta_w_blanks,
                   id_column = "sample_barcode")

# drop blanks
blanks <- grep("blank", colnames(chem_abund_w_blanks), value = T)

chem_abund <- chem_abund_w_blanks[ , c(blanks):=NULL ]

# subset metadata to match cleaned abundance
chem_meta <- chem_meta[match(colnames(chem_abund), chem_meta$sample_barcode)]




# call function to make phyloseq of metabolite data
chem_phy <- make_chem_phyloseq(chem_abund = chem_abund,
                               chem_meta = chem_meta, 
                               id_column = "sample_barcode")

# convert to relative abundance
chem_phy <- transform_sample_counts(chem_phy , function(x) x / sum(x))


# set color scheme for sample types
sample_type_cols <- c(CCA = "#6469ed", Coral = "#e49c4c", Limu = "#7cc854" )

# 4. Plot metabolite NMDS and calculate permanovas -----------------------------------------------------

# initalize chem distances list
chem_dist <- list()

# initialize plots list
chem_p <- list()

# NMDS for all samples, by type
chem_dist$canberra <- vegdist(veganifyOTU(chem_phy), method = "canberra")
chem_dist$bray <- vegdist(veganifyOTU(chem_phy), method = "bray")
chem_dist$euclidean <- vegdist(veganifyOTU(chem_phy), method = "euclidean")

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
  ) + scale_color_manual(values = sample_type_cols)+ 
    geom_text()
  
}

g <-arrangeGrob(chem_p[[1]], chem_p[[2]], chem_p[[3]],
                nrow = 1, ncol = 3)

ggsave(paste0("output/NMDS/all_metabolites_NMDS_dists.pdf"),  plot = g, width = 15, height = 5)


# run NMDS on metabolites with blanks

dist_blanks <- vegdist(veganifyOTU(chem_phy_w_blanks), method = "canberra")

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

ggsave("output/NMDS/all_metabolites_NMDS_w_blanks.pdf", plot = blank_p)

# chem permanovas 
# permanova for all samples, by type and site

# initalize anova table
chem_anovas <- list()

chem_anovas$all <-
  do_permanova(
    chem_phy,
    var_name = "sample_type*site_name",
    dist_obj = chem_dist[[1]],
    description = "All microbe samples"
  )

# 5. Metabolite by sample type NMDS and permanova --------------------------------------------
sample_types <- unique(chem_phy@sam_data$sample_type)

chem_type_p <- list()

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

write.csv(bind_rows(chem_anovas), "output/NMDS/metabolite_anova_results.csv")


g <-arrangeGrob(chem_type_p[[1]], chem_type_p[[3]], chem_type_p[[5]],
                chem_type_p[[2]], chem_type_p[[4]], chem_type_p[[6]],
                nrow = 2, ncol = 3)


ggsave(paste0("output/NMDS/sample_type_metabolites_NMDS.pdf"),  plot = g, width = 15, height = 5)



# Load Saved  Microbial Data----------------------------------------------------------------
final_marine_phy <- readRDS("data/processed/final_marine_phy.rds")

final_unifrac <- readRDS("data/processed/final_unifrac.rds")



# 6. All microbe NMDS and permanova ---------------------------------------------------------

# phyloseq object
micro_phy <- final_marine_phy

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


# permanova on all microbe sample type separation
micro_anovas <- list()
micro_anovas$all <- do_permanova(micro_phy, var_name = "sample_type", micro_dist, description = "all microbe samples")


# 7. Microbe by sample type NMDS and permanovas ---------------------------

# initialize list of plots for microbial samples types
micro_type_p <- list()

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


ggsave(paste0("output/NMDS/sample_type_microbes_NMDS.pdf"),  plot = g, width = 15, height = 5)


write.csv(bind_rows(micro_anovas), "output/NMDS/microbe_anova_results.csv")


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
pair_chem_phy <- subset_samples(chem_phy, sample_barcode %in% micro_phy@sam_data$sample_barcode)

pair_micro_phy  <- subset_samples(micro_phy, sample_barcode %in% chem_phy@sam_data$sample_barcode)

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

## check that names match
if (all( sample_names(pair_micro_phy) %in% sample_names(pair_chem_phy) ) ){
  message("chem and micro sample names match")
} else (
  warning("chem and micro sample names are different")
)


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

# 11. All microbe heatmap ------------------------------------------------------

plot_heatmap(phyloseq_obj = micro_phy, description = "microbe")


# 12. Metabolite eulerr diagram -----------------------------------------------

# merge samples by group
chem_type_merge <- merge_by_group(chem_phy, "sample_type")

chem_type_merge[chem_type_merge >0] <- 1

chem_type_merge <- t(chem_type_merge)


type_eul <- euler(chem_type_merge)

pdf("output/Euler/Metabolite_Sample_Type_Euler.pdf")
plot(type_eul, fills = sample_type_cols, quantities = T)
dev.off()
