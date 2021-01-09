## Load libraries --------------------------
# General utility
library(phyloseq)
library(data.table)
library(tidyr)
library(biomformat)
# Plotting and analysis
library(eulerr)
library(vegan)
library(pairwiseAdonis)

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
# Save Point: Load Cleaned Microbe and Metabolite Data----------------------------------------------------------------
final_marine_phy <- readRDS("data/processed/final_marine_phy.rds")
micro_phy <- final_marine_phy

final_unifrac <- readRDS("data/processed/final_unifrac.rds")

chem_phy <- readRDS("data/processed/chem_phy.rds")


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

