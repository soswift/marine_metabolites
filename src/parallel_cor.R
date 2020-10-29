# Run a bajillion correlations in parallel
library(parallel)

metabolite_file <- "data/processed/table_exports/paired_marine_metabolite_abundance_flat_table.csv"
microbe_file <- "data/processed/table_exports/paired_marine_microbe_abundance_flat_table.csv"

n_cores <- 4

# read in
chem_abund  <- read.csv(metabolite_file, header = T, row.names = 1)
micro_abund <- read.csv(microbe_file, header = T, row.names = 1)

# order by samples
chem_abund <- chem_abund[ , order(colnames(chem_abund))]
micro_abund <- micro_abund[ , order(colnames(micro_abund))]

# check samples match
all(colnames(micro_abund) == colnames(chem_abund))

# change each microbe into a list
micro_list <- data.frame(t(micro_abund))

# get_hardcor() runs correlations for a microbe on each row of a metabolite matrix
get_hardcor <- function(a_microbe) {
  apply(X = chem_abund, MARGIN = 1, function(x) {
    cor(x, a_microbe, method = "spearman")
  })
}

# run on all microbes
all_cors <- mclapply2(X = micro_list,
                      FUN = get_hardcor,
                      mc.cores = n_cores)


