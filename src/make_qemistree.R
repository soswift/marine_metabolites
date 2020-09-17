## use ggtree to reconstruct qemistree and highlight features
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(tidytree)
library(ape)
library(data.table)
library(ggstance)
library(ggtreeE)
# Read Tree ---------------------------------------------------

# read in tree file 
tree_raw <- read.tree(file = "data/raw/qemistree/qemistree.tree")


# Read Tables -------------------------------------------------

class_tips <- fread("data/raw/qemistree/labels.tsv", key = "id")
  
color_tips <- fread("data/raw/qemistree/colors.tsv", key = "id")

barplot_tips <- fread("data/raw/qemistree/barplots.tsv", key = "id")

tip_to_feature <- fread("data/raw/qemistree/Fingerprints to features.tsv", key = "id")

tip_data <- merge(class_tips, color_tips)
tip_data <- merge(tip_data, barplot_tips)
tip_data <- merge(tip_to_feature, tip_data)

long_tip_data <- melt(
  tip_data,
  id.vars = c("id"),
  measure.vars = c("CCA", "Coral", "Limu","Blank"),
  variable.name = "organism",
  value.name = "relabund"
) %>% as.data.frame()

tip_data <- as.data.frame(tip_data)
row.names(tip_data) <- tip_data$id

# Annotate Tree -----------------------------------------------


# try annotating with stacked barcharts

# try stacked barcharts, no tree (won't work because of row.names requirement)


ggplot(long_tip_data) +
  geom_bar(aes(x= id, y = relabund, fill = organism), stat = "identity", position = "stack")

# tree with heatmap of relative abundance for each organism

p1 <- ggtree(tree_raw, layout = "fan", open.angle = 30) %<+% tip_data +
  aes(color = superclass, show.legend = F)+
  theme_void()

# add heatmap
gheatmap(p1, tip_data[ , c("CCA","Coral","Limu","Blank")], offset = 1, width=0.6, 
         colnames=T, legend_title="Relative Abundance", colnames_offset_y = .25)


ggsave("output/tree/test_tree.png", width = 15, height = 20, units = "in")
            
# bargraph (for flat tree)  
facet_plot(p1, panel = "bargraph", data = long_tip_data, geom = geom_barh,
           mapping = aes(x = relabund, fill = organism, color = organism), stat = "identity", width = .1) 

