library(here)
library(treeio)
library(ggtree)
library(readxl)
library(tidyverse)
library(ggtreeExtra)
library(ggnewscale)
library(ggpubr)
library(ggtext)
library(ggstar)
#library(tidytree)

here()

# Set colors ----
# colorblind palletes
colorBlindBlack8  <- c("#000000", "#0072B2", "#56B4E9", "#CC79A7", "#D55E00",
                       "#009E73", "#E69F00", "#F0E442")
# https://bioconductor.org/packages/release/bioc/vignettes/ggtreeExtra/inst/doc/ggtreeExtra.html#add-multiple-layers-on-the-same-position.

# Read in tree file ----
tree_name <- "withArrA_clustalo-I20220114-212957-0331-26091192-p2m.raxml.bestTree"
tree <- read.newick(file = here("data", tree_name))
attributes(tree)
str(tree)

## preview the tree ----
ggtree(tree) +
  geom_tiplab(size=3) +
  xlim(0,5)

## get bootstrap values if present and round to two decimals ----
d1 <- tree$node.label
if (!is.null(d1)) {
  d2 <- round(as.numeric(tree$node.label),2)
  d <- as.data.frame(cbind(label=d1, nlabel=d2))
  d[is.na(d)] = ""
  tree$node.label<-d[[2]][match(tree$node.label, d[[1]])]
}

# Read in metadata ----
meta_file <- "table_S2_SuppInfoFigure07.xlsx"
annotations <- read_xlsx(path = here("data", meta_file))


str(tree$tip.label)
tree_labels <- tree$tip.label
annotations_names <- annotations$name

## make sure tree tip names are the same as those in the annotation file ----
setdiff(tree_labels, annotations$name)

# Join tree data to annotations ----
tree_data <- as.treedata(tree)@phylo %>% as_tibble()

# Format the annotations using tree data ----
# Must have a column called "node"
tree_data_annotations <-
  right_join(tree_data, annotations,  by = c("label" = "name")) %>%
  mutate(short_name = paste0(short_name, sep = " ", number)) %>%
  column_to_rownames("label") %>%
  replace_na(list(Oxygen = "unknown", Water = "unknown")) %>%
  select(-tip_id, -gene_JGI_protein_NCBI, -amino_acid_sequence)

# Prepare annotations for geom_fruit ----
tree_data_annotations_fruit <- tree_data_annotations %>%
  mutate(arsenic_metab = case_when(Arsenotrophy == 1 ~ "ArxA-chemotrophy",
                                   Arsenotrophy == 2 ~ "ArxA-phototrophy",
                                   Arsenotrophy == 3 ~ "ArrA",
                                   TRUE ~ "ArxA-unknown")) %>%
  mutate(Arsenotrophy = factor(Arsenotrophy, levels = c(0,1,2,3)),
         taxonomy = as.factor(taxonomy),
         metagenome = as.factor(metagenome),
         Oxygen = factor(x=Oxygen, levels = c("aerobic", "micro", "anaerobic", "unknown")),
         Water = factor(x=Water, levels = c("freshwater","hot spring", "saline", "wastewater", "unknown")),
         bar = 1) %>%
  rownames_to_column(var="ID") %>%
  select(-number)

# Prepare annotation data tables ----
tree_data_annotations_labels <- tree_data_annotations %>%
  select(node, short_name, taxonomy) %>%
  rename(tax = taxonomy,
         shortname = short_name)

# Save tree data, this includes the node labels both terminal and internal labels
#table_out <- "table_S2.csv"
#tree_data_annotations_labels %>%
#  rownames_to_column("tip.label") %>%
#  write.csv(., file = here(here("data", table_out)), row.names = FALSE)


# Draw the tree ----
g1 <- ggtree(tree, size=.25) %<+% tree_data_annotations_labels +
  ylim(0, 70) +
  geom_tiplab(aes(label=shortname, color=tax),
              size=2.5,
              align=TRUE) +
  scale_color_manual(
    name="Taxonomy",
    values=colorBlindBlack8,
    guide="none"
  ) +
  scale_x_ggtree() +
  geom_fruit(
    data = tree_data_annotations_fruit,
    geom = geom_star,
    mapping = aes(y=ID, x=bar, fill=arsenic_metab, starshape=arsenic_metab),
    starstroke=0.25,
    size=2,
    pwidth=0.18,
    offset = 1.6,
    axis.params = list(axis="x", title="Arsenotrphy", title.angle=90)
    )  +
  scale_fill_manual(
    name="Arsenotrophy",
    values=c("#D81B60", "#1E88E5", "yellow", "grey95"),
    guide=guide_legend(keywidth=0.5, keyheight=0.5, order=1,
                       override.aes=list(size=2,
                                         starshape=c("ArrA"=1,
                                                     "ArxA-chemotrophy"=13,
                                                     "ArxA-phototrophy"=15,
                                                     "ArxA-unknown"=28))
    ),
    na.translate=FALSE,
  ) +
  scale_starshape_manual(
    values=c(1, 13, 15, 28),
    guide="none"
  )
g1
# set a new scale fill before making another geom_fruit
g2 <- g1 + new_scale_fill() +
  geom_fruit(
    data = tree_data_annotations_fruit,
    geom = geom_bar,
    mapping = aes(y=ID, x=bar, fill=taxonomy),
    orientation="y",
    stat="identity",
    color="white",
    pwidth=0.1,
    offset = 0.1,
    axis.params = list(axis="x", title="Taxonomy", title.angle=90)
  ) +
  scale_fill_manual(
    name="Taxonomy",
    values=colorBlindBlack8,
    guide=guide_legend(keywidth=0.7, keyheight=.6, order=2,
                       override.aes=list(size=0.25))
  )
# set a new scale fill before making another geom_fruit
g3 <- g2 + new_scale_fill() +
  geom_fruit(
    data = tree_data_annotations_fruit,
    geom = geom_bar,
    mapping = aes(y=ID, x=bar, fill=metagenome),
    orientation="y",
    stat="identity",
    color="white",
    pwidth=0.1,
    axis.params = list(axis="x", title="Source", title.angle=90)
  ) +
  scale_fill_viridis_d(name="Source",
                       option = "plasma",
                       guide=guide_legend(keywidth=0.7, keyheight=.6, order=3,
                                          override.aes=list(size=0.25))
  )
# set a new scale fill before making another geom_fruit
g4 <- g3 + new_scale_fill() +
  geom_fruit(
    data = tree_data_annotations_fruit,
    geom = geom_bar,
    mapping = aes(y=ID, x=bar, fill=Oxygen),
    orientation="y",
    stat="identity",
    color="white",
    pwidth=0.1,
    axis.params = list(axis="x", title="Oxygen", title.angle=90)
  ) +
  scale_fill_manual(name="Oxygen",
                       values = c("grey10", "grey60", "grey80", "grey95"),
                       guide=guide_legend(keywidth=0.7, keyheight=.6, order=4,
                                          override.aes=list(size=0.25))
  )
# set a new scale fill before making another geom_fruit
g5 <- g4 + new_scale_fill() +
  geom_fruit(
    data = tree_data_annotations_fruit,
    geom = geom_bar,
    mapping = aes(y=ID, x=bar, fill=Water),
    orientation="y",
    stat="identity",
    color="white",
    pwidth=0.1,
    axis.params = list(axis="x", title="Water", title.angle=90)
  ) +
  scale_fill_manual(name="Water",
                    values = c("#56B4E9", "#CC79A7", "#E69F00", "#009E73", "grey95"),
                    guide=guide_legend(keywidth=0.7, keyheight=.6, order=5,
                                       override.aes=list(size=0.25))
  )
# move legends to the left
g6 <- g5 + theme(legend.text=element_text(size=7),
           legend.title=element_text(size=8),
           legend.position = c(.12, .6),
           legend.background = element_blank()) +
  geom_treescale(x=0.1, y=12, width=0.2, fontsize = 2.5) # add scale bar

# g6 show the final tree ----
g6

# Print/save tree ----
ggsave(filename = here("results/figure_7_ORIO_genome_paper.svg"), plot=g6, height = 6.5, width = 5)
