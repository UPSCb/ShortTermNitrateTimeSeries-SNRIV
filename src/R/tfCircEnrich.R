library(stringr)
library(here)
library(igraph)
library(dplyr)
library(tidyverse)
library(asnipe)
library(purrr)
library(readxl)

g <- readRDS(here("data/seidr/clustering/graph.rds"))
degree(g)

set.seed(10)
de <- igraph::degree(g)
st <- igraph::strength(g)
be <- readRDS(here("data/seidr/clustering/betweeness.rds"))

names=V(g)$name

d=data.frame(node.name=names, degree=de, strength=st, betweeness=be) 
tfs <- read_tsv("data/seidr/nitResponsiveTFPotraID.txt", col_names = F)
tfs1 <- tfs %>% select(2) %>% unique()
colnames(tfs1) <- "Potra"

# tfs2 <- tfs %>% mutate(Sequence_Acc = str_remove(Sequence_Acc, "\\.\\d+$")) %>% 
#   distinct(Sequence_Acc)

tfInG <- tfs1 %>% left_join(d, by = c("Potra"= "node.name"))

tfInBb9 <- tfInG %>%
  left_join(dplyr::select(bb9, Source, Source_Cluster), by=c("Potra" = "Source")) %>% 
  distinct()

write_tsv(tfInBb9, here("data/tfInBb9Cluster.tsv"))

# restart R
# circadian gene enrichment in Clusters
library(dplyr)
library(purrr)
library(tibble)
library(readr)

bb9 <- read_tsv(here("data/seidr/clustering/filtered_backbone-9-percent.tsv"),
                col_names = TRUE, col_types = cols(.default = col_character()),
                show_col_types = FALSE)

cluster_info <- rbind(
  bb9[, c("Source", "Source_Cluster")] %>% setNames(c("geneID", "Cluster")),
  bb9[, c("Target", "Target_Cluster")] %>% setNames(c("geneID", "Cluster"))
) %>% distinct()

gene_to_cluster <- setNames(cluster_info$Cluster, cluster_info$geneID)

# first_degree_neighbor_genes <- readRDS("/pfs/stor10/users/home/s/shruti/project/ShortTermNitrateTimeSeries-SNRIV/data/seidr/clustering/first_degree_neighbor_genes.rds")
# 
# neighbors <- map(first_degree_neighbor_genes, function(timepoint) {
#   c(timepoint$up_neighbors, timepoint$down_neighbors)}) %>% unlist() %>% unique()
# 
# neighbors_cluster_df <- tibble(
#   gene = neighbors,
#   cluster = gene_to_cluster[neighbors]
# ) %>% filter(!is.na(cluster))

gene_with_cluster <- cluster_info$geneID

circClock <- read_excel("data/nitResponsiveTFnCircadian.xlsx", sheet = 3)

results <- neighbors_cluster_df %>%
  count(cluster) %>%
  rowwise() %>%
  mutate(
    a = sum(neighbors_cluster_df$gene %in% circClock$Potra & neighbors_cluster_df$cluster == cluster),         
    b = sum(!neighbors_cluster_df$gene %in% circClock$Potra & neighbors_cluster_df$cluster == cluster),       
    c = sum(circClock$Potra %in% gene_with_cluster & !(circClock$Potra %in% neighbors_cluster_df$gene)), 
    d = length(gene_with_cluster) - a - b - c,                                                           
    pval = fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value,
    odds_ratio = fisher.test(matrix(c(a, b, c, d), nrow = 2))$estimate
  ) %>%
  ungroup() %>%
  arrange(pval)

print(results)
