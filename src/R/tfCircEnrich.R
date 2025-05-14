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
circReg <- read_excel("data/nitResponsiveTFnCircadian.xlsx", sheet = 4)
TFgene <- read_excel("data/nitResponsiveTFnCircadian.xlsx", sheet = 2)

cluster_info_ext <- cluster_info %>%
  mutate(IsCircClock = geneID %in% circClock$Potra,
         IsCircReg = geneID %in% circReg$Potra,
         IsTF = geneID %in% TFgene$Potra)
sum(cluster_info_ext$IsCircClock)
sum(cluster_info_ext$IsCircReg)
sum(cluster_info_ext$IsTF)

table(cluster_info_ext$Cluster,cluster_info_ext$IsCircClock)
table(cluster_info_ext$Cluster,cluster_info_ext$IsCircReg)
table(cluster_info_ext$Cluster,cluster_info_ext$IsTF)
table(cluster_info_ext$IsCircReg,cluster_info_ext$IsTF)
table(cluster_info_ext$IsCircClock,cluster_info_ext$IsTF)

FisherCircClock <- map(unique(cluster_info$Cluster), function(c){
  dat <- cluster_info_ext %>%
    mutate(IsCluster = Cluster %in% c)
  mat <- matrix(c(sum(dat$IsCluster & dat$IsCircClock), 
                  sum(dat$IsCluster & !dat$IsCircClock), 
                  sum(!dat$IsCluster & dat$IsCircClock), 
                  sum(!dat$IsCluster & !dat$IsCircClock)), 
                nrow = 2)
  return(broom::tidy(fisher.test(mat)))
}) %>%
  bind_rows() %>%
  mutate(Cluster = unique(cluster_info$Cluster),
         Sig = ifelse(p.value < 0.05, "Sig","ns"))

FisherCircReg <- map(unique(cluster_info$Cluster), function(c){
  dat <- cluster_info_ext %>%
    mutate(IsCluster = Cluster %in% c)
  mat <- matrix(c(sum(dat$IsCluster & dat$IsCircReg), 
                  sum(dat$IsCluster & !dat$IsCircReg), 
                  sum(!dat$IsCluster & dat$IsCircReg), 
                  sum(!dat$IsCluster & !dat$IsCircReg)), 
                nrow = 2)
  return(broom::tidy(fisher.test(mat)))
}) %>%
  bind_rows() %>%
  mutate(Cluster = unique(cluster_info$Cluster),
         Sig = ifelse(p.value < 0.001, "Sig","ns"))

FisherTF <- map(unique(cluster_info$Cluster), function(c){
  dat <- cluster_info_ext %>%
    mutate(IsCluster = Cluster %in% c)
  mat <- matrix(c(sum(dat$IsCluster & dat$IsTF), 
                  sum(dat$IsCluster & !dat$IsTF), 
                  sum(!dat$IsCluster & dat$IsTF), 
                  sum(!dat$IsCluster & !dat$IsTF)), 
                nrow = 2)
  return(broom::tidy(fisher.test(mat)))
}) %>%
  bind_rows() %>%
  mutate(Cluster = unique(cluster_info$Cluster),
         Sig = ifelse(p.value < 0.001, "Sig","ns"))

# Clus1  maybeCircClock, circReg, nitResTF
# Clus2  
# Clus3  circClock, circReg, nitResTF       4hrUP
# Clus4  nitResTF                           4hrDOWN
# Clus5  circReg, nitResTF                  4hrUP, 8hrDOWN
# Clus6  circReg                            2hrUP, 4hrUP
# Clus7  nitResTF
# Clus8  
# Clus9                                     4hrUP
# Clus10 

# 1. Start from DEG-firstDegreeNeighbor graphml
# 2. Filter for TF's firstDegreeeNeighbor from 1.
# 3. Check how many TFs left? For TFs that left, how many edges they have?
# 4. Check if those TFs that are left are CircReg or CircClock?

# 0. Cluster-wise GO from bb9


# results <- neighbors_cluster_df %>%
#   count(cluster) %>%
#   rowwise() %>%
#   mutate(
#     a = sum(neighbors_cluster_df$gene %in% circClock$Potra & neighbors_cluster_df$cluster == cluster),         
#     b = sum(!neighbors_cluster_df$gene %in% circClock$Potra & neighbors_cluster_df$cluster == cluster),       
#     c = sum(circClock$Potra %in% gene_with_cluster & !(circClock$Potra %in% neighbors_cluster_df$gene)), 
#     d = length(gene_with_cluster) - a - b - c,                                                           
#     pval = fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value,
#     odds_ratio = fisher.test(matrix(c(a, b, c, d), nrow = 2))$estimate
#   ) %>%
#   ungroup() %>%
#   arrange(pval)

print(results)
