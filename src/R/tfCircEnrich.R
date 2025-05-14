library(stringr)
library(here)
library(igraph)
library(dplyr)
library(tidyverse)
library(asnipe)
library(purrr)
library(readxl)
library(tibble)
library(readr)

# circadian and TF enrichment in Clusters
bb9 <- read_tsv(here("data/seidr/clustering/filtered_backbone-9-percent.tsv"),
                col_names = TRUE, col_types = cols(.default = col_character()),
                show_col_types = FALSE)

cluster_info <- rbind(
  bb9[, c("Source", "Source_Cluster")] %>% setNames(c("geneID", "Cluster")),
  bb9[, c("Target", "Target_Cluster")] %>% setNames(c("geneID", "Cluster"))
) %>% distinct()

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

write.table(FisherTF, "data/fisherTF.txt", quote = F, sep = "\t", 
            row.names = F, col.names = T)
write.table(FisherCircReg, "data/fisherCircReg.txt", quote = F, sep = "\t", 
            row.names = F, col.names = T)
write.table(FisherCircClock, "data/fisherCircClock.txt", quote = F, sep = "\t", 
            row.names = F, col.names = T)

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

# 0. Cluster-wise GO from bb9
# 1. Start from DEG-firstDegreeNeighbor graphml
# 2. Filter for TF's firstDegreeeNeighbor from 1.
# 3. Check how many TFs left? For TFs that left, how many edges they have?
# 4. Check if those TFs that are left are CircReg or CircClock?

library(igraph)
library(readxl)

g <- read_graph("data/seidr/clustering/firstDegreeNeighbours_4h.graphml",
                format = "graphml")
TFgene <- unique(read_excel("data/nitResponsiveTFnCircadian.xlsx", sheet = 2)$Potra)
deg4h <- unique((read_excel("data/analysis/DE/allDeg.xlsx", sheet = "S1B_4h")$Gene_Id))

edge_list <- ends(g, E(g), names = TRUE)

tf_deg_indices <- which(
  (edge_list[,1] %in% tf_list & edge_list[,2] %in% deg_list) |
    (edge_list[,1] %in% deg_list & edge_list[,2] %in% tf_list)
)

g_tf_deg <- subgraph.edges(g, E(g)[tf_deg_indices], delete.vertices = FALSE)

