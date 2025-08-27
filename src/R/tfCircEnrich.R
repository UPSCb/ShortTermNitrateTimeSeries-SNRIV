# circadian and TF enrichment in Clusters
library(tidyverse)
library(here)
library(igraph)
library(asnipe)
library(readxl)

bb9 <- read_tsv(here("data/seidr/clustering/filtered_backbone-9-percent.tsv"),
                col_names = TRUE, col_types = cols(.default = col_character()),
                show_col_types = FALSE)

cluster_info <- rbind(
  bb9[, c("Source", "Source_Cluster")] %>% setNames(c("geneID", "Cluster")),
  bb9[, c("Target", "Target_Cluster")] %>% setNames(c("geneID", "Cluster"))
) %>% distinct()

# AT homologs for circadian clock and circadian regulation in Potra
circClock <- read_excel("data/enrichment/nitResponsiveTFnCircadian.xlsx", sheet = 3)
circReg <- read_excel("data/enrichment/nitResponsiveTFnCircadian.xlsx", sheet = 4)

# all Potra TFs from plantgenie but remove the circadian clock genes and
# the best homolog of nitrate responsive # TFs from Varala et al 2018 and Vidal et al., 2020
TFgene <- read_excel("data/enrichment/nitResponsiveTFnCircadian.xlsx", sheet = 2)

# Nitrogen metabolism genes
nitGene <- read_excel("data/enrichment/nitResponsiveTFnCircadian.xlsx", sheet = 5)

cluster_info_ext <- cluster_info %>%
  mutate(IsCircClock = geneID %in% circClock$Potra,
         IsCircReg = geneID %in% circReg$Potra,
         IsTF = geneID %in% TFgene$Potra,
         IsNitGene = geneID %in% nitGene$potra)

sum(cluster_info_ext$IsCircClock)
sum(cluster_info_ext$IsCircReg)
sum(cluster_info_ext$IsTF)
sum(cluster_info_ext$IsNitGene)

table(cluster_info_ext$Cluster,cluster_info_ext$IsCircClock)
table(cluster_info_ext$Cluster,cluster_info_ext$IsCircReg)
table(cluster_info_ext$Cluster,cluster_info_ext$IsTF)
table(cluster_info_ext$IsCircReg,cluster_info_ext$IsTF)
table(cluster_info_ext$IsCircClock,cluster_info_ext$IsTF)
table(cluster_info_ext$Cluster,cluster_info_ext$IsNitGene)
table(cluster_info_ext$IsCircReg,cluster_info_ext$IsNitGene)

FisherCircClock <- map(unique(cluster_info$Cluster), function(c){
  dat <- cluster_info_ext %>%
    mutate(IsCluster = Cluster %in% c)
  mat <- matrix(c(sum(dat$IsCluster & dat$IsCircClock), 
                  sum(dat$IsCluster & !dat$IsCircClock), 
                  sum(!dat$IsCluster & dat$IsCircClock), 
                  sum(!dat$IsCluster & !dat$IsCircClock)), 
                nrow = 2)
  return(broom::tidy(fisher.test(mat)) %>% mutate(Cluster = c))
}) %>%
  bind_rows() %>%
  mutate(Sig = ifelse(p.value < 0.05, "Sig","ns"))

FisherCircReg <- map(unique(cluster_info$Cluster), function(c){
  dat <- cluster_info_ext %>%
    mutate(IsCluster = Cluster %in% c)
  mat <- matrix(c(sum(dat$IsCluster & dat$IsCircReg), 
                  sum(dat$IsCluster & !dat$IsCircReg), 
                  sum(!dat$IsCluster & dat$IsCircReg), 
                  sum(!dat$IsCluster & !dat$IsCircReg)), 
                nrow = 2)
  return(broom::tidy(fisher.test(mat)) %>% mutate(Cluster = c))
}) %>%
  bind_rows() %>%
  mutate(Sig = ifelse(p.value < 0.001, "Sig","ns"))

FisherTF <- map(unique(cluster_info$Cluster), function(c){
  dat <- cluster_info_ext %>%
    mutate(IsCluster = Cluster %in% c)
  mat <- matrix(c(sum(dat$IsCluster & dat$IsTF), 
                  sum(dat$IsCluster & !dat$IsTF), 
                  sum(!dat$IsCluster & dat$IsTF), 
                  sum(!dat$IsCluster & !dat$IsTF)), 
                nrow = 2)
  return(broom::tidy(fisher.test(mat)) %>% mutate(Cluster = c))
}) %>%
  bind_rows() %>%
  mutate(Sig = ifelse(p.value < 0.001, "Sig","ns"))

FisherNitGene <- map(unique(cluster_info$Cluster), function(c){
  dat <- cluster_info_ext %>%
    mutate(IsCluster = Cluster %in% c)
  mat <- matrix(c(sum(dat$IsCluster & dat$IsNitGene), 
                  sum(dat$IsCluster & !dat$IsNitGene), 
                  sum(!dat$IsCluster & dat$IsNitGene), 
                  sum(!dat$IsCluster & !dat$IsNitGene)), 
                nrow = 2)
  return(broom::tidy(fisher.test(mat)) %>% mutate(Cluster = c))
}) %>%
  bind_rows() %>%
  mutate(Sig = ifelse(p.value < 0.001, "Sig","ns"))

write.table(FisherTF, "data/enrichment/fisherTF.txt", quote = F, sep = "\t", 
            row.names = F, col.names = T)
write.table(FisherCircReg, "data/enrichment/fisherCircReg.txt", quote = F, sep = "\t", 
            row.names = F, col.names = T)
write.table(FisherCircClock, "data/enrichment/fisherCircClock.txt", quote = F, sep = "\t", 
            row.names = F, col.names = T)
write.table(FisherNitGene, "data/enrichment/fisherNitGene.txt", quote = F, sep = "\t", 
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

# Restart R 
# Check for the cluster with interesting TFs in different timepoints
library(igraph)
library(readxl)

TFgene <- unique(read_excel("data/enrichment/nitResponsiveTFnCircadian.xlsx", sheet = 2)$Potra)
deg2h <- unique((read_excel("data/analysis/DE/allDeg.xlsx", sheet = "S1A_2h")$Gene_Id))
deg4h <- unique((read_excel("data/analysis/DE/allDeg.xlsx", sheet = "S1B_4h")$Gene_Id))
deg8h <- unique((read_excel("data/analysis/DE/allDeg.xlsx", sheet = "S1C_8h")$Gene_Id))

get_neighbors <- function(goi, g) {
  ego_list <- make_ego_graph(g, order = 1, nodes = V(g)[name %in% goi])
  return(ego_list)}

# g <- read_graph("data/seidr/clustering/firstDegreeNeighbours_2h.graphml", format = "graphml")
# g <- read_graph("data/seidr/clustering/firstDegreeNeighbours_4h.graphml", format = "graphml")
# g <- read_graph("data/seidr/clustering/firstDegreeNeighbours_8h.graphml", format = "graphml")

# TF_neighbor_graph_list <- get_neighbors(unique(TFgene), g)
# TF_neighbor_graph <- Reduce("%u%",TF_neighbor_graph_list)
 
# saveRDS(TF_neighbor_graph, "data/enrichment/TFsubgraph_2h.rds")
 
# TF_deg_graph <- subgraph(TF_neighbor_graph, vids = V(TF_neighbor_graph)[name %in% c(TFgene,deg4h)])
# degree(TF_deg_graph)[names(degree(TF_deg_graph)) %in% TFgene] %>%
#   as.data.frame()

timepoints <- c("2h", "4h", "8h")
deg_vars <- list(deg2h = deg2h, deg4h = deg4h, deg8h = deg8h)

library(dplyr)
library(tidyr)
library(purrr)

TF_degree_dfs <- map(timepoints, function(tp) {
  g <- read_graph(paste0("data/seidr/clustering/firstDegreeNeighbours_", tp, ".graphml"), format = "graphml")
  
  TF_neighbors <- get_neighbors(unique(TFgene), g)
  TF_graph <- Reduce("%u%", TF_neighbors)
  
  # saveRDS(TF_graph, paste0("data/enrichment/TFsubgraph_", tp, ".rds"))
  
  deg_var <- deg_vars[[tp]]
  TF_deg_graph <- subgraph(TF_graph, vids = V(TF_graph)[name %in% c(TFgene, deg_var)])
  degs <- degree(TF_deg_graph)
  
  data.frame(
    TF = names(degs)[names(degs) %in% TFgene],
    degree = degs[names(degs) %in% TFgene],
    source = tp
  )
})

all_degrees <- bind_rows(TF_degree_dfs)
common_TFs <- reduce(map(TF_degree_dfs, ~ .x$TF), intersect)

common_degrees_wide <- all_degrees %>% filter(TF %in% common_TFs) %>%
  pivot_wider(names_from = source, values_from = degree)

subgraph_TFs <- intersect(V(TFsubgraph_8h)$name, TFgene)

circClock <- unique(read_excel("data/enrichment/nitResponsiveTFnCircadian.xlsx", sheet = 3)$Potra)
circReg <- unique(read_excel("data/enrichment/nitResponsiveTFnCircadian.xlsx", sheet = 4)$Potra)

summary(subgraph_TFs %in% circReg)
#     FALSE TRUE 
# 4h  329   66 
# 2h  173   45 
# 8h  176   38 

summary(subgraph_TFs %in% circClock)
#     FALSE  TRUE
# 2h  213 5 
# 4h 391 4 
# 8h  210 4

edge_list <- ends(TFsubgraph_2h, E(TFsubgraph_2h), names = TRUE)
tf_edge_list <- edge_list[edge_list[,1] %in% TFgene | edge_list[,2] %in% TFgene,]
nrow(tf_edge_list)

