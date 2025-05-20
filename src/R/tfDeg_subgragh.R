library(tidyverse)
library(here)
library(igraph)
library(asnipe)
library(readxl)

# bb9 graph
g <- readRDS("data/seidr/clustering/graph.rds")

# ref gene IDs
circClock <- read_excel("data/enrichment/nitResponsiveTFnCircadian.xlsx", sheet = 3)
circReg <- read_excel("data/enrichment/nitResponsiveTFnCircadian.xlsx", sheet = 4)
TFgene <- read_excel("data/enrichment/nitResponsiveTFnCircadian.xlsx", sheet = 2)
nitGene <- read_excel("data/enrichment/nitResTFnCircadian.xlsx", sheet = 2)

get_neighbors <- function(goi, g) {
  ego_list <- make_ego_graph(g, order = 1, nodes = V(g)[name %in% goi])
  names(ego_list) <- names(V(g)[name %in% goi])
  return(ego_list)}

TF_neighbor_graph_list <- get_neighbors(unique(TFgene$Potra), g)
TF_neighbor_graph_neighbors <- map(TF_neighbor_graph_list, ~ names(V(.x)))

deg2h <- unique((read_excel("data/analysis/DE/allDeg.xlsx", sheet = "S1A_2h")$Gene_Id))
deg4h <- unique((read_excel("data/analysis/DE/allDeg.xlsx", sheet = "S1B_4h")$Gene_Id))
deg8h <- unique((read_excel("data/analysis/DE/allDeg.xlsx", sheet = "S1C_8h")$Gene_Id))

TF_neighbor_deg_count <- map(TF_neighbor_graph_neighbors, function(n){
  return(c("2h" = sum(n %in% deg2h),
           "4h" = sum(n %in% deg4h),
           "8h" = sum(n %in% deg8h),
           "allNei" = length(n)))
}) %>% bind_rows(.id = "TF")

TF_neighbor_deg_count_filter <- TF_neighbor_deg_count %>%
  select(-allNei) %>%
  pivot_longer(cols = c("2h","4h","8h")) %>%
  filter(value > 0) %>%
  add_count(TF) %>%
  filter(n > 1) %>%
  select(-n) %>%
  pivot_wider()

TFinterest <- TF_neighbor_deg_count %>%
  filter(`2h` >= 10) %>%
  pull(TF)

graph_2hdeg_10neighbor <- Reduce("%u%", keep_at(TF_neighbor_graph_list,
                                      at = TFinterest)) %>%
  subgraph(vids = V(.)[name %in% c(TFinterest, deg2h)])

# plot(graph_2hdeg_10neighbor, layout=layout_with_fr(graph_2hdeg_10neighbor), vertex.label="", vertex.color="gold", edge.color="grey30",)
write_graph(graph_2hdeg_10neighbor,format = "graphml",file="data/seidr/graph_2hdeg_10neighbor.graphml")
write_graph(graph_2hdeg_10neighbor,format = "edgelist",file="data/seidr/graph_2hdeg_10neighbor.test")
