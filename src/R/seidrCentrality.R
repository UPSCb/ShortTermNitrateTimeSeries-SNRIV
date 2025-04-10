library(here)
# couldn't get the follwing to run from seidrPageRank.R:
# library(igraph)
# library(readr)
# 
# sf <- read_tsv(here("data/seidr/clustering/filtered_backbone-9-percent.tsv"),
#                col_names=T,col_types=cols(.default=col_character()),
#                show_col_types=FALSE)
# 
# d.graf <- graph.edgelist(as.matrix(sf[sf$X3=="Directed",1:2]),directed=TRUE)
# I don't know what is happening here, don't get this to work:
# u1.graf <- graph.edgelist(as.matrix(sf[sf$X3=="Undirected",1:2]),directed=TRUE)
# u2.graf <- graph.edgelist(as.matrix(sf[sf$X3=="Undirected",2:1]),directed=TRUE)
# graf <- union(d.graf,u1.graf,u2.graf)
# 
# pr <- page_rank(graf)$vector
# pr

library(igraph)
library(dplyr)
library(tidyverse)
library(tidyr)
library(purr)
library(tibble)

bb9 <- read.delim("data/seidr/clustering/filtered_backbone-9-percent.tsv", 
                  stringsAsFactors = FALSE)
bb9_beaut <- bb9 %>%
  separate(irp_score.irp_rank, into = c("irp_score", "irp_rank"), sep = ";", 
           convert = TRUE)

# have been already filtered out before, still run this
bb9_beaut <- bb9_beaut %>% filter(!is.na(irp_score))

# g <- graph_from_data_frame(bb9_beaut, directed = F)
g <- graph_from_data_frame(bb9_beaut, directed = T)

# IRP = inverse of weights
E(g)$weight <- bb9_beaut$irp_score

E(g)$inv_weight<- 1 / E(g)$weight

# Following measures to calculate:
# Degree= how many direct connections 
# Strength=Weighted degree (e.g., based on IRP score)
# Betweenness	= node lies on shortest paths or gene in between clusters (e.g. "switches"),
# Closeness	of a node is to all others
# Eigenvector	Node’s influence based on neighbor quality
# PageRank? or strength to find highly connected genes (e.g. "modulators").
#  any other?

centrality <- data.frame(
  node = V(g)$name,
  in_degree = degree(g, mode = "in"),
  out_degree = degree(g, mode = "out"),
  total_degree = degree(g, mode = "all"),
  strength_in = strength(g, mode = "in", weights = E(g)$weight),
  strength_out = strength(g, mode = "out", weights = E(g)$weight),
  
  betweenness = betweenness(g, directed = T, weights = E(g)$inv_weight, normalized = T),
  closeness_in = closeness(g, mode = "in", weights = E(g)$inv_weight, normalized = T),
  closeness_out = closeness(g, mode = "out", weights = E(g)$inv_weight, normalized = T),
  
  pagerank = page_rank(g, directed = TRUE, weights = E(g)$weight)$vector,
  # optional
  eigenvector = eigen_centrality(g, directed = TRUE, weights = E(g)$weight)$vector 
  #there is a warning when I calculate with directed=T for eigencentrality, 
  # but eigen we already have 
)

first_neighbors_list <- lapply(V(g)$name, function(v) {
  list(node = v, in_neighbors = names(neighbors(g, v, mode = "in")),
       out_neighbors = names(neighbors(g, v, mode = "out")))
})

write.csv(centrality, "data/seidr/clustering/bb9_centrality_directed.csv", row.names = FALSE)
# write.csv(centrality, "data/seidr/clustering/bb9_centrality_undirected.csv", row.names = FALSE)

neighbor_df <- map_dfr(first_neighbors_list, function(x) {
  tibble(node = x$node, direction = c(rep("in", length(x$in_neighbors)), 
                                      rep("out", length(x$out_neighbors))),
         neighbor = c(x$in_neighbors, x$out_neighbors))})

# what criteria is good for hubgenes, page rank?

