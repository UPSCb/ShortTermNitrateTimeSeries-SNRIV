library(here)
library(igraph)
library(dplyr)
library(tidyverse)
library(tidyr)
library(Publish)
library(tibble)
library(readr)
library(readxl)
 
bb9 <- read_tsv(here("data/seidr/clustering/filtered_backbone-9-percent.tsv"),
               col_names=T,col_types=cols(.default=col_character()),
               show_col_types=FALSE)

bb9 <- bb9 %>%
  separate(irp_score.irp_rank, into = c("irp_score", "irp_rank"), sep = ";", 
           convert = TRUE)

# have been already filtered out before, still run this
bb9 <- bb9 %>% filter(!is.na(irp_score))

d.graf <- graph_from_edgelist(as.matrix(bb9[bb9$Type=="Directed",1:2]),directed=TRUE)
u1.graf <- graph_from_edgelist(as.matrix(bb9[bb9$Type=="Undirected",1:2]),directed=TRUE)
u2.graf <- graph_from_edgelist(as.matrix(bb9[bb9$Type=="Undirected",2:1]),directed=TRUE)
g <- igraph::union(d.graf,u1.graf,u2.graf)
saveRDS(g,"data/seidr/clustering/graph.rds")

pr <- page_rank(g)$vector
saveRDS(pr,"data/seidr/clustering/pageRank.rds")

# Run this once you go through the tutorial, example:
# https://dshizuka.github.io/networkanalysis/index.html

# IRP = inverse of weights
# E(g)$weight <- bb9$irp_score
# 
# #Check that we do not need to use "-"
# E(g)$inv_weight<- 1 / E(g)$weight

#   in_degree = degree(g, mode = "in"),
#   out_degree = degree(g, mode = "out"),
#   total_degree = degree(g, mode = "all"),
#   strength_in = strength(g, mode = "in", weights = E(g)$weight),
#   strength_out = strength(g, mode = "out", weights = E(g)$weight),
#   betweenness = betweenness(g, directed = T, weights = E(g)$inv_weight, normalized = T),
#   closeness_in = closeness(g, mode = "in", weights = E(g)$inv_weight, normalized = T),
#   closeness_out = closeness(g, mode = "out", weights = E(g)$inv_weight, normalized = T),
#   eigenvector = eigen_centrality(g, directed = TRUE, weights = E(g)$weight)$vector 

# Degree centrality =number of edges connected to a given node
de = igraph::degree(g)

# node strength=sum of the weights of edges connected to the node. 
st= strength(g)

# betweenness centrality= number of shortest paths thorugh a node. Nodes with 
# high betweenness might be influential in a network if, eg, they capture 
# most amount of information # flowing through the network because information 
# tends to flow through them.
be=betweenness(g, normalized=T)
saveRDS(be,"data/seidr/clustering/betweeness.rds")

# Assemble the node-level measures 
names=V(g)$name
d=data.frame(node.name=names, degree=de, strength=st, betweenness=be, pagerank=pr) 
saveRDS(d,"data/seidr/clustering/nodeMeasures.rds")

# first_neighbors_list <- lapply(V(g)$name, function(v) {
#   list(node = v, in_neighbors = names(neighbors(g, v, mode = "in")),
#        out_neighbors = names(neighbors(g, v, mode = "out")))})
# 
# neighbor_df <- map_dfr(first_neighbors_list, function(x) {
#   tibble(node = x$node, direction = c(rep("in", length(x$in_neighbors)),
#                                       rep("out", length(x$out_neighbors))),
#          neighbor = c(x$in_neighbors, x$out_neighbors))})

# Resrtart R
# different Timepoint is in separate sheet
g <- readRDS("/pfs/stor10/users/home/s/shruti/project/ShortTermNitrateTimeSeries-SNRIV/data/seidr/clustering/graph.rds")

deg_file <- "data/deg/allDeg.xlsx"
timepoint_map <- c("S1A_2h" = "2h", "S1B_4h" = "4h", "S1C_8h" = "8h", 
                   "S1D_12h" = "12h")

deg_split_list <- imap(timepoint_map, ~ {
  df <- read_excel(deg_file, sheet = .y) %>%
    select(Gene_Id, Log2_Fold_Change)
  list(up = df %>% filter(Log2_Fold_Change > 0),
       down = df %>% filter(Log2_Fold_Change < 0)
  )
}) %>% set_names(timepoint_map) 

#' Extract first degree neighbours of DEGs genes from network
#' 
# Use your dataframe or this
# subgrafs <- map(subgraf <- make_ego_graph(g,1, get.vertex.attribute(graf,"name") %in% goi))

# ego_list <- make_ego_graph(g, order = 1, nodes = V(g)[name %in% deg_split_list[["2h"]]$up$Gene_Id])

library(purrr)
get_neighbors <- function(goi, g) {
  ego_list <- make_ego_graph(g, order = 1, nodes = V(g)[name %in% goi])
  return(ego_list)
}

deg_neighbor_graph <- purrr::map(deg_split_list, ~ {
  list(up_neighbors = get_neighbors(.x$up$Gene_Id, g),
    down_neighbors = get_neighbors(.x$down$Gene_Id, g))})

first_degree_neighbour_genes <- map(deg_neighbor_graph, function(t) {
  map(t, function(n) {
    map(n, function(g) {
      names(V(g))
    }) %>%
      unlist() %>%
      unique()
  })
})

saveRDS(deg_neighbor_graph, "data/seidr/clustering/deg_neighbor_graph.rds")
saveRDS(first_degree_neighbour_genes, "data/seidr/clustering/first_degree_neighbor_genes.rds")
# purrr::walk(barplot(table(sapply(lapply(deg_neighbors,clusters),"[[","csize")),
#         labs=2,main="DEG cluster size",ylab="occurence",xlab="csize"))

# another code from Edoardo:
# extract_firstDegreeNeighbours <-  function(goi, graf) {
#   subgraf <- make_ego_graph(graf,1,
#                             get.vertex.attribute(graf,"name") %in% goi)
#   
#   gene_vector <- map(subgraf, ~ V(.x)$name) %>% unlist() %>% unique()
#   
# }

# list_of_DEGsAndNieghbours_vectors <- map(list_of_DEGs_vectors, ~ extract_firstDegreeNeighbours(.x, g))

#' combine all these networks together
combined_networks <- map(deg_neighbor_graph
                         Reduce("%u%",deg_neighbors))

#' Look at how many clusters we get and how many nodes are involved
# clusters(fdn)

#' Let's export the data for visualisation
walk(
write_graph(fdn,format = "graphml",file="firstDegreeNeighbour.graphml")
)

# what criteria is good for hubgenes, page rank?

# Then you see the single timepoint networks with cytoscape
# Then check if the genes with highest betweenness and pagre rank are in the DEGs

