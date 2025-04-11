library(here)
library(igraph)
library(dplyr)
library(tidyverse)
library(tidyr)
library(purrr)
library(tibble)
library(readr)
 
bb9 <- read_tsv(here("data/seidr/clustering/filtered_backbone-9-percent.tsv"),
               col_names=T,col_types=cols(.default=col_character()),
               show_col_types=FALSE)

bb9 <- bb9 %>%
  separate(irp_score.irp_rank, into = c("irp_score", "irp_rank"), sep = ";", 
           convert = TRUE)

# have been already filtered out before, still run this
bb9 <- bb9 %>% filter(!is.na(irp_score))

# d.graf <- graph.edgelist(as.matrix(bb9[bb9$Type=="Directed",1:2]),directed=TRUE)

# u1.graf <- graph.edgelist(as.matrix(bb9[bb9$Type=="Undirected",1:2]),directed=TRUE)
# u2.graf <- graph.edgelist(as.matrix(bb9[bb9$Type=="Undirected",2:1]),directed=TRUE)
# g <- union(d.graf,u1.graf,u2.graf)

g <- graph.edgelist(as.matrix(bb9[bb9$Type=="Directed",1:2]),directed=TRUE)

pr <- page_rank(g)$vector
pr

# IRP = inverse of weights
E(g)$weight <- bb9$irp_score

#Check that we do not need to use "-"
E(g)$inv_weight<- 1 / E(g)$weight



# Following measures to calculate:
# Degree= how many direct connections 
# Strength=Weighted degree (e.g., based on IRP score)
# Betweenness	= node lies on shortest paths or gene in between clusters (e.g. "switches"),
# Closeness	of a node is to all others
# Eigenvector	Node’s influence based on neighbor quality
# PageRank? or strength to find highly connected genes (e.g. "modulators").
#  any other?

# Ask Nico to check this
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
  # Nico's function does not use the weights option for page rank.Should we use irp_score?
  pagerank = page_rank(g, directed = TRUE, weights = E(g)$weight)$vector,
  # optional
  eigenvector = eigen_centrality(g, directed = TRUE, weights = E(g)$weight)$vector 
  #there is a warning when I calculate with directed=T for eigencentrality, 
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


# Load DE genes for every Timepoint, separately
DEG_by_T <- list()
# goi <- sub("\\.1","",read.csv("DE-genes.csv",as.is = TRUE)[,1])

#' Extract the first degree neighbours of these genes from the network
#Use your dataframe or this
subgrafs <- map(
subgraf <- make_ego_graph(graf,1,
                          get.vertex.attribute(graf,"name") %in% goi)
)

walk(
barplot(table(sapply(lapply(subgraf,clusters),"[[","csize")),
        las=2,main="Gene of interest cluster size",
        ylab="occurence",xlab="csize")
)

#' combine all these networks together
combined_networks <- map(
fdn <- Reduce("%u%",subgraf)
)

#' Look at how many clusters we get and how many nodes are involved
# clusters(fdn)

#' Let's export the data for visualisation
walk(
write_graph(fdn,format = "graphml",file="firstDegreeNeighbour.graphml")
)

# what criteria is good for hubgenes, page rank?

# Then you see the single timepoint networks with cytoscape
# Then check if the genes with highest betweenness and pagre rank are in the DEGs

