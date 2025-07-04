library(tidyverse)
library(here)
library(igraph)
library(asnipe)
library(readxl)

# par(mar = c(1, 1, 1, 1))
g <- readRDS("data/seidr/clustering/graph.rds")

is_tf <- function(x) grepl("TF", x)
is_deg <- function(x) grepl("DEG", x)

TFgene <- unique(read_excel("data/enrichment/nitResponsiveTFnCircadian.xlsx", 
                            sheet = "tfCatPotra")$potra)

get_neighbors <- function(goi, g) {
  ego_list <- make_ego_graph(g, order = 1, nodes = V(g)[name %in% goi])
  names(ego_list) <- V(g)[name %in% goi]$name
  return(ego_list)
}

TF_neighbor_graph_list <- get_neighbors(TFgene, g)
TF_neighbor_graph_neighbors <- map(TF_neighbor_graph_list, ~ names(V(.x)))

deg2h <- unique(read_excel("data/analysis/DE/allDeg.xlsx", sheet = "S1A_2h")$Gene_Id)
deg4h <- unique(read_excel("data/analysis/DE/allDeg.xlsx", sheet = "S1B_4h")$Gene_Id)
deg8h <- unique(read_excel("data/analysis/DE/allDeg.xlsx", sheet = "S1C_8h")$Gene_Id)

TF_neighbor_deg_count <- map(TF_neighbor_graph_neighbors, function(n) {
  c("2h" = sum(n %in% deg2h),
    "4h" = sum(n %in% deg4h),
    "8h" = sum(n %in% deg8h),
    "allNei" = length(n))
}) %>% bind_rows(.id = "TF")

# Select TFs with >=10 neighbors DE at 2h
TFinterest <- TF_neighbor_deg_count %>% filter(`2h` >= 10) %>% pull(TF)
TFinterest <- TF_neighbor_deg_count %>% filter(`4h` >= 10) %>% pull(TF)
TFinterest <- TF_neighbor_deg_count %>% filter(`8h` >= 10) %>% pull(TF)

# Select TFs with >=1 neighbors DE at 2h, just save it don't use it
# TFinterest <- TF_neighbor_deg_count %>% filter(`2h` >= 1) %>% pull(TF)
# TFinterest <- TF_neighbor_deg_count %>% filter(`4h` >= 1) %>% pull(TF)
# TFinterest <- TF_neighbor_deg_count %>% filter(`8h` >= 1) %>% pull(TF)

# or Select the top 20 TFs' neighbors DE at 2h
# TFinterest <- TF_neighbor_deg_count %>% arrange(desc(`2h`)) %>% slice_head(n = 20) %>% pull(TF)

graph_2hdeg_10neighbor <- Reduce("%u%", keep_at(TF_neighbor_graph_list, 
                                                at = TFinterest)) %>%
  subgraph(vids = V(.)[name %in% c(TFinterest, deg2h)])

graph_4hdeg_10neighbor <- Reduce("%u%", keep_at(TF_neighbor_graph_list, 
                                                at = TFinterest)) %>%
  subgraph(vids = V(.)[name %in% c(TFinterest, deg4h)])

graph_8hdeg_10neighbor <- Reduce("%u%", keep_at(TF_neighbor_graph_list, 
                                                at = TFinterest)) %>%
  subgraph(vids = V(.)[name %in% c(TFinterest, deg8h)])

subgraph <- graph_2hdeg_10neighbor
subgraph <- graph_4hdeg_10neighbor
subgraph <- graph_8hdeg_10neighbor

annotations <- read.table("data/enrichment/classGene", header = TRUE, sep = "\t", 
                          stringsAsFactors = FALSE)
V(subgraph)$type <- annotations$Class[match(V(subgraph)$name, annotations$Potra)]

# Filter edges: TF → TF or DEG
edge_ends <- ends(subgraph, E(subgraph))
source_types <- V(subgraph)$type[match(edge_ends[, 1], V(subgraph)$name)]
target_types <- V(subgraph)$type[match(edge_ends[, 2], V(subgraph)$name)]
keep_edges <- E(subgraph)[is_tf(source_types) & (is_tf(target_types) | is_deg(target_types))]
subgraph <- subgraph.edges(subgraph, eids = keep_edges, delete.vertices = TRUE)

V(subgraph)$color <- ifelse(is_tf(V(subgraph)$type) & is_deg(V(subgraph)$type), "yellow",
                            ifelse(is_tf(V(subgraph)$type), "steelblue",
                                   ifelse(is_deg(V(subgraph)$type), "firebrick", "gray")))
V(subgraph)$shape <- ifelse(is_tf(V(subgraph)$type) & is_deg(V(subgraph)$type), "csquare",
                            ifelse(is_tf(V(subgraph)$type), "circle",
                                   ifelse(is_deg(V(subgraph)$type), "square", "triangle")))
par(mar = c(1, 1, 1, 1)) 
plot(subgraph, layout = layout_as_tree(subgraph), vertex.label = V(subgraph)$name,
     vertex.size = degree(subgraph) * 0.04 + 0.08, vertex.label.cex = 0.2,
     vertex.color = V(subgraph)$color, vertex.shape = V(subgraph)$shape,
     edge.arrow.size = 0.025)

# Create pure TF → DEG bipartite graph
edge_ends <- ends(subgraph, E(subgraph))
source_types <- V(subgraph)$type[match(edge_ends[,1], V(subgraph)$name)]
target_types <- V(subgraph)$type[match(edge_ends[,2], V(subgraph)$name)]
keep_edges <- E(subgraph)[is_tf(source_types) & is_deg(target_types)]
bipartite_graph <- subgraph.edges(subgraph, eids = keep_edges, delete.vertices = TRUE)

# Annotate roles
V(bipartite_graph)$is_tf <- is_tf(V(bipartite_graph)$type)
V(bipartite_graph)$is_deg <- is_deg(V(bipartite_graph)$type)

V(bipartite_graph)$color <- ifelse(V(bipartite_graph)$is_tf & V(bipartite_graph)$is_deg, "purple",
                                   ifelse(V(bipartite_graph)$is_tf, "firebrick", "forestgreen"))
V(bipartite_graph)$shape <- ifelse(V(bipartite_graph)$is_tf & V(bipartite_graph)$is_deg, "csquare",
                                   ifelse(V(bipartite_graph)$is_tf, "circle", "square"))

# Radial layout for bipartite TF → DEG
TFs <- V(bipartite_graph)[V(bipartite_graph)$is_tf]
layout_coords <- matrix(NA, nrow = vcount(bipartite_graph), ncol = 2)
rownames(layout_coords) <- V(bipartite_graph)$name

angle_increment <- 2 * pi / length(TFs)
for (i in seq_along(TFs)) {
  tf_name <- TFs[i]$name
  angle <- (i - 1) * angle_increment
  x <- cos(angle)
  y <- sin(angle)
  layout_coords[tf_name, ] <- c(x, y)
  
  degs <- neighbors(bipartite_graph, tf_name, mode = "out")
  degs <- degs[is_deg(V(bipartite_graph)[degs]$type)]
  
  for (j in seq_along(degs)) {
    deg_name <- V(bipartite_graph)[degs[j]]$name
    r <- 1.4 + j * 0.15
    layout_coords[deg_name, ] <- c(r * cos(angle), r * sin(angle))
  }
}

# Plot radial TF–DEG graph
plot(
  bipartite_graph,
  layout = layout_coords,
  vertex.label = V(bipartite_graph)$name,
  vertex.label.cex = 0.1,
  vertex.size = degree(bipartite_graph) * 0.2 + 0.4,
  vertex.color = V(bipartite_graph)$color,
  vertex.shape = V(bipartite_graph)$shape,
  edge.arrow.size = 0.02,
  main = "Radial TF–DEG Regulatory Network"
)

write_graph(subgraph,format = "graphml",file="data/enrichment/graph_TF10nei2hdeg.graphml")
write_graph(bipartite_graph,format = "graphml",file="data/enrichment/graph_TF10nei2hdegbipart.graphml")

write_graph(subgraph,format = "graphml",file="data/enrichment/graph_TF10nei4hdeg.graphml")
write_graph(bipartite_graph,format = "graphml",file="data/enrichment/graph_TF10nei4hdegbipart.graphml")

write_graph(subgraph,format = "graphml",file="data/enrichment/graph_TF10nei8hdeg.graphml")
write_graph(bipartite_graph,format = "graphml",file="data/enrichment/graph_TF10nei8hdegbipart.graphml")

plot(subgraph, layout = layout_as_tree(subgraph), vertex.label = V(subgraph)$name,
     vertex.size = degree(subgraph) * 0.04 + 0.08, vertex.label.cex = 0.2,
     vertex.color = V(subgraph)$color, vertex.shape = V(subgraph)$shape,
     edge.arrow.size = 0.025)

plot(subgraph, layout = layout_as_tree(subgraph), vertex.label = V(subgraph)$name,
     vertex.size = degree(subgraph) *  0.04 + 0.08, vertex.label.cex = 0.2,
     vertex.color = "steelblue", edge.arrow.size = 0.025)

plot(subgraph, layout = layout_nicely(subgraph), vertex.label = V(subgraph)$name,
     vertex.size = degree(subgraph) *  0.04 + 0.08, vertex.label.cex = 0.2,
     vertex.color = "steelblue", edge.arrow.size = 0.025)

plot(subgraph, layout = layout_on_sphere(subgraph), vertex.label = V(subgraph)$name,
     vertex.size = degree(subgraph) *  0.04 + 0.08, vertex.label.cex = 0.2,
     vertex.color = "steelblue", edge.arrow.size = 0.025)

plot(subgraph, layout = layout_randomly(subgraph), vertex.label = V(subgraph)$name,
     vertex.size = degree(subgraph) *  0.04 + 0.08, vertex.label.cex = 0.2,
     vertex.color = "steelblue", edge.arrow.size = 0.025)

# for hierarchical 
subgraph1 <- as.directed(subgraph, mode = "mutual")  # or "arbitrary"
layout <- layout_with_sugiyama(subgraph1)$layout
plot(subgraph1, layout = layout, vertex.label = V(subgraph)$name,
     vertex.size = degree(subgraph) *  0.04 + 0.08, vertex.label.cex = 0.2,
     vertex.color = "steelblue", edge.arrow.size = 0.025)

# layout <- layout_with_fr(subgraph)  # Fruchterman-Reingold for general network
# layout_with_kk(subgraph) #Kamada-Kawai for smaller graph
# layout_with_drl(subgraph) # Distributed Recursive Layout for large graphs.
# layout_in_circle (subgraph) # circle.
# layout_as_tree(graph, root = V(graph)[1]) #Tree layout (hierarchical) top-down structure.
# layout_on_grid(graph)# grid/ matrix-style visualizations).
# layout_with_lgl — Large Graph Layout # Fast for very large graphs.
# layout_with_graphopt # similar to force-directed, allows tuning with parameters.
# layout_with_sugiyama(graph)$layout # directed or layered networks
layouts <- list(fr = layout_with_fr(subgraph),
                kk = layout_with_kk(subgraph),
                dr= layout_with_drl(subgraph),
                grid= layout_on_grid(subgraph),
                circle = layout_in_circle(subgraph),
                tree = layout_as_tree(subgraph),
                sugiyama = layout_with_sugiyama(subgraph)$layout)

par(mfrow = c(3, 3))  # 2 rows, 3 columns
for (name in names(layouts)) {
  plot(subgraph,
       layout = layouts[[name]],
       main = paste("Layout:", name),
       vertex.label = NA,
       vertex.size = 5)
}