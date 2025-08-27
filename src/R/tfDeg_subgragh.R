library(tidyverse)
library(here)
library(igraph)
library(asnipe)
library(readxl)

# bb9 graph
g <- readRDS("data/seidr/clustering/graph.rds")

# all Potra TFs from plantgenie but remove the circadian clock genes and
# the best homolog of nitrate responsive # TFs from Varala et al 2018 and Vidal et al., 2020
tf <- unique(read_excel("data/enrichment/nitResponsiveTFnCircadian.xlsx", 
                            sheet = "tfCatPotra")$potra)

get_neighbors <- function(goi, g) {
  ego_list <- make_ego_graph(g, order = 1, nodes = V(g)[name %in% goi])
  names(ego_list) <- names(V(g)[name %in% goi])
  return(ego_list)}

TF_neighbor_graph_list <- get_neighbors(tf, g)
TF_neighbor_graph_neighbors <- map(TF_neighbor_graph_list, ~ names(V(.x)))

deg2h <- unique((read_excel("data/analysis/DEallDeg.xlsx", sheet = "S1A_2h")$Gene_Id))
deg4h <- unique((read_excel("data/analysis/DE/allDeg.xlsx", sheet = "S1B_4h")$Gene_Id))
deg8h <- unique((read_excel("data/analysis/DE/allDeg.xlsx", sheet = "S1C_8h")$Gene_Id))

# if you want the graphs time wise
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

TFinterest <- TF_neighbor_deg_count %>% filter(`2h` >= 10) %>%
  pull(TF)

# TFinterest <- TF_neighbor_deg_count %>% arrange(desc(`2h`)) %>%
#   slice_head(n = 20) %>% pull(TF)

graph_2hdeg_10neighbor <- Reduce("%u%", keep_at(TF_neighbor_graph_list,
                                                at = TFinterest)) %>%
  subgraph(vids = V(.)[name %in% c(TFinterest, deg2h)])

# plot(graph_2hdeg_10neighbor, layout=layout_with_fr(graph_2hdeg_10neighbor), 
# vertex.label="", vertex.color="gold", edge.color="grey30",)
write_graph(graph_2hdeg_10neighbor,format = "graphml",file="data/seidr/graph_2hdeg_10neighbor.graphml")

# Plot TFs in graph at different times
library(ggplot2)
library(pheatmap)

# summarize the number of tf-degs time-wise and tair family wise into plots
df <- read_excel("data/enrichment/tf-degInNetwork.xlsx",  sheet = 1)
df <- read_excel("data/enrichment/tf-degInNetwork.xlsx",  sheet = 2)
df <- read_excel("data/enrichment/tf-degInNetwork.xlsx",  sheet = 3)

df_long <- df %>% 
  # select(tair, `2h`, `4h`, `8h`) %>%
  select(family, `2h`, `4h`, `8h`) %>%
  pivot_longer(cols = c(`2h`, `4h`, `8h`), names_to = "Time", 
               values_to = "Number_DEGs") %>%
  group_by(family, Time) %>%
  summarize(Number_DEGs = sum(Number_DEGs, na.rm = TRUE), .groups = "drop")

d_filtered <- df_long %>%
  group_by(family) %>%
  filter(sum(Number_DEGs, na.rm = TRUE) > 0)

ggplot(df_filtered, aes(x = family, y = Number_DEGs, fill = Time)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme_minimal() +
  geom_text(aes(label = round(Number_DEGs, 2)), 
            position = position_dodge(width = 0.9), 
            vjust = -0.3, size = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")


heatmap_data <- df_long %>%
  pivot_wider(names_from = Time, values_from = Number_DEGs)

heatmap_data1 <- heatmap_data %>% select(`2h`,`4h`,`8h`)
rownames(heatmap_data1) <- heatmap_data$family

result <- pheatmap(heatmap_data1, cluster_rows = T, cluster_cols = F,
                   display_numbers = T, fontsize_row = 6, fontsize_number = 4,
                   color = colorRampPalette(c("white", "skyblue", "navy"))(50))

row_order <- rownames(heatmap_data1)[result$tree_row$order]
# write.csv(data.frame(Row_Clustering_Order = row_order), "row_clustering_order.csv", row.names = FALSE,
          # quote = F)

