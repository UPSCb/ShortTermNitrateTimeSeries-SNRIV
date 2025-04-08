library(tidyverse)
library(here)



filter_data <- function(file_path) {
  # Read the TSV file
  data <- read_tsv(file_path) %>% 
    filter(!Source == Target)
  
  # Group by Source_Cluster and summarize counts
  result_Source <- data %>%
    group_by(Source_Cluster) %>%
    summarize(
      count_same = sum(Cluster == Source_Cluster),
      count_minus1 = sum(Cluster == -1)
    ) %>%
    mutate(
      to_keep = ifelse(count_same > count_minus1, "Keep", "Remove")
    ) %>% 
    filter(to_keep == "Remove")
  
  # Group by Target_Cluster and summarize counts
  result_Target <- data %>%
    group_by(Target_Cluster) %>%
    summarize(
      count_same = sum(Cluster == Target_Cluster),
      count_minus1 = sum(Cluster == -1)
    ) %>%
    mutate(
      to_keep = ifelse(count_same > count_minus1, "Keep", "Remove")
    ) %>% 
    filter(to_keep == "Remove")
  
  # Create the to_remove vector
  to_remove <- c(result_Source$Source_Cluster, result_Target$Target_Cluster) %>% unique()
  
  # Print the length of the to_remove vector
  print(paste("The number of removed clusters is ", length(to_remove)))
  
  # Filter the data
  data_filtered <- data %>% 
    filter(!Source_Cluster %in% to_remove) %>% 
    filter(!Target_Cluster %in% to_remove)
  
  # Print the number of removed rows
  print(paste("The number of removed edges is :", nrow(data) - nrow(data_filtered)))
  
  return(data_filtered)
}

data_1 <- filter_data(here("data/seidr/clustering/backbone-1-percent.tsv"))
  
data_5 <- filter_data(here("data/seidr/clustering/backbone-5-percent.tsv"))  

data_9 <- filter_data(here("data/seidr/clustering/backbone-9-percent.tsv"))  

data_1 %>%
  count(Cluster) %>% nrow()

data_5 %>%
  count(Cluster) %>% nrow()

data_9 %>%
  count(Cluster) %>% nrow()

data <- bind_rows(
  mutate(data_1, Dataset="data1"),
  mutate(data_5, Dataset="data5"),
  mutate(data_9, Dataset="data9")
)

data %>% 
  group_by(Dataset) %>% 
  summarize(unique_combined = n_distinct(union(Source, Target)))


data %>% 
  group_by(Dataset) %>% 
  summarize(unique_sources = n_distinct(union(Source_Cluster, Target_Cluster)))


# Backbone 5 has more genes (17782), but a high numbr of clusters (155)
# Backbone 9 has less genes (10962), but a low numbr of clusters (18)

write_tsv(data_5, here("data/seidr/clustering/filtered_backbone-5-percent.tsv"))
write_tsv(data_9, here("data/seidr/clustering/filtered_backbone-9-percent.tsv"))


# Calculate eigenegenes
library(WGCNA)
load(here("data/analysis/DE/vst-aware.rda"))

temp_Source <- data_9 %>% select(Source, Source_Cluster) %>% distinct() %>% 
  dplyr::rename(Gene = Source,
                Cluster = Source_Cluster)

temp_Target <- data_9 %>% select(Target, Target_Cluster) %>% distinct() %>% 
  dplyr::rename(Gene = Target,
                Cluster = Target_Cluster)

temp <- bind_rows(temp_Source, temp_Target) %>% distinct()

cluster_color <- temp %>% select(Cluster) %>% distinct()
cluster_color$Color <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdYlGn"))(nrow(cluster_color))

temp <- temp %>% left_join(cluster_color, by= "Cluster")

colors <- temp$Color
names(colors) <- temp$Gene

vst_filtered <- vst[rownames(vst) %in% temp$Gene,]


wgcna <- moduleEigengenes(t(vst_filtered), colors)

eigengenes <- t(wgcna$eigengenes)
rownames(eigengenes) <- gsub("^ME", "", rownames(eigengenes))

rownames(eigengenes) <- cluster_color$Cluster[match(rownames(eigengenes),cluster_color$Color)]

pheatmap(eigengenes)



# Check transcription factors of past interest because of obsession

interest <- c("Potra2n15c29002", "Potra2n12c23983")
data %>% filter(Source %in% interest) %>% 
  select(Dataset, Source, Type)  %>% distinct()

data %>% filter(Source %in% interest) %>% 
  group_by(Dataset)  %>% 
  summarize(
    edges = n()
  )



