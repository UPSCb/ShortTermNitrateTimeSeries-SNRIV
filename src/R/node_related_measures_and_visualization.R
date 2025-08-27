library(here)
library(igraph)
library(dplyr) 
library(tidyverse)
library(asnipe)
library(purrr)

g <- readRDS(here("data/seidr/clustering/graph.rds"))
g <- readRDS(here("data/seidr/clustering/graph.rds"))
degree(g)

set.seed(10)
de <- igraph::degree(g)
st <- igraph::strength(g)
be <- readRDS(here("data/seidr/clustering/betweeness.rds"))

names=V(g)$name

d=data.frame(node.name=names, degree=de, strength=st, betweeness=be) 
head(d) #display first 6 lines of data

plot(strength~degree, data=d)
# In our case strength=betweenness

plot(betweeness~degree, data=d)

# Get interesting genes

genes_interest <- d %>% 
  filter(betweeness > 0.0010 | degree > 700)

bb9 <- read_tsv(here("data/seidr/clustering/filtered_backbone-9-percent.tsv"),
                col_names=T,col_types=cols(.default=col_character()),
                show_col_types=FALSE)

genes_interest <- genes_interest %>%
  left_join(dplyr::select(bb9, Source, Source_Cluster), by=c("node.name" = "Source")) %>% 
  distinct()

write_tsv(genes_interest, here("data/seidr/clustering/genes_interest_centrality_measures/genes_interest_whole_network.tsv"))  

genes_all <- d %>%
  left_join(dplyr::select(bb9, Source, Source_Cluster), by=c("node.name" = "Source")) %>% 
  distinct()

write_tsv(genes_all, here("data/seidr/clustering/all_genes_network_stats.tsv"))  

# Now do single time point analysis
get_genes_interest <- function(graphml_file, output_tsv) {

gr <- igraph::read_graph(graphml_file,
                         format = "graphml")

set.seed(10)
de <- igraph::degree(gr)

names=V(gr)$name

d=data.frame(node.name=names, degree=de) 

violin <- ggplot(data = d, aes(x = basename(graphml_file), y= degree)) +
  geom_violin() +
  labs(x = "", y = "Degree")

genes_interest <- d %>%
  mutate(
    percentile99 = degree > quantile(d$degree, 0.995),
    Top20 = degree %in% head(base::sort(d$degree, decreasing=TRUE), n=20)
  ) %>% 
  filter(percentile99 & Top20)

# p <- plot(gr, vertex.label="", vertex.color="gold", edge.color="slateblue", vertex.size=de*2)

print(violin)
# print(p)

write_tsv(genes_interest, output_tsv)

return(genes_interest)

}

timepoints <- c("2h", "4h", "8h", "12h")
genes_interests <- purrr::map(timepoints, ~ get_genes_interest(
  here(paste0("data/seidr/clustering/firstDegreeNeighbours_", .x, ".graphml")),
       here(paste0("data/seidr/clustering/genes_interest_centrality_measures/genes_interest_firstDegreeNeighbours_", .x, ".graphml"))
))