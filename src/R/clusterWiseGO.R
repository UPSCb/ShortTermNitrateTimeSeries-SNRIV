library(readr)
library(topGO)
library(dplyr)
library(tidyverse)
library(tidyr)
library(tibble)
library(here)

background <- readLines("data/deg/reducedBackgroundGenes.csv")

mapping <- read.delim("data/seidr/clustering/gene_go.txt", stringsAsFactors = F)
prepAnnot <- function(annot){
  # Remove transcript version suffixes (e.g., .1, .2)
  annot$Sequence.Name <- sub("\\..*$", "", annot$Sequence.Name)
  # Split GO terms by semicolon
  geneID2GO <- strsplit(annot$Annotation.GO.ID, ";")
  names(geneID2GO) <- annot$Sequence.Name
  return(geneID2GO)
}

goannot <- prepAnnot(mapping)

topGO_combined <- function(set,background,annotation,
                           ontology=c("BP","CC","MF"),
                           algorithm="parentchild",
                           statistic="fisher",
                           p.adjust=sort(p.adjust.methods),
                           alpha=0.05,
                           getgenes=FALSE){
  
  p.adjust <- match.arg(p.adjust)
  
  # create the allGenes
  allGenes <- factor(as.integer(background %in% set))
  names(allGenes) <- background
  
  # iterate over the ontologies
  lst <- lapply(ontology, function(o,g,a){
    GOdata <- new("topGOdata", 
                  ontology = o, 
                  allGenes = g,
                  annot = annFUN.gene2GO, 
                  gene2GO = a)
    results <- runTest(GOdata,algorithm=algorithm,statistic=statistic)
    
    n <- ifelse(p.adjust=="none",
                sum(score(results) <= alpha),
                sum(p.adjust(score(results),method=p.adjust) <= alpha))
    
    if(n==0){
      return(NULL)
    } else{
      #I added "numChar=1000 to not trim the GO
      resultTable <- as_tibble(GenTable(GOdata,
                                        results,
                                        numChar=1000,
                                        topNodes=n)) %>% 
        rename_with(function(sel){"FDR"},.cols=last_col())
      if(getgenes){
        resultTable <- resultTable %>%
          mutate(allgenes = map(GO.ID,function(x){allGO[[x]]})) %>%
          mutate(siggenes = map(allgenes,function(x){unlist(x)[unlist(x) %in% set]}))  %>%
          mutate(allgenes = map(allgenes,function(x){paste(x,collapse = "|")}) %>% unlist(use.names = F),
                 siggenes = map(siggenes,function(x){paste(x,collapse = "|")}) %>% unlist(use.names = F))
      }
      resultTable
    }
  },allGenes,annotation)
  names(lst) <- ontology
  return(bind_rows(lst, .id = "GO category"))
}

# A. For all the cluster specific genes 
bb9 <- read.delim("data/seidr/clustering/filtered_backbone-9-percent.tsv", 
                  stringsAsFactors = F)

# cluster_genes <- split(c(bb9$Source, bb9$Target), c(bb9$Source_Cluster, bb9$Target_Cluster))
# res.list <- lapply(cluster_genes, unique)
bb9_df1 <- bb9 %>% dplyr::select(Source, Source_Cluster) %>% rename(Gene=Source, Cluster=Source_Cluster)
bb9_df2 <- bb9 %>% dplyr::select(Target, Target_Cluster) %>% rename(Gene=Target, Cluster=Target_Cluster)
combined_df <- bind_rows(bb9_df1,bb9_df2) %>% distinct()
cluster_g <- split(combined_df$Gene, combined_df$Cluster)

res <- list()
for (name in names(cluster_g)) {
  res[[paste0(name)]] <- cluster_g[[name]]
  print(res)
}

res <- Filter(function(x) length(x) > 1, res)
results <- map(res, ~ topGO_combined(.x, background, goannot, alpha = 0.05, 
                                     p.adjust = "BH"))
results <- discard(results, ~ nrow(.x) == 0)
walk2(results, names(results), ~ write_tsv(.x, file = here(paste0("data/seidr/clustering/GO_cluster/GO-enr_", .y, ".tsv"))))
saveRDS(results, "data/seidr/clustering/clusterGO.rds")

ntop <- 20

for (name in names(results)) {
  suppressMessages(
    data <- results[[name]] %>%
      dplyr::mutate(FDR = ifelse(FDR == "< 1e-30", "1e-30", FDR)) %>%
      dplyr::mutate(scores = -log10(as.numeric(FDR))) %>% 
      dplyr::select(Term, scores))
  data <- data %>% dplyr::arrange(desc(scores)) %>% dplyr::slice(1:ntop)
  data$Term <- factor(data$Term, levels = data$Term[order(data$scores, decreasing = TRUE)])
  
  p <- ggplot(data, aes(x =Term, y = GeneRatio, color = FDR, size = Count)) +
    geom_point() + scale_color_gradient(low = "red", high = "blue") +
    theme_bw() + ylab("GeneRatio") + xlab("") +
    facet_grid(`GO category` ~ ., scales = "free", space = "free") +
    ggtitle(paste0(name," GO enrichment")) + coord_flip()
  
  plot(p)
  
  ggsave(
    filename = file.path("data/seidr/clustering/GO_cluster/", paste0("GO-enr_", name, ".png")),
    plot = p, width = 12, height = 8.5)
}

# B. for the first degree neighbours of degs
first_degree_neighbour_genes <- load("data/seidr/clustering/first_degree_neighbor_genes.rds")
res.list <- map(first_degree_neighbour_genes, ~ Filter(function(x) length(x) > 1, .x))
results <- map(res.list, ~ map(.x, ~ topGO_combined(.x, background, goannot, 
                                                    alpha = 0.05, 
                                                    p.adjust = "BH")))
results <- map(results, ~ discard(.x, ~ nrow(.x) == 0))
saveRDS(results, "data/seidr/clustering/firstNeighborGO.rds")

ntop <- 20

# p-value in the circles (optional)
# for (name in names(results)) {
  for(nei in names(results[[name]])) {
  suppressMessages(
    data <- results[[name]][[nei]] %>%
      dplyr::mutate(FDR = ifelse(FDR == "< 1e-30", "1e-30", FDR)) %>%
      dplyr::mutate(scores = -log10(as.numeric(FDR))) %>% 
      dplyr::select(Term, scores)
  )
  
  # Add this line to filter the top terms
  data <- data %>% dplyr::arrange(desc(scores)) %>% dplyr::slice(1:ntop)
  
  # Convert Term to a factor and specify the levels to be in the order of descending scores
  data$Term <- factor(data$Term, levels = data$Term[order(data$scores, decreasing = TRUE)])
  
  p <- ggplot(data,
              aes(x = Term, y = scores, size = scores, fill = scores)) +
    expand_limits(y = 1) +
    geom_point(shape = 21) +
    scale_size(range = c(2.5,12.5), name="-log10(FDR)") +
    scale_fill_continuous(low = 'royalblue', high = 'red4', name="-log10(FDR)") +
    xlab('') +
    ylab('-log10(FDR)') +
    labs(
      title = paste0("Enriched GOs for result ", name),
      subtitle = paste('Top', ntop, 'terms ordered by adjustes pvalue'),
      caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001') +
    geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
               linetype = c("dotted", "longdash", "solid"),
               colour = c("black", "black", "black"),
               size = c(0.5, 1.5, 3)) +
    theme_bw(base_size = 24) +
    theme(
      legend.position = 'right',
      legend.background = element_rect(),
      plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
      plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
      plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
      
      axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
      axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
      axis.title = element_text(size = 12, face = 'bold'),
      axis.title.x = element_text(size = 12, face = 'bold'),
      axis.title.y = element_text(size = 12, face = 'bold'),
      axis.line = element_line(colour = 'black'),
      
      #Legend
      legend.key = element_blank(), # removes the border
      legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
      legend.text = element_text(size = 14, face = "bold"), # Text size
      title = element_text(size = 14, face = "bold")) +
    coord_flip()
  
  print(p)
}
}

# gene ratio in the circles
for (name in names(results)) {
  for(nei in names(results[[name]])) {
    write_tsv(results[[name]][[nei]], 
              file = here(paste0("data/seidr/clustering/GO/GO-enr_", name, "_",nei ,".tsv")))
      dat <- results[[name]][[nei]] %>%
      mutate(FDR = parse_double(sub("<","",FDR)),
             GeneRatio = Significant/Annotated,
             Count = as.numeric(Significant)) %>%
    # slice(1:ntop) %>%
    mutate(Term = factor(Term, levels = Term[order(GeneRatio)]))
    p <- ggplot(dat, aes(x =Term, y = GeneRatio, color = FDR, size = Count)) +
      geom_point() + scale_color_gradient(low = "red", high = "blue") +
      theme_bw() + ylab("GeneRatio") + xlab("") +
      facet_grid(`GO category` ~ ., scales = "free", space = "free") +
      ggtitle(paste0(nei, ": ", name," GO enrichment")) + coord_flip()
    plot(p)
    ggsave(
      filename = file.path("data/seidr/clustering/GO/", paste0("GO-enr_", name, "_", nei, ".png")),
      plot = p, width = 10, height = 6)
  }
}
