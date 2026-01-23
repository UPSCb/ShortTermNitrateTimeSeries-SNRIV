#  This is for the GO enrichments in the gene coexpression networks and can be 
# adapted for another dataset
library(readr)
library(topGO)
library(dplyr)
library(tidyverse)
library(tidyr)
library(tibble)
library(here)

background <- readLines(here("data/analysis/DE/bg.csv"))

mapping <- read.delim(here("data/seidr/clustering/gene_go.txt"), stringsAsFactors = F)
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

# Go to A, B or C specific GOs
# A. For all the cluster specific genes 
bb9 <- read.delim(here("data/seidr/clustering/filtered_backbone-9-percent.tsv"), 
                  stringsAsFactors = F)

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

# with all vst expressed as background
results1 <- map(res, ~ topGO_combined(.x, background, goannot, alpha = 0.05, 
                                     p.adjust = "BH"))
results1 <- discard(results1, ~ nrow(.x) == 0)

dir.create("data/seidr/clustering/GO_cluster/")
walk2(results1, names(results1), ~ write_tsv(.x, file = here(paste0("data/seidr/clustering/GO_cluster/GO-enr_", .y, ".tsv"))))
saveRDS(results1, "data/seidr/clustering/GO_cluster/clusterGO.rds")

ntop <- 20

for (name in names(results1)) {
  suppressMessages(
    data <- results1[[name]] %>%
      dplyr::mutate(FDR = ifelse(FDR == "< 1e-30", "1e-30", FDR)) %>%
      dplyr::mutate(scores = -log10(as.numeric(FDR))) %>% 
      dplyr::select(Term, scores))
  data <- data %>% dplyr::arrange(desc(scores)) %>% dplyr::slice(1:ntop)
  data$Term <- factor(data$Term, levels = data$Term[order(data$scores, decreasing = TRUE)])
  
  p <- ggplot(data,
              aes(x = Term, y = scores, size = scores, fill = scores)) +
    expand_limits(y = 1) + geom_point(shape = 21) +
    scale_size(range = c(2.5,12.5), name="-log10(FDR)") +
    scale_fill_continuous(low = 'blue', high = 'red', name="-log10(FDR)") +
    xlab('') + ylab('-log10(FDR)') +
    labs(
      title = paste0("Enriched GOs for ", name),
      subtitle = paste('Top', ntop, 'terms ordered by adjustes pvalue'),
      caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001') +
    geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
               linetype = c("dotted", "longdash", "solid"),
               colour = c("black", "black", "black"),
               size = c(0.5, 1.5, 3)) +
    theme_bw(base_size = 24) +
    theme(
      legend.position = 'right', legend.background = element_rect(),
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
  
  ggsave(
    filename = file.path(here("data/seidr/clustering/GO_cluster/"), paste0("GO-enr_", name, ".png")),
    plot = p, width = 12, height = 8.5)
}

files <- list.files(path = "data/seidr/clustering/GO_cluster/", pattern = "\\.tsv$", full.names = TRUE)

combined_df <- lapply(files, function(file) {
  df <- read.delim(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  if ("FDR" %in% names(df)) {
    df$FDR <- suppressWarnings(as.numeric(df$FDR))
  }
  df$Cluster <- basename(file)
  
  return(df)
}) %>%
  bind_rows()

write.table(combined_df, "data/seidr/clustering/GO_cluster/combined_GO.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#A.2 with bb9 all genes as background
bb9bg <- unique(combined_df$Gene)
results2 <- map(res, ~ topGO_combined(.x, bb9bg, goannot, alpha = 0.05, 
                                      p.adjust = "BH"))
results2 <- discard(results2, ~ nrow(.x) == 0)

dir.create("data/seidr/clustering/GO_cluster_bb9bg/")
walk2 (results2, names(results2), ~ write_tsv(.x, file = here(paste0("data/seidr/clustering/GO_cluster_bb9bg/GO-enr_", .y, ".tsv"))))
saveRDS (results2, "data/seidr/clustering/GO_cluster_bb9bg/clusterGO.rds")

#merge GOs into one file
files <- list.files(path = "data/seidr/clustering/GO_cluster_bb9bg/", pattern = "\\.tsv$", full.names = TRUE)

combined_df <- lapply(files, function(file) {
  df <- read.delim(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    if ("FDR" %in% names(df)) {
    df$FDR <- suppressWarnings(as.numeric(df$FDR))
  }
  df$SourceFile <- basename(file)
  
  return(df)
}) %>%
  bind_rows()

write.table(combined_df, "data/seidr/clustering/GO_cluster_bb9bg/combined_GO.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# B. for the first degree neighbours of degs
first_degree_neighbour_genes <- load(here("data/seidr/clustering/first_degree_neighbor_genes.rds"))
res.list <- map(first_degree_neighbour_genes, ~ Filter(function(x) length(x) > 1, .x))
results <- map(res.list, ~ map(.x, ~ topGO_combined(.x, background, goannot, 
                                                    alpha = 0.05, 
                                                    p.adjust = "BH")))
results <- map(results, ~ discard(.x, ~ nrow(.x) == 0))
dir.create(here("data/seidr/clustering/GO_first_neighbors/"))
saveRDS(results, here("data/seidr/clustering/GO_first_neighbors/firstNeighborGO.rds"))

#merge GOs into one file
files <- list.files(path = "data/seidr/clustering/GO_first_neighbors/", pattern = "\\.tsv$", full.names = TRUE)

combined_df <- lapply(files, function(file) {
  df <- read.delim(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  if ("FDR" %in% names(df)) {
    df$FDR <- suppressWarnings(as.numeric(df$FDR))
  }
  df$SourceFile <- basename(file)
  
  return(df)
}) %>%
  bind_rows()

write.table(combined_df, "data/seidr/clustering/GO_first_neighbors/combined_GO.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

# gene ratio in the circles
for (name in names(results)) {
  for(nei in names(results[[name]])) {
    write_tsv(results[[name]][[nei]], 
              file = here(paste0(here("data/seidr/clustering/GO_first_neighbors/GO-enr_"), name, "_",nei ,".tsv")))
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
      filename = file.path(here("data/seidr/clustering/GO_first_neighbors/"), paste0("GO-enr_", name, "_", nei, ".png")),
      plot = p, width = 10, height = 6)
  }
}

# C. For specific genes of a timepoint of a specific cluster
# Restart R
library(topGO)
library(dplyr)
library(tidyverse)
library(tidyr)
library(here)
library(readxl)
library(dplyr)
library(stringr)
library(tibble)
library(purrr)

background <- readLines(here("data/analysis/DE/bg.csv"))

mapping <- read.delim(here("data/seidr/clustering/gene_go.txt"), stringsAsFactors = F)
prepAnnot <- function(annot){
  # Remove transcript version suffixes (e.g., .1, .2)
  annot$Sequence.Name <- sub("\\..*$", "", annot$Sequence.Name)
  # Split GO terms by semicolon
  geneID2GO <- strsplit(annot$Annotation.GO.ID, ";")
  names(geneID2GO) <- annot$Sequence.Name
  return(geneID2GO)
}

goannot <- prepAnnot(mapping)

bb9 <- read.delim(here("data/seidr/clustering/filtered_backbone-9-percent.tsv"), 
                  stringsAsFactors = F)

bb9_df1 <- bb9 %>% dplyr::select(Source, Source_Cluster) %>% rename(Gene=Source, Cluster=Source_Cluster)
bb9_df2 <- bb9 %>% dplyr::select(Target, Target_Cluster) %>% rename(Gene=Target, Cluster=Target_Cluster)
combined_df <- bind_rows(bb9_df1,bb9_df2) %>% distinct()
cluster_g <- split(combined_df$Gene, combined_df$Cluster)

deg_file <- here("data/analysis/DE/allDeg.xlsx")
timepoint_map <- c("S1A_2h" = "2h", "S1B_4h" = "4h", "S1C_8h" = "8h", 
                   "S1D_12h" = "12h")

deg_split_list <- imap(timepoint_map, ~ {
  df <- read_excel(deg_file, sheet = .y) %>%
    select(Gene_Id, Log2_Fold_Change)
  list(up = df %>% filter(Log2_Fold_Change > 0),
       down = df %>% filter(Log2_Fold_Change < 0)
  )
}) %>% set_names(timepoint_map) 

targets <- tribble(~Time, ~Regulation, ~Cluster,
  "2h", "up", "6", "2h", "up", "9", "2h", "down", "9",
  "4h", "up", "6", "4h", "up", "9", "4h", "up", "5", "4h", "down", "4",
  "8h", "down", "8")

matched_genes <- function(cluster_genes, deg_genes) {
  cluster_short <- str_sub(cluster_genes, 1, 12)
  deg_short <- str_sub(deg_genes, 1, 12)
  deg_genes[deg_short %in% cluster_short]}

target_deg_genes <- list()

for (i in seq_len(nrow(targets))) {
  time <- targets$Time[i]
  reg <- targets$Regulation[i]
  cluster <- targets$Cluster[i]
  
  cluster_name <- paste0("Cluster", cluster)
  genes_in_cluster <- cluster_g[[cluster_name]]
  deg_df <- deg_split_list[[time]][[reg]]
  
  matched_gene_ids <- matched_genes(genes_in_cluster, deg_df$Gene_Id)
  filtered_genes <- deg_df %>% filter(Gene_Id %in% matched_gene_ids)
  
  name <- paste(time, toupper(reg), cluster_name, sep = "_")
  target_deg_genes[[name]] <- filtered_genes
}

target_deg_genes <- Filter(function(x) length(x) > 1, target_deg_genes)

topGO_combined <- function(set,background,annotation,
                           ontology=c("BP","CC","MF"),
                           algorithm="parentchild",
                           statistic="fisher",
                           p.adjust=sort(p.adjust.methods),
                           alpha=0.05,
                           getgenes=FALSE){
  
  p.adjust <- match.arg(p.adjust)
  
  # create the allGenes
  allGenes <- factor(as.integer(background %in% set), levels = c(0, 1))
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

results3 <- map(target_deg_genes, ~ topGO_combined(.x, background, goannot, 
                                                   alpha = 0.05, 
                                                   p.adjust = "BH"))

#no signiifcant terms found