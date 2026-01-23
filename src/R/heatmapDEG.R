#' ---
#' title: "Differential Expression"
#' author: "Shruti"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#' ---
#' # Setup

#' * Libraries
setwd("/mnt/picea/home/schoudhary/shruti/ShortTermNitrateTimeSeries-SNRIV/")
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(gplots)
  library(ggplot2)
  library(here)
  library(hyperSpec)
  library(RColorBrewer)
  library(tidyverse)
  library(VennDiagram)
  library(readxl)
})

#' * Graphics
pal=brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

# Define your contrasts as a list of vectors
dds <- readRDS(here("data/dds2.rds"))
samples <- read.table(here("doc/sampleTremula.csv"), header = T)
load(here("data/analysis/DE/vst-aware.rda"))

# Take average across all replicates per sample
treatmentTime <- c("None_0h", "KCL_2h", "KNO3_2h","KCL_4h", "KNO3_4h", 
                   "KCL_8h", "KNO3_8h", "KCL_12h", "KNO3_12h", 
                   "KCL_24h","KNO3_24h", "KCL_48h", "KNO3_48h")

allAve <- sapply(treatmentTime, function(level) as.matrix(rowMeans(vst[, dds$treatmentTime == level])))
rownames(allAve) <- rownames(vst)
allAve <- allAve[rowSums(allAve) != 0, ]

# add gene symbols
atnnotation <- read.table(here("doc/potraATgenes.txt"), header = T, sep = "\t")
notation <- atnnotation[match(rownames(allAve), atnnotation$Potra),]
all(rownames(allAve) == notation$Potra)
rownames(allAve) <- paste(rownames(allAve), notation$AT_best_hit, sep = " ")

# selective genes heatmap
ligbiosyn <- c("Potra2n19c33285","Potra2n1c351","Potra2n10c20411","Potra2n9c19275",
               "Potra2n9c19246","Potra2n8c17885","Potra2n8c16946","Potra2n7c16495",
               "Potra2n6c15127","Potra2n6c14201","Potra2n5c11765","Potra2n3c7945",
               "Potra2n3c6847","Potra2n3c6817","Potra2n3c6783","Potra2n1c307",
               "Potra2n1c2649","Potra2n1c1541","Potra2n16c30091","Potra2n16c29966",
               "Potra2n13c24959","Potra2n12c23766")

dir.create(here("data/intGeneHeatmap"),showWarnings=FALSE)

hmap <- function(selGene, file_name) {
  selGene <- ifelse(selGene %in% notation$Potra, paste(selGene, notation$AT_best_hit[match(selGene, notation$Potra)]), selGene)
  vst2 <- allAve[rownames(allAve) %in% selGene, ]
  svg(file.path(".",paste0(file_name,".svg")), pointsize = 8)
  heatmap.2(t(scale(t(vst2))), distfun = pearson.dist,
            hclustfun = function(X){hclust(X,method="ward.D2")},
            trace="none", col=hpal, margins =c(12,12), cexCol = 1,
            cexRow = 0.3, main = file_name, key = TRUE, keysize = 1,
            Colv = FALSE, Rowv = T
            # labRow = paste(rownames(vst2), nmet$Gene.family[match(rownames(vst2), selGene)])
  )
  dev.off()
}

hmap(ligbiosyn, "lig")

# for time Bins heatmap in suppl figure
timeBins <- read.csv(here("doc/timeBins_fdr0.05_lfc1.0.csv"), header = T)

gene_set_timeBin <- list(
  uptimeBin2h = timeBins %>% filter(type == "upregulated" & filename == "KNO3_2h_vs_KCL_2h") %>% pull(genes),
  uptimeBin4h = timeBins %>% filter(type == "upregulated" & filename == "KNO3_4h_vs_KCL_4h") %>% pull(genes),
  uptimeBin8h = timeBins %>% filter(type == "upregulated" & filename == "KNO3_8h_vs_KCL_8h") %>% pull(genes),
  uptimeBin12h = timeBins %>% filter(type == "upregulated" & filename == "KNO3_12h_vs_KCL_12h") %>% pull(genes),
  dntimeBin2h = timeBins %>% filter(type == "downregulated" & filename == "KNO3_2h_vs_KCL_2h") %>% pull(genes),
  dntimeBin4h = timeBins %>% filter(type == "downregulated" & filename == "KNO3_4h_vs_KCL_4h") %>% pull(genes),
  dntimeBin8h = timeBins %>% filter(type == "downregulated" & filename == "KNO3_8h_vs_KCL_8h") %>% pull(genes)
)

hmap4 <- function(selGene, file_name) {
  selGene <- ifelse(selGene %in% notation$Potra, 
                    paste(selGene, notation$AT_best_hit[match(selGene, notation$Potra)]), selGene)
  vst2 <- allAve[rownames(allAve) %in% selGene, ]
  
  output_dir <- here("data/intGeneHeatmap")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  short_file_name <- substr(file_name, 1, 100)
  output_file <- file.path(output_dir, paste0(short_file_name, ".svg"))
  
  svg(output_file, pointsize = 8)
  
  heatmap_result <- heatmap.2(
    t(scale(t(vst2))), distfun = pearson.dist,
    hclustfun = function(X) hclust(X, method = "ward.D2"),
    trace = "none", col = hpal, margins = c(12, 12), cexCol = 1,
    cexRow = 0.3, main = file_name, key = TRUE, keysize = 1,
    Colv = FALSE, Rowv = TRUE, dendrogram = "row")
  
  gene_order <- rownames(vst2)[rev(heatmap_result$rowInd)]
  write.table(gene_order, file = file.path(output_dir, paste0(short_file_name, "_gene_order.txt")),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  dev.off()
}

lapply(names(gene_set_timeBin), function(name) hmap4(gene_set_timeBin[[name]], name))

# Heatmaps for aspwood data of selective DEGs in suppl figure
aspwoodtpm <- read.table("/mnt/ada/projects/aspseq/htuominen/SNR-results/publisheddatasets/AspWood_tpm.txt",
                         header = TRUE)
aspwoodtpm <- dcast(aspwoodtpm, gene_id ~ sample_name)
orderaspwood <- c("T1-Phloem-01","T1-Phloem-02","T1-Phloem-03","T1-Phloem-04",
                  "T1-Phloem-05","T1-Cambium-06","T1-Cambium-07","T1-Cambium-08",
                  "T1-Cambium-09","T1-Cambium-10","T1-Cambium-11","T1-Cambium-12",
                  "T1-Expanding-xylem-13","T1-Expanding-xylem-14","T1-Expanding-xylem-15",
                  "T1-Expanding-xylem-16","T1-Expanding-xylem-17","T1-Expanding-xylem-18",
                  "T1-Expanding-xylem-19","T1-Lignified-xylem-20","T1-Lignified-xylem-21",
                  "T1-Lignified-xylem-22","T1-Lignified-xylem-23","T1-Lignified-xylem-24",
                  "T1-Lignified-xylem-25","T2-Phloem-01","T2-Phloem-02","T2-Phloem-03",
                  "T2-Phloem-04","T2-Phloem-05","T2-Cambium-06","T2-Cambium-07",
                  "T2-Cambium-08","T2-Cambium-09","T2-Cambium-10","T2-Cambium-11",
                  "T2-Expanding-xylem-12","T2-Expanding-xylem-13","T2-Expanding-xylem-14",
                  "T2-Expanding-xylem-15","T2-Expanding-xylem-16","T2-Expanding-xylem-17",
                  "T2-Expanding-xylem-18","T2-Expanding-xylem-19","T2-Lignified-xylem-20",
                  "T2-Lignified-xylem-21","T2-Lignified-xylem-22","T2-Lignified-xylem-23",
                  "T2-Lignified-xylem-24","T2-Lignified-xylem-25","T2-Lignified-xylem-26",
                  "T3-Phloem-01","T3-Phloem-02","T3-Phloem-03","T3-Phloem-04",
                  "T3-Phloem-05","T3-Cambium-06","T3-Cambium-07","T3-Cambium-08",
                  "T3-Cambium-09","T3-Cambium-10","T3-Cambium-11","T3-Cambium-12",
                  "T3-Cambium-13","T3-Cambium-14","T3-Expanding-xylem-15",
                  "T3-Expanding-xylem-16","T3-Expanding-xylem-17","T3-Expanding-xylem-18",
                  "T3-Expanding-xylem-19","T3-Expanding-xylem-20","T3-Expanding-xylem-21",
                  "T3-Lignified-xylem-22","T3-Lignified-xylem-23","T3-Lignified-xylem-24",
                  "T3-Lignified-xylem-25","T3-Lignified-xylem-26","T3-Lignified-xylem-27",
                  "T3-Lignified-xylem-28","T4-Phloem-01","T4-Phloem-02","T4-Phloem-03",
                  "T4-Phloem-04","T4-Phloem-05","T4-Cambium-06","T4-Cambium-07",
                  "T4-Cambium-08","T4-Cambium-09","T4-Cambium-10","T4-Cambium-11",
                  "T4-Cambium-12","T4-Expanding-xylem-13","T4-Expanding-xylem-14",
                  "T4-Expanding-xylem-15","T4-Expanding-xylem-16","T4-Expanding-xylem-17",
                  "T4-Expanding-xylem-18","T4-Expanding-xylem-19","T4-Expanding-xylem-20",
                  "T4-Lignified-xylem-21","T4-Lignified-xylem-22","T4-Lignified-xylem-23",
                  "T4-Lignified-xylem-24","T4-Lignified-xylem-25","T4-Lignified-xylem-26",
                  "T4-Lignified-xylem-27","T4-Lignified-xylem-28","gene_id")

aspwoodtpm <- aspwoodtpm[orderaspwood]
tree1data <- subset(aspwoodtpm, select = grep("T1-*", colnames(aspwoodtpm)))
rownames(tree1data)<- aspwoodtpm$gene_id

# for all time bins degs
hmap2 <- function(selGene, file_name) {
  tres1 <- tree1data[rownames(tree1data) %in% selGene, ]
  tres1 <- tres1[rowSums(tres1 != 0) > 0, ]
  
  output_dir <- here("data/intGeneHeatmap/")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  short_file_name <- substr(file_name, 1, 100)
  output_file <- file.path(output_dir, paste0(short_file_name, ".svg"))
  order_file <- file.path(output_dir, paste0(short_file_name, "_gene_order.txt"))
  
  heatmap_result <- heatmap.2(
    t(scale(t(tres1))),
    distfun = pearson.dist,
    hclustfun = function(X) hclust(X, method = "ward.D2"),
    trace = "none", col = hpal, margins = c(12, 12),
    cexCol = 1, cexRow = 0.1, key = TRUE, keysize = 1,
    main = file_name, Colv = FALSE, Rowv = TRUE, dendrogram = "row"
  )
  
  gene_order <- rownames(tres1)[rev(heatmap_result$rowInd)]
  write.table(gene_order, file = order_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  # Save the heatmap as SVG
  svg(output_file, width = 10, height = 10, pointsize = 12, family = "Arial")
  heatmap.2(
    t(scale(t(tres1))),
    distfun = pearson.dist,
    hclustfun = function(X) hclust(X, method = "ward.D2"),
    trace = "none", col = hpal, margins = c(12, 12),
    cexCol = 1, cexRow = 0.1, key = TRUE, keysize = 1,
    main = file_name, Colv = FALSE, Rowv = TRUE, dendrogram = "row"
  )
  dev.off()
}
for (name in names(gene_set_timeBin)) {hmap2(gene_set_timeBin[[name]], name)}

# for all degs -all up together and all down together
deg1 <- read_excel("../singleCellTimeSeriesRNASeqManuscript/finalResultsDeseq/miscDegs/intersection_fdr0.05_lfc1.0.xlsx", 
                   sheet = "combined")
gene_sets1 <- list(upregulated = deg1 %>%
                     filter(`level in KNO3` == "Upregulated") %>% select(ID) %>% pull(ID),
                   downregulated = deg1 %>%
                     filter(`level in KNO3` == "Downregulated") %>% select(ID) %>% pull(ID))

hmap1 <- function(selGene, file_name) {
  tres <- tree1data[rownames(tree1data) %in% selGene, ]
  tres1 <- tres[rowSums(tres != 0) > 0, ]
  
  output_dir <- here("data/intGeneHeatmap/")
  if (!dir.exists(output_dir)) {dir.create(output_dir, recursive = TRUE)}
  
  short_file_name <- substr(file_name, 1, 100)
  output_file <- file.path(output_dir, paste0(short_file_name, ".svg"))
  
  svg(output_file, pointsize = 12, family = "Arial", width = 20, height = 20)
  par(family = "Arial", cex = 1)
  
  heatmap_result <- heatmap.2(
    t(scale(t(tres1))),
    distfun = pearson.dist,
    hclustfun = function(X) hclust(X, method = "ward.D2"),
    trace = "none", col = hpal, margins = c(12, 12), cexCol = 1,
    cexRow = 0.5, key = TRUE, keysize = 1, main = file_name,
    Colv = FALSE, Rowv = TRUE, dendrogram = "row"
  )
  
  gene_order <- rownames(tres1)[rev(heatmap_result$rowInd)]
  write.table(gene_order, file = file.path(output_dir, paste0(short_file_name, "_gene_order.txt")),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  dev.off()
}

for (name in names(gene_set_timeBin)) {
  # Open the SVG device with appropriate size
  svg_file <- here("data/intGeneHeatmap", paste0(substr(name, 1, 100), ".svg"))
  svg(svg_file, width = 10, height = 10, pointsize = 12, family = "Arial")
  
  # Call the heatmap function
  hmap2(gene_set_timeBin[[name]], name)
  
  # Close the device to ensure each plot is saved correctly
  dev.off()
  
  # Reset the plotting parameters (optional but recommended)
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1)
}

subset_up <- tree1data[rownames(tree1data) %in% gene_sets1$upregulated, ]
subset_dn <- tree1data[rownames(tree1data) %in% gene_sets1$downregulated, ]
highest_expression_tissue <- apply(subset_dn, 1, function(x) names(x)[which.max(x)])
result <- data.frame(Gene = rownames(subset_dn), Highest_Expression_Tissue = highest_expression_tissue)
write.csv(result, "highest_expression_tissue.csv", row.names = FALSE, quote = F)

# for all degs -all up together and all down together time wise
deg3 <- read.csv("doc/deg_KNOvsKCL_fdr0.05_lfc0.5.csv")
gene_sets3 <- list(upAll2h = deg3 %>%
                     filter(timepoint == "KNO3_2h_vs_KCL_2h" & level.in.KNO3 == "Upregulated") %>%
                     select(X) %>% pull(X),
                   upAll4h = deg3 %>%
                     filter(timepoint == "KNO3_4h_vs_KCL_4h" & level.in.KNO3 == "Upregulated") %>%
                     select(X) %>% pull(X),
                   upAll8h = deg3 %>%
                     filter(timepoint == "KNO3_8h_vs_KCL_8h" & level.in.KNO3 == "Upregulated") %>%
                     select(X) %>% pull(X),
                   upAll12h = deg3 %>%
                     filter(timepoint == "KNO3_12h_vs_KCL_12h" & level.in.KNO3 == "Upregulated") %>%
                     select(X) %>% pull(X),
                   dnAll2h = deg3 %>%
                     filter(timepoint == "KNO3_2h_vs_KCL_2h" & level.in.KNO3 == "Downregulated") %>%
                     select(X) %>% pull(X),
                   dnAll4h = deg3 %>%
                     filter(timepoint == "KNO3_4h_vs_KCL_4h" & level.in.KNO3 == "Downregulated") %>%
                     select(X) %>% pull(X),
                   dnAll8h = deg3 %>%
                     filter(timepoint == "KNO3_8h_vs_KCL_8h" & level.in.KNO3 == "Downregulated") %>%
                     select(X) %>% pull(X))

for (name in names(gene_sets3)) {hmap1(gene_sets3[[name]], name)}

#' # Session Info 
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```