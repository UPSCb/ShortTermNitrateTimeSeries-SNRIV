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
suppressPackageStartupMessages({
    library(data.table)
    library(DESeq2)
    library(gplots)
    # library(ggplots2)
    library(here)
    library(hyperSpec)
    library(RColorBrewer)
    library(tidyverse)
    library(VennDiagram)
})

#' * Helper files
suppressMessages({
    # source(here("UPSCb-common/templates/plotEnrichedTreemap.R"))
    source(here("UPSCb-common/src/R/featureSelection.R"))
    source(here("UPSCb-common/src/R/volcanoPlot.R"))
    # souxrce(here("UPSCb-common/src/R/gopher.R"))
})

#' * Graphics
pal=brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' * Data
#' ```{r load, echo=FALSE,eval=FALSE}
#' CHANGEME - here you are meant to load an RData object
#' that contains a DESeqDataSet object. If you ran the 
#' biological QA template, you need not change anything
#' ```
load(here("data/analysis/salmon/dds.rds"))

#' ## Normalisation for visualisation
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)
dir.create(here("data/analysis/DE"),showWarnings=FALSE)
save(vst,file=here("data/analysis/DE/vst-aware.rda"))
write_delim(as.data.frame(vst) %>% rownames_to_column("ID"),
            here("data/analysis/DE/vst-aware.tsv"))

#' ## Gene of interests
#' ```{r goi, echo=FALSE,eval=FALSE}
#' CHANGEME - Here, you can plot the expression pattern of your gene of
#' interest. You need to have a list of genes in a text file, one geneID per line
#' The ID should exist in your vst data.
#' Note that for the plot to work, you also need to edit the first function (`line_plot`)
#' at the top of this file
#' ```
goi <- read_lines(here("doc/goi.txt"))
stopifnot(all(goi %in% rownames(vst)))
dev.null <- lapply(goi,line_plot,dds=dds,vst=vst)

#' ## Differential Expression
dds <- DESeq(dds)

#' * Dispersion estimation
#' The dispersion estimation is adequate
plotDispEsts(dds)

#' Check the different contrasts
resultsNames(dds)
dds$treatmentTime <- relevel(dds$treatmentTime, "None_0h")
dds <- DESeq(dds)
saveRDS(dds,file=here("data/dds2.rds"))

#' ## Results
#' ## Differential Expression

#'  1. design: dds: effect of Treatment on the Time effect
#'  = ~ Treatment + Time+ Treatment: Time

#'  2. design1: dds1: effect of Treatment:  = ~ Treatment

#'  3. design2: dds2: effect of Time:  = ~ Time

#'  4. design3: dds3: effect of Time, while regressing out variation due 
#'  to Treatment = ~ Treatment + Time
#'  dds 3 means that there is no interaction- i.e, the effect of treatment is
#'  the same in all time points

#'  5. design4: dds4: effect of Time on Treatment effects (same as dds)
#'  = ~ Time + Treatment + Time:Treatment

#'  6. design5: dds5: effect of Time, while regressing out variation due 
#'  to Treatment = ~ Time + Treatment and then use groups
#'
#' design = ~ Treatment + Time+ Treatment: Time) #dds
#' design = ~ Treatment) #dds1
#' design = ~ Time) #dds2
#' design = ~ Treatment + Time) #dds3
#' design = ~ Time + Treatment + Time:Treatment) #dds4
#' design = ~ Time + Treatment) #dds5
#'
#' 0h 2h 4h 8h 12h 24h 48h
#' 1. KCL/None dds, KNO3/None dds, design ~time for both
#' 2. KCL/KNO3 dds, design ~time*treatment- effect of time is removed
#' 3. KCL/KNO3/None/None dds (change treatment name so None doesn't exist anymore),
#' design ~time*treatment
#' 
#' For KCL versus KNO3 dds excluding 0 time, design ~time*treatment
#' dds was generated like so:
# dds <- DESeqDataSetFromTximport (txi=txi, colData = samples, 
#                                  design =~Time*Treatment)
# dds$Time <- relevel(dds$Time,ref = "2h")
# dds <- DESeq(dds)
# Time4h.TreatmentKNO3= all DEGs which are not in other time points

#' the second way is by using LRT:
# full_model <- ~ Time + Treatment + Treatment:Time
# For LRT test, provide a reduced model, that is the full model without
# treatment:time term:
# reduced_model <- ~ Treatment + Time
# generate dds like the following:
# dds1 <- DESeqDataSetFromTximport(txi =txi, colData = samples,
#                                  design = ~ Treatment + Time + Treatment:Time)
# dds_lrt_time <- DESeq(dds1, test="LRT", reduced = ~ Treatment + Time)
# 
# clusters <- degPatterns(cluster_rlog, metadata = samples, time="Time", 
#                         col="Treatment")

# From Edoardo:
# Define your contrasts as a list of vectors
dds <- readRDS("~/shruti/ShortTermNitrateTimeSeries-SNRIV/data/dds2.rds")
samples <- read.table(here("doc/sampleTremula.csv"), header = T)
colnames(dds) <- paste(samples$sample, samples$treatmentTime,
                       samples$replicate, sep = "_")
load("~/shruti/ShortTermNitrateTimeSeries-SNRIV/data/analysis/DE/vst-aware.rda")

# Define contrasts
contrast2h = c("treatmentTime", "KNO3_2h", "KCL_2h", "None_0h")
contrast4h = c("treatmentTime", "KNO3_4h", "KCL_4h", "None_0h")
contrast8h = c("treatmentTime", "KNO3_8h", "KCL_8h", "None_0h")
contrast12h = c("treatmentTime", "KNO3_12h", "KCL_12h", "None_0h")
contrast24h = c("treatmentTime", "KNO3_24h", "KCL_24h", "None_0h")
contrast48h = c("treatmentTime", "KNO3_48h", "KCL_48h", "None_0h")

# Initialize lists to store results
results_list <- list()
bins_up <- list()
bins_dn <- list()

# Define the function to process each contrast group
process_contrast <- function(contrast, dds, vst, bins_up, bins_dn, results_list) {
  description <- contrast[1]
  kno_timepoint <- contrast[2]
  kcl_timepoint <- contrast[3]
  baseline <- contrast[4]
  
  # Generate contrast names for file naming
  contrast_name1 <- paste(description, kno_timepoint, "vs", kcl_timepoint, sep = "_")
  contrast_1 <- c(description, kno_timepoint, kcl_timepoint)
  
  contrast_name2 <- paste(description, kno_timepoint, "vs", baseline, sep = "_")
  contrast_2 <- c(description, kno_timepoint, baseline)
  
  contrast_name3 <- paste(description, kcl_timepoint, "vs", baseline, sep = "_")
  contrast_3 <- c(description, kcl_timepoint, baseline)
  
  # Extract results for KNO3_Xh vs KCL_Xh
  results_KNO3_vs_KCL <- extract_results(dds, vst, contrast_1, plot = FALSE, default_prefix = contrast_name1)
  Upregulated_KNO3_vs_KCL <- results_KNO3_vs_KCL$up
  Downregulated_KNO3_vs_KCL <- results_KNO3_vs_KCL$dn
  
  # Extract results for KNO3_Xh vs None_0h
  results_KNO3_vs_0h <- extract_results(dds, vst, contrast_2, plot = FALSE, default_prefix = contrast_name2)
  Upregulated_KNO3_vs_0h <- results_KNO3_vs_0h$up
  Downregulated_KNO3_vs_0h <- results_KNO3_vs_0h$dn
  
  # Extract results for KCL_Xh vs None_0h
  results_KCL_vs_0h <- extract_results(dds, vst, contrast_3, plot = FALSE, default_prefix = contrast_name3)
  Upregulated_KCL_vs_0h <- results_KCL_vs_0h$up
  Downregulated_KCL_vs_0h <- results_KCL_vs_0h$dn
  
  # Find intersection and unique genes
  common_upreg <- intersect(Upregulated_KNO3_vs_0h, Upregulated_KNO3_vs_KCL)
  true_upreg <- setdiff(common_upreg, Upregulated_KCL_vs_0h)
  intersection_all_three_upreg <- Reduce(intersect, list(Upregulated_KNO3_vs_0h, Upregulated_KNO3_vs_KCL, Upregulated_KCL_vs_0h))
  final_upreg <- union(true_upreg, intersection_all_three_upreg)
  
  common_dnreg <- intersect(Downregulated_KNO3_vs_0h, Downregulated_KNO3_vs_KCL)
  true_dnreg <- setdiff(common_dnreg, Downregulated_KCL_vs_0h)
  intersection_all_three_dnreg <- Reduce(intersect, list(Downregulated_KNO3_vs_0h, Downregulated_KNO3_vs_KCL, Downregulated_KCL_vs_0h))
  final_dnreg <- union(true_dnreg, intersection_all_three_dnreg)
  
  # Find new genes upregulated and downregulated at KNO3_Xh for the first time
  all_genes_already_in_bins_up <- unlist(bins_up)
  new_genes_up <- setdiff(final_upreg, all_genes_already_in_bins_up)
  bins_up[[contrast_name1]] <- new_genes_up
  
  all_genes_already_in_bins_dn <- unlist(bins_dn)
  new_genes_down <- setdiff(final_dnreg, all_genes_already_in_bins_dn)
  bins_dn[[contrast_name1]] <- new_genes_down
  
  # Store results in results_list
  # results_list[[contrast_name1]] <- list(
  #   upregulated = common_upreg,
  #   downregulated = common_dnreg
  # )
  
  results_list[[kno_timepoint]] <- list(
    upregulated = common_upreg,
    downregulated = common_dnreg
  )
  
  # Save intersection genes to CSV
  # intersection_genes <- data.frame(upregulated = I(list(common_upreg)), downregulated = I(list(common_dnreg)))
  # write.csv(intersection_genes, file = paste0("data/analysis/DE/", contrast_name1, "_intersection_degs.csv"), row.names = FALSE)
  # 
  true_degs <- data.frame(upregulated = I(list(final_upreg)), downregulated = I(list(common_upreg)))
  write.csv(true_degs, file = paste0("data/analysis/DE/", contrast_name1, "_true_degs.csv"), row.names = FALSE)
  
  # Handle empty vectors and save bins_up and bins_dn to CSV
  if (length(new_genes_up) > 0) {
    bins_up_df <- data.frame(genes = new_genes_up, type = "upregulated")
  } else {
    bins_up_df <- data.frame(genes = character(0), type = character(0))
  }
  
  if (length(new_genes_down) > 0) {
    bins_dn_df <- data.frame(genes = new_genes_down, type = "downregulated")
  } else {
    bins_dn_df <- data.frame(genes = character(0), type = character(0))
  }
  
  bins_combined <- rbind(bins_up_df, bins_dn_df)
  write.csv(bins_combined, file = paste0("data/analysis/DE/", contrast_name1, "_bins.csv"), row.names = FALSE)
  
  return(list(results_list = results_list, bins_up = bins_up, bins_dn = bins_dn))
}

# List of contrasts
contrasts <- list(contrast2h, contrast4h, contrast8h, contrast12h, contrast24h, 
                  contrast48h)

# Initialize lists to store results
results_list <- list()
bins_up <- list()
bins_dn <- list()

# Loop through each contrast and process
for (contrast in contrasts) {
  res <- process_contrast(contrast, dds, vst, bins_up, bins_dn, results_list)
  results_list <- res$results_list
  bins_up <- res$bins_up
  bins_dn <- res$bins_dn
}

# Optionally, you can print or inspect the results_list, bins_up, and bins_dn
print(results_list)
print(bins_up)
print(bins_dn)

  # Optionally, you can store the results for further processing
  # results_list[[paste0(contrast[2], "_vs_", contrast[3])]] <- result

# If you want to clean vst file of zero rows
vst_matrix_clean <- vst[rowSums(is.na(vst) | vst == 0) < ncol(vst), ]
nrow(vst_matrix_clean)
save(vst_matrix_clean,file=here("data/analysis/DE/vst-noZero.rda"))

# If you want to filter genes from vst with NA values
library(dplyr)
file_list <- list.files(path = "data/analysis/DE", pattern = "*results.csv", full.names = TRUE)
gene_list <- list()

for (file in file_list) {
  data <- read.csv(file_list[1], stringsAsFactors = FALSE)
  
  if ("pvalue" %in% colnames(data)) {
    filtered_genes <- data %>% filter(!is.na(pvalue)) %>% select(X) 
    gene_list <- append(gene_list, filtered_genes$X)
  }
}

unique_genes <- unique(unlist(gene_list))
unique_genes <- as.character(unique_genes)

vst_subset <- vst[rownames(vst) %in% unique_genes, , drop = FALSE]
nrow(vst_subset)
save(vst_subset,file=here("data/analysis/DE/vst-noNA.rda"))

#' # Session Info 
#'  ```{r session info, echo=FALSE}
#'  sessionInfo()
#'  ```


