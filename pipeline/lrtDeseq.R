suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(gplots)
  library(here)
  library(hyperSpec)
  # library(RColorBrewer)
  library(tidyverse)
  # library(gplots)
  # library(DEGreport)
  # library(SummarizedExperiment)
  library(magrittr)
  library(parallel)
  library(plotly)
  library(pvclust)
  library(tximport)
  # library(pheatmap)
  library(readr)
  library(vsn)
})

extract_DEGs <- function(result_name, dds_object, alpha = 0.05) {
  results_object <- results(dds_object, name = result_name)
  sel_effect <- results_object$padj <= alpha & !is.na(results_object$padj)
  DEGs <- data.frame(
    Gene = rownames(results_object)[sel_effect],
    # baseMean = results_object$baseMean[sel_effect],
    log2FoldChange = results_object$log2FoldChange[sel_effect],
    # lfcSE = results_object$lfcSE[sel_effect],
    # stat = results_object$stat[sel_effect],
    padj = results_object$padj[sel_effect]
  )
  return(DEGs)
}
output_dir <- "doc/"

# https://github.com/hbctraining/DGE_workshop_salmon_online/blob/master/lessons/08a_DGE_LRT_results.md  
 
# ## Exploring results from the Likelihood ratio test (LRT)
# DESeq2 also offers LRT as an alternative when evaluating expression change 
# across more than two levels. significant genes are the ones that change in
# expression in any direction across the different factor levels.

# Generally, this test will result in a larger number of genes than the 
# individual pair-wise comparisons. While the LRT is a test of significance
# for differences of any level(s) of the factor, one should not expect it to be 
# exactly equal to the union of sets of genes using Wald tests (although we do 
# expect a high degree of overlap).

# run QC with desqQC.R script

# keep the main effect in the last position of your design formula which will be 
# default for results() 
# dds <- DESeqDataSetFromTximport(txi =txi, colData = samples,
#                                  design = ~ Treatment + Time + Treatment:Time)

# suggested by Edoardo: for this to work, modify the sample file first
# to exclude or duplicate the metadata with 0h or you will get an error
# dds1 <- DESeqDataSetFromTximport(txi=txi,
#                                  colData = samples, design =~Time*Treatment)
# dds_lrt_time <- DESeq(dds, test="LRT", reduced = ~Time*Treatment)

# Case1a: exclude genes at zero h
# dds_lrt_time <- DESeq(dds, test="LRT", reduced = ~ Treatment + Time)

# Case1b: include genes which are DE at time 0- don't know what it means here
# dds_lrt_timeOnly <- DESeq(dds, test="LRT", reduced = ~ Time)

# Case1c: include genes which are DE at treatment 0- don't know what it means here
# dds_lrt_TreatmentOnly <- DESeq(dds, test="LRT", reduced = ~ Treatment)

# no need for contrasts since we are not making a pair-wise comparison.

### Why are fold changes reported for an LRT test?
# For analyses using the likelihood ratio test, the p-values are determined solely
# by the difference in deviance between the full and reduced model formula.
# A single log2 fold change is printed in the results table for consistency 
# with other results table outputs, but is not associated with the actual test.**

# Columns relevant to the LRT test:**
# baseMean`: mean of normalized counts for all samples
# stat`: the difference in deviance between reduced and full model
# pvalue`: stat value is compared to a chi-squared distribution to generate a pvalue
# padj`: BH adjusted p-values
# log2FoldChange`: log2 fold change
# lfcSE`: standard error

# use the `name` argument the results table to get these results for two groups
# like in LRT, p-value and Lfc will vary for each, but all DEGs remain same

# Define time points and corresponding result names
# time_points <- c("0h", "2h", "4h", "8h", "12h", "24h", "48h")
# treatment <- c("KNO3", "KCL", "None")

# For tremula
# case 4.1: combine treatment and time in to one variable
samples <- read_table(here("doc/sampleTremula.csv"))
tx2gene <- suppressMessages(read_delim(here("../single_cell_analysis_poplar/reference/annotation/tx2gene.tsv.gz"),
                                       delim="\t", col_names=c("TXID","GENE")))
txi <- suppressMessages(tximport(files = samples$file,type = "salmon",
                                 tx2gene=tx2gene))
dds2 <- DESeqDataSetFromTximport(txi=txi, colData = samples, 
                                 design =~treatmentTime)
dds2 <- DESeq(dds2)

# correct the level with 0h as the base, 
# there is no need to do a separate relevel. Just give the contrasts
# i.e., the dds2$treatmentTime <- relevel(dds2$treatmentTime, "KCL_8h")
# will give the same reslts for KCL_8h vs KNO3_8h as in the relevel with 0h
dds2$treatmentTime <- relevel(dds2$treatmentTime, "None_0h")
dds2 <- DESeq(dds2)
saveRDS(dds2,file=here("data/dds2.rds"))

# counts <- txi$counts
# colnames(counts) <- samples$sample

# 4.1.a get DEGs of each treatment.time versus Time 0h
dds2 <- readRDS("~/shruti/ShortTermNitrateTimeSeries-SNRIV/data/dds2.rds")
result_names <- resultsNames(dds2)
DEGs_list <- list()

for (i in seq_along(result_names)) {
  result_name <- result_names[i]
  DEGs <- extract_DEGs(result_name, dds2)
  DEGs_list[[result_name]] <- DEGs
  output_file <- file.path(output_dir, paste0(result_name, ".csv"))
  write.csv(DEGs, file = output_file, row.names = FALSE)
}

# 4.1.b. between the treatments at specific time points
kno_kcl_2h <- results(dds2, contrast = c("treatmentTime","KNO3_2h", "KCL_2h"))
kno_kcl_4h <- results(dds2, contrast = c("treatmentTime","KNO3_4h", "KCL_4h"))
kno_kcl_8h <- results(dds2, contrast = c("treatmentTime","KNO3_8h", "KCL_8h"))
kno_kcl_12h <- results(dds2, contrast = c("treatmentTime","KNO3_12h", "KCL_12h"))
kno_kcl_24h <- results(dds2, contrast = c("treatmentTime","KNO3_24h", "KCL_24h"))
kno_kcl_48h <- results(dds2, contrast = c("treatmentTime","KNO3_48h", "KCL_48h"))

sel2h <- kno_kcl_2h$padj <= 0.05 & ! is.na(kno_kcl_2h$padj) 
sel4h <- kno_kcl_4h$padj <= 0.05 & ! is.na(kno_kcl_4h$padj)
sel8h <- kno_kcl_8h$padj <= 0.05 & ! is.na(kno_kcl_8h$padj)
sel12h <- kno_kcl_12h$padj <= 0.05 & ! is.na(kno_kcl_12h$padj)
sel24h <- kno_kcl_24h$padj <= 0.05 & ! is.na(kno_kcl_24h$padj)
sel48h <- kno_kcl_48h$padj <= 0.05 & ! is.na(kno_kcl_48h$padj)

write.csv(kno_kcl_2h[sel2h,],file = here("doc/KNO2h_vs_KCL2h.csv"), quote = F)
write.csv(kno_kcl_4h[sel4h,],file = here("doc/KNO4h_vs_KCL4h.csv"), quote = F)
write.csv(kno_kcl_8h[sel8h,],file = here("doc/KNO8h_vs_KCL8h.csv"), quote = F)
write.csv(kno_kcl_12h[sel12h,],file = here("doc/KNO12h_vs_KCL12h.csv"),quote = F)
write.csv(kno_kcl_24h[sel24h,],file = here("doc/KNO24h_vs_KCL24h.csv"), quote = F)
write.csv(kno_kcl_48h[sel48h,],file = here("doc/KNO48h_vs_KCL48h.csv"), quote = F)

# Idea 1: Firstly, filtered the genes based on FDR<0.01 and lfc>1
# To find the true Difference After 2h Of KNO3 treatment in plants, 
# pick DEGs that belong to and are unique to kno2h_vs_0h among intersection of 
# kno2h_vs_kcl_2h, kno2h_vs_0h and kcl2h_vs_0h. This is done to discard the
# DEGs which might be there due to influence of kcl2h treatment. In this, DEGs
# due to time or treatment alone have not been considered.

# Idea2: if I want to see DEGs that are different in 2h vs 4h or between 
# different timepoints, I should try another approach

# 4.1.c. with treatment only
dds2 <- DESeqDataSetFromTximport(txi=txi, colData = samples, design =~treatment)
dds2 <- DESeq(dds2)
dds2$treatment <- relevel(dds2$treatment, "None")
dds2 <- DESeq(dds2)

result_names <- resultsNames(dds2)
DEGs_list <- list()

for (i in seq_along(result_names)) {
  result_name <- result_names[i]
  DEGs <- extract_DEGs(result_name, dds2)
  DEGs_list[[result_name]] <- DEGs
  output_file <- file.path(output_dir, paste0(result_name, ".csv"))
  write.csv(DEGs, file = output_file, row.names = FALSE)
}

# 4.1.d. with time only
dds2 <- DESeqDataSetFromTximport(txi=txi, colData = samples, design =~time)
dds2 <- DESeq(dds2)
dds2$time <- relevel(dds2$time, "0h")
dds2 <- DESeq(dds2)

result_names <- resultsNames(dds2)
DEGs_list <- list()

for (i in seq_along(result_names)) {
  result_name <- result_names[i]
  DEGs <- extract_DEGs(result_name, dds2)
  DEGs_list[[result_name]] <- DEGs
  output_file <- file.path(output_dir, paste0(result_name, ".csv"))
  write.csv(DEGs, file = output_file, row.names = FALSE)
}

# 4.2. remove 0h sample
samples <- read_table(here("doc/sampleTremula.csv"))
tx2gene <- suppressMessages(read_delim(here("../single_cell_analysis_poplar/reference/annotation/tx2gene.tsv.gz"),
                                       delim="\t", col_names=c("TXID","GENE")))
txi <- suppressMessages(tximport(files = samples$file, type = "salmon",
                                 tx2gene=tx2gene))

dds <- DESeqDataSetFromTximport(txi =txi, colData = samples,
                                 design = ~ treatment + time + treatment:time)

# 4.2.a: exclude genes at zero h- 
# gives the difference in the effect of treatment over time (treatment:time)
# LRT- we remove the treatment-specific differences over time. 
# Genes with small p values from this test are those which at one or more time
# points after time 0 showed a treatment-specific effect. Note therefore that 
# this will not give small p values to genes that moved up or down over time 
# in the same way in both treatments
dds_lrt_treattime <- DESeq(dds, test="LRT", reduced = ~ treatment + time)

dds_lrt_treattime$time <- relevel(dds_lrt_treattime$time,"2h")
dds_lrt_treattime <- DESeq(dds_lrt_treattime)

kno_kcl <- results (dds_lrt_treattime, name = "treatment_KNO3_vs_KCL")
selknkcl <- kno_kcl$padj <= 0.05 & ! is.na(kno_kcl$padj)
write.table(kno_kcl[selknkcl,],file = here("doc/kno_kcl.txt"), sep="\t",
            col.names=T, row.names=T, quote=F)

# 4.2.b: include genes which are DE at time 0h
dds_lrt_timeOnly <- DESeq(dds, test="LRT", reduced = ~ time)
dds_lrt_timeOnly$time <- relevel(dds_lrt_timeOnly$time,"2h")

dds_lrt_timeOnly$time <- relevel(dds_lrt_timeOnly$time,"2h")
dds_lrt_timeOnly <- DESeq(dds_lrt_timeOnly)

kno_kcl <- results (dds_lrt_timeOnly, name = "treatment_KNO3_vs_KCL")
selknkcl <- kno_kcl$padj <= 0.05 & ! is.na(kno_kcl$padj)
write.table(kno_kcl[selknkcl,],file = here("doc/kno_kcl.txt"), sep="\t",
            col.names=T, row.names=T, quote=F)

# 4.2.c: include genes which are DE at time 0- don't know what it means here
dds_lrt_TreatmentOnly <- DESeq(dds, test="LRT", reduced = ~ treatment)

dds_lrt_TreatmentOnly$time <- relevel(dds_lrt_TreatmentOnly$time, "2h")
dds_lrt_TreatmentOnly  <- DESeq(dds_lrt_TreatmentOnly)

kno_kcl <- results (dds_lrt_TreatmentOnly, name = "treatmentKNO3.time4h")
selknkcl <- kno_kcl$padj <= 0.05 & ! is.na(kno_kcl$padj)
write.table(kno_kcl[selknkcl,],file = here("doc/kno_kcl.txt"), sep="\t",
            col.names=T, row.names=T, quote=F)


# 4.3: samples at zero h as KNO3- 
samples <- read_table(here("doc/knoNone.csv"))
tx2gene <- suppressMessages(read_delim(here("../single_cell_analysis_poplar/reference/annotation/tx2gene.tsv.gz"),
                                       delim="\t", col_names=c("TXID","GENE")))
txi <- suppressMessages(tximport(files = samples$file, type = "salmon",
                                 tx2gene=tx2gene))

dds <- DESeqDataSetFromTximport(txi =txi, colData = samples, design = ~ time)
dds <- DESeq(dds)
dds$time <- relevel(dds$time,"0h")
dds <- DESeq(dds)

result_names <- resultsNames(dds)

DEGs_list <- list()

for (i in seq_along(result_names)) {
  result_name <- result_names[i]
  DEGs <- extract_DEGs(result_name, dds)
  DEGs_list[[result_name]] <- DEGs
  output_file <- file.path(output_dir, paste0(result_name, ".csv"))
  write.csv(DEGs, file = output_file, row.names = FALSE)
}

# 4.4: samples at zero h as KCL- 
samples <- read_table(here("doc/kclNone.csv"))
tx2gene <- suppressMessages(read_delim(here("../single_cell_analysis_poplar/reference/annotation/tx2gene.tsv.gz"),
                                       delim="\t", col_names=c("TXID","GENE")))
txi <- suppressMessages(tximport(files = samples$file, type = "salmon",
                                 tx2gene=tx2gene))

dds <- DESeqDataSetFromTximport(txi =txi, colData = samples, design = ~ time)
dds <- DESeq(dds)
dds$time <- relevel(dds$time,"0h")
dds <- DESeq(dds)

result_names <- resultsNames(dds)

DEGs_list <- list()

for (i in seq_along(result_names)) {
  result_name <- result_names[i]
  DEGs <- extract_DEGs(result_name, dds)
  DEGs_list[[result_name]] <- DEGs
  output_file <- file.path(output_dir, paste0(result_name, ".csv"))
  write.csv(DEGs, file = output_file, row.names = FALSE)
}

# 4.5: duplicate samples for 0h into KCL and KNO3
samples <- read_table(here("doc/kclkno0hDupl.csv"))
tx2gene <- suppressMessages(read_delim(here("../single_cell_analysis_poplar/reference/annotation/tx2gene.tsv.gz"),
                                       delim="\t", col_names=c("TXID","GENE")))
txi <- suppressMessages(tximport(files = samples$file, type = "salmon",
                                 tx2gene=tx2gene))

dds <- DESeqDataSetFromTximport(txi =txi, colData = samples,
                                design = ~ treatment + time + treatment:time)

# 4.5.a: exclude genes at zero h- 
# gives the difference in the effect of treatment over time (treatment:time)
# LRT- we remove the treatment-specific differences over time. 
# Genes with small p values from this test are those which at one or more time
# points after time 0 showed a treatment-specific effect. Note therefore that 
# this will not give small p values to genes that moved up or down over time 
# in the same way in both treatments
dds_2 <- DESeq(dds, test="LRT", reduced = ~ treatment + time)
saveRDS(dds_lrt_treattime,file=here("data/dds1.rds"))
result_names <- resultsNames(dds_lrt_treattime)
DEGs_list <- list()
for (i in seq_along(result_names)) {
  result_name <- result_names[i]
  DEGs <- extract_DEGs(result_name, dds_lrt_treattime)
  DEGs_list[[result_name]] <- DEGs
  output_file <- file.path(output_dir, paste0(result_name, ".csv"))
  write.csv(DEGs, file = output_file, row.names = FALSE)
}

# 4.5.b: include genes which are DE at time 0h
dds_lrt_timeOnly <- DESeq(dds, test="LRT", reduced = ~ time)
saveRDS(dds_lrt_timeOnly,file=here("data/dds1.rds"))
result_names <- resultsNames(dds_lrt_timeOnly)
DEGs_list <- list()
for (i in seq_along(result_names)) {
  result_name <- result_names[i]
  DEGs <- extract_DEGs(result_name, dds_lrt_timeOnly)
  DEGs_list[[result_name]] <- DEGs
  output_file <- file.path(output_dir, paste0(result_name, ".csv"))
  write.csv(DEGs, file = output_file, row.names = FALSE)
}

# 4.5.c: include genes which are DE at time 0- don't know what it means here
dds_lrt_TreatmentOnly <- DESeq(dds, test="LRT", reduced = ~ treatment)
saveRDS(dds_lrt_TreatmentOnly,file=here("data/dds1.rds"))
result_names <- resultsNames(dds_lrt_TreatmentOnly)
DEGs_list <- list()
for (i in seq_along(result_names)) {
  result_name <- result_names[i]
  DEGs <- extract_DEGs(result_name, dds_lrt_TreatmentOnly)
  DEGs_list[[result_name]] <- DEGs
  output_file <- file.path(output_dir, paste0(result_name, ".csv"))
  write.csv(DEGs, file = output_file, row.names = FALSE)
}

# Go to Mfuzz.R for clustering and GO