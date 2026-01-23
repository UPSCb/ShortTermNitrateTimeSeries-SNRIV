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
#' 
suppressPackageStartupMessages({
    library(data.table)
    library(DESeq2)
    library(gplots)
    library(ggplots2)
    library(here)
    library(hyperSpec)
    library(RColorBrewer)
    library(tidyverse)
    library(VennDiagram)
})
#'
#' * Helper files
suppressMessages({
    source(here("UPSCb-common/src/R/featureSelection.R"))
    source(here("UPSCb-common/src/R/volcanoPlot.R"))
})
#'
#'#' * Graphics
pal=brewer.pal(8,"Dark2")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' * Functions
#' 1. plot specific gene expression
line_plot <- function(dds=dds,vst=vst,gene_id=gene_id){
  message(paste("Plotting",gene_id))
  sel <- grepl(gene_id,rownames(vst))
  stopifnot(sum(sel)==1)
  
  p <- ggplot(bind_cols(as.data.frame(colData(dds)),
                        data.frame(value=vst[sel,])),
              aes(x=sample,y=value,col=treatment,group=time)) +
    geom_point() + geom_smooth() +
    scale_y_continuous(name="VST expression") + 
    ggtitle(label=paste("Expression for: ",gene_id))
  
  suppressMessages(suppressWarnings(plot(p)))
  return(NULL)
}

#' 2. extract the DE results. Default cutoffs are
#' from Schurch _et al._, RNA, 2016
extract_results <- function(dds,vst,contrast,padj=0.05,lfc=1, plot=TRUE,
                            verbose=TRUE, export=TRUE,
                            default_dir=here("data/analysis/DE/"),
                            default_prefix="DE-",
                            labels=colnames(dds), sample_sel=1:ncol(dds),
                            expression_cutoff=0, cexCol=.9, 
                            debug=FALSE,filter=c("median",NULL),...){
  
  # get the filter
  if(!is.null(match.arg(filter))){
    filter <- rowMedians(counts(dds,normalized=TRUE))
    message("Using the median normalized counts as default, set filter=NULL to revert to using the mean")
  }
  
  # validation
  if(length(contrast)==1){
    res <- results(dds,name=contrast,filter = filter)
  } else {
    res <- results(dds,contrast=contrast,filter = filter)
  }
  
  stopifnot(length(sample_sel)==ncol(vst))
  
  if(plot){
    par(mar=c(5,5,5,5))
    volcanoPlot(res)
    par(mar=mar)
  }
  
  # a look at independent filtering
  if(plot){
    plot(metadata(res)$filterNumRej,
         type="b", ylab="number of rejections",
         xlab="quantiles of filter")
    lines(metadata(res)$lo.fit, col="red")
    abline(v=metadata(res)$filterTheta)
  }
  
  if(verbose){
    message(sprintf("The independent filtering cutoff is %s, removing %s of the data",
                    round(metadata(res)$filterThreshold,digits=5),
                    names(metadata(res)$filterThreshold)))
    
    max.theta <- metadata(res)$filterNumRej[which.max(metadata(res)$filterNumRej$numRej),"theta"]
    message(sprintf("The independent filtering maximises for %s %% of the data, corresponding to a base mean expression of %s (library-size normalised read)",
                    round(max.theta*100,digits=5),
                    round(quantile(counts(dds,normalized=TRUE),probs=max.theta),digits=5)))
  }
  
  if(plot){
    qtl.exp=quantile(counts(dds,normalized=TRUE),probs=metadata(res)$filterNumRej$theta)
    dat <- data.frame(thetas=metadata(res)$filterNumRej$theta,
                      qtl.exp=qtl.exp,
                      number.degs=sapply(lapply(qtl.exp,function(qe){
                        res$padj <= padj & abs(res$log2FoldChange) >= lfc & 
                          ! is.na(res$padj) & res$baseMean >= qe
                      }),sum))
    if(debug){
      plot(ggplot(dat,aes(x=thetas,y=qtl.exp)) + 
             geom_line() + geom_point() +
             scale_x_continuous("quantiles of expression") + 
             scale_y_continuous("base mean expression") +
             geom_hline(yintercept=expression_cutoff,
                        linetype="dotted",col="red"))
      
      p <- ggplot(dat,aes(x=thetas,y=qtl.exp)) + 
        geom_line() + geom_point() +
        scale_x_continuous("quantiles of expression") + 
        scale_y_log10("base mean expression") + 
        geom_hline(yintercept=expression_cutoff,
                   linetype="dotted",col="red")
      suppressMessages(suppressWarnings(plot(p)))
      
      plot(ggplot(dat,aes(x=thetas,y=number.degs)) + 
             geom_line() + geom_point() +
             geom_hline(yintercept=dat$number.degs[1],linetype="dashed") +
             scale_x_continuous("quantiles of expression") + 
             scale_y_continuous("Number of DE genes"))
      
      plot(ggplot(dat,aes(x=thetas,y=number.degs[1] - number.degs),aes()) + 
             geom_line() + geom_point() +
             scale_x_continuous("quantiles of expression") + 
             scale_y_continuous("Cumulative number of DE genes"))
      
      plot(ggplot(data.frame(x=dat$thetas[-1],
                             y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
             geom_line() + geom_point() +
             scale_x_continuous("quantiles of expression") + 
             scale_y_continuous("Number of DE genes per interval"))
      
      plot(ggplot(data.frame(x=dat$qtl.exp[-1],
                             y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
             geom_line() + geom_point() +
             scale_x_continuous("base mean of expression") + 
             scale_y_continuous("Number of DE genes per interval"))
      
      p <- ggplot(data.frame(x=dat$qtl.exp[-1],
                             y=diff(dat$number.degs[1] - dat$number.degs)),aes(x,y)) + 
        geom_line() + geom_point() +
        scale_x_log10("base mean of expression") + 
        scale_y_continuous("Number of DE genes per interval") + 
        geom_vline(xintercept=expression_cutoff,
                   linetype="dotted",col="red")
      suppressMessages(suppressWarnings(plot(p)))
    }
  }
  
  sel <- res$padj <= padj & abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) & 
    res$baseMean >= expression_cutoff
  
  if(verbose){
    message(sprintf(paste(
      ifelse(sum(sel)==1,
             "There is %s gene that is DE",
             "There are %s genes that are DE"),
      "with the following parameters: FDR <= %s, |log2FC| >= %s, base mean expression > %s"),
      sum(sel),padj,
      lfc,expression_cutoff))
  }
  
  # proceed only if there are DE genes
  if(sum(sel) > 0){
    val <- rowSums(vst[sel,sample_sel,drop=FALSE])==0
    if (sum(val) >0){
      warning(sprintf(paste(
        ifelse(sum(val)==1,
               "There is %s DE gene that has",
               "There are %s DE genes that have"),
        "no vst expression in the selected samples"),sum(val)))
      sel[sel][val] <- FALSE
    } 
    
    if(export){
      if(!dir.exists(default_dir)){
        dir.create(default_dir,showWarnings=FALSE,recursive=TRUE,mode="0771")
      }
      write.csv(res,file=file.path(default_dir,paste0(default_prefix,"results.csv")))
      write.csv(res[sel,],file.path(default_dir,paste0(default_prefix,"genes.csv")))
    }
    if(plot & sum(sel)>1){
      heatmap.2(t(scale(t(vst[sel,sample_sel]))),
                distfun = pearson.dist,
                hclustfun = function(X){hclust(X,method="ward.D2")},
                trace="none",col=hpal,labRow = FALSE,
                labCol=labels[sample_sel],...
      )
    }
  }
  return(list(all=rownames(res[sel,]),
              up=rownames(res[sel & res$log2FoldChange > 0,]),
              dn=rownames(res[sel & res$log2FoldChange < 0,])))
}


#' Data
load(here("data/analysis/salmon/dds.rds"))

#' Normalisation for visualisation
vsd <- varianceStabilizingTransformation(dds,blind=FALSE)
vst <- assay(vsd)
vst <- vst - min(vst)
dir.create(here("data/analysis/DE"),showWarnings=FALSE)
save(vst,file=here("data/analysis/DE/vst-aware.rda"))
write_delim(as.data.frame(vst) %>% rownames_to_column("ID"),
            here("data/analysis/DE/vst-aware.tsv"))

#' line plot for any Gene of interests
goi <- read_lines(here("doc/goi.txt"))
stopifnot(all(goi %in% rownames(vst)))
dev.null <- lapply(goi,line_plot,dds=dds,vst=vst)

#' ## Differential Expression
dds <- DESeq(dds)

#' * Dispersion estimation
#' The dispersion estimation is adequate
plotDispEsts(dds)

#' Check the different contrasts and relevel the dds object to none_0h sample
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

# Define your contrasts as a list of vectors
dds <- readRDS("~/shruti/ShortTermNitrateTimeSeries-SNRIV/data/dds2.rds")
samples <- read.table(here("doc/sampleTremula.tsv"), header = T)
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
  # results_KNO3_vs_0h <- extract_results(dds, vst, contrast_2, plot = FALSE, default_prefix = contrast_name2)
  # Upregulated_KNO3_vs_0h <- results_KNO3_vs_0h$up
  # Downregulated_KNO3_vs_0h <- results_KNO3_vs_0h$dn
  # 
  # # Extract results for KCL_Xh vs None_0h
  # results_KCL_vs_0h <- extract_results(dds, vst, contrast_3, plot = FALSE, default_prefix = contrast_name3)
  # Upregulated_KCL_vs_0h <- results_KCL_vs_0h$up
  # Downregulated_KCL_vs_0h <- results_KCL_vs_0h$dn
  
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

# Loop through each contrast and process
for (contrast in contrasts) {
  res <- process_contrast(contrast, dds, vst, bins_up, bins_dn, results_list)
  results_list <- res$results_list
  bins_up <- res$bins_up
  bins_dn <- res$bins_dn
}

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


