#' ---
#' title: "Short term nitrate time series SNRIV- Biological QA"
#' * Libraries
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(gplots)
  library(here)
  library(magrittr)
  library(parallel)
  library(plotly)
  library(pvclust)
  library(tidyverse)
  library(tximport)
  library(readr)
  library(vsn)
})

#' * Helper functions
source(here("../ERF85GeneExp/UPSCb-common/src/R/featureSelection.R"))

#' * Graphics
hpal <- colorRampPalette(c("blue","white","red"))(100)

# The palette with black:
# hpal <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# To use for fills, add
# scale_fill_manual(values=cbPalette)
# 
# # To use for line and point colors, add
# scale_colour_manual(values=cbPalette)

#' # Data
#' * Sample information
#' ```{r Instructions1,eval=FALSE,echo=FALSE}
#' # The csv file should contain the sample information, including the sequencing file name, 
#' # any relevant identifier, and the metadata of importance to the study design
#' # as columns, e.g. the SamplingTime for a time series experiment
#'  ```
                   
# samples <- read_table(here("doc/kclknoZeroExcl.csv"))
samples <- read_table(here("doc/sampleTremula.csv"))

#' * tx2gene translation table
#' ```{r Instructions2,eval=FALSE,echo=FALSE}
#' # This file is necessary if your species has more than one transcript per gene.
#' #
#' # It should then contain two columns, tab delimited, the first one with the transcript
#' # IDs and the second one the corresponding
#' #
#' # If your species has only one transcript per gene, e.g. Picea abies v1, then
#' # comment the next line
#' ```
# tx2gene <- suppressMessages(read_delim(here("data/salmonT89/tx2gene.tsv"),
#                                        delim="\t", col_names=c("TXID","GENE")))
tx2gene <- suppressMessages(read_delim(here("../single_cell_analysis_poplar/reference/annotation/tx2gene.tsv.gz"),
                                       delim="\t", col_names=c("TXID","GENE")))

#' Read the expression at the gene level
#' ```{r CHANGEME4,eval=FALSE,echo=FALSE}
#' If the species has only one transcript per gene, or if you are conducting QA 
#' in transcript level replace with the following:

# txi <- suppressMessages(tximport(files = samples$Filename, type = "salmon",
#                                  tx2gene=tx2gene))
txi <- suppressMessages(tximport(files = samples$file,type = "salmon",
                                 tx2gene=tx2gene))

counts <- txi$counts
colnames(counts) <- samples$NGI_Id

#' 
#' # Quality Control
#' ## "Not expressed" genes
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))
#T89: For all samples: "8.3% percent (2462) of 29747 genes are not expressed"
#tremula: For all samples: "20.2% percent (7478) of 37075 genes are not expressed"

#' ## Sequencing depth
#' * Let us take a look at the sequencing depth, colouring by treatment
#' ```{r Clevel,eval=FALSE,echo=FALSE}
#' ```
dat <- tibble(x=colnames(counts),y=colSums(counts)) %>% 
  bind_cols(samples)

ggplot(dat,aes(x,y,fill=treatment)) + 
  geom_col() + 
  scale_y_continuous(name="reads") +
  facet_grid(~ factor(treatment), scales = "free") +
  theme_bw() + 
  theme(axis.text.x=element_text(angle=90,size=4),
        axis.title.x=element_blank())

#' ## per-gene mean expression
#' 
#' _i.e._ the mean raw count of every gene across samples is calculated
#' and displayed on a log10 scale.
#' 
ggplot(data.frame(value=log10(rowMeans(counts))),aes(x=value)) + 
  geom_density(na.rm = TRUE) +
  ggtitle("gene mean raw counts distribution") +
  scale_x_continuous(name="mean raw counts (log10)") + 
  theme_bw()

#' ```{r CHANGEME6,eval=FALSE,echo=FALSE}
#' # In the following, the second mutate also needs changing, I kept it 
#' # as an example to illustrate the first line. SampleID would be 
#' # a column in the samples object (the metadata) that uniquely identify
#' # the samples.
#' # If you have only a single metadata, then remove the second mutate call
#' # If you have more, add them as needed.
#' ```
#' 
#' ## Per-sample expression
dat <- as.data.frame(log10(counts)) %>%
  utils::stack() %>%
  mutate(SampleID=samples$NGI_Id[match(ind,samples$NGI_Id)])

ggplot(dat,aes(x=values,group=ind,col=SampleID)) + 
  geom_density(na.rm = TRUE) + 
  ggtitle("sample raw counts distribution") +
  scale_x_continuous(name="per gene raw counts (log10)") + 
  theme_bw()

#' 
#' * Export raw expression data
write.csv(counts,file=here("data/raw-unormalised-gene-exp.csv"))

#' 
#' # Data normalisation 
#' ## Preparation
#' For visualization, the data is submitted to a variance stabilization
#' transformation using _DESeq2_. The dispersion is estimated independently
#' of the sample tissue and replicate. 
#'  
#'  ```{r CHANGEME7,eval=FALSE,echo=FALSE}
#'  # In the following, we provide the expected expression model, based on the study design.
#'  # It is technically irrelevant here, as we are only doing the quality assessment of the data, 
#'  # but it does not harm setting it correctly for the differential expression analyses that may follow.
#'  ```
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
#' because we have onlly one set of 0h, so the model 'time + treat+time*treatment
#' would not work directly unless we either change the sample file information:
#' 0h 2h 4h 8h 12h 24h 48h- time points
#' KCL/None dds, KNO3/None dds, design ~time for both
#' KCL/KNO3 dds, design ~time*treatment or §group
#' KCL/KNO3/None/None dds (change treatment name so None doesn't exist anymore),
#'  design ~time*treatment
#'
#'  Scenario 1: KCL vs KNO3 dds excluding 0h time, then use design ~time*treatment
# For T89:
dds <- DESeqDataSetFromTximport(txi=txi,
                                colData = samples, design =~Time*Treatment)

#' For T89: use LRT with the same settings as KCL vs KNO3 dds excluding 0h time,
#' then use design ~time*treatment
# full_model <- ~ Time + Treatment + Treatment:Time
# For LRT test, provide a reduced model, that is the full model without
# treatment:time term:
# reduced_model <- ~ Treatment + Time

dds <- DESeqDataSetFromTximport(txi =txi, colData = samples,
                                 design = ~ Treatment + Time + Treatment:Time)
# Run in lrtDesq.R
# dds_lrt_time <- DESeq(dds1, test="LRT", reduced = ~ Treatment + Time)
# 
# clusters <- degPatterns(cluster_rlog, metadata = meta, time="Time", 
#                         col="Treatment")

#'  Scenario 2: merge the two variables time and treat and do pairwise compare
#' For tremula
dds <- DESeqDataSetFromTximport(
  txi=txi, colData = samples, design =~treatmentTime)

#'
#' Check the size factors (_i.e._ the sequencing library size effect)
#' 
dds <- estimateSizeFactors(dds)
boxplot(normalizationFactors(dds),
        main="Sequencing libraries size factor",
        las=2,log="y")

#' and without outliers:
boxplot(normalizationFactors(dds),
        main="Sequencing libraries size factor without outlier",
        las=2,log="y", outline=FALSE)
abline(h=1, col = "Red", lty = 3)

#' ## Variance Stabilising Transformation
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)
#'
#' ## Validation
#' 
#' let's look at standard deviations before and after VST normalization. 
#' This plot is to see whether there is a dependency of SD on the mean. 
#' 
#' Before:  
meanSdPlot(log1p(counts(dds))[rowSums(counts(dds))>0,])

#' After VST normalization, the red line is almost horizontal which indicates
#' no dependency of variance on mean (homoscedastic).

meanSdPlot(vst[rowSums(vst)>0,])

#' ## QC on the normalised data
#' ### PCA
pc <- prcomp(t(vst))
percent <- round(summary(pc)$importance[2,]*100)

#' ### Scree plot
#' We define the number of variable of the model: 
nvar = 1

#' An the number of possible combinations
#' ```{r CHANGEME8,eval=FALSE,echo=FALSE}
#' This needs to be adapted to your study design. Add or drop variables aas needed.
#' ```
nlevel=nlevels(dds$Treatment) *nlevels(dds$Time)
nlevel=nlevels(dds$Time)
nlevel=nlevels(dds$treatmentTime)

#' We plot the percentage explained by different components
#' 
#' * the red line represents number of variables in the model  
#' * the orange line represents number of variable combinations.
#' 
ggplot(tibble(x=1:length(percent),y=cumsum(percent)),aes(x=x,y=y)) +
  geom_line() + scale_y_continuous("variance explained (%)",limits=c(0,100)) +
  scale_x_continuous("Principal component") + 
  geom_vline(xintercept=nvar,colour="red",linetype="dashed",size=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nvar],colour="red",linetype="dashed",size=0.5) +
  geom_vline(xintercept=nlevel,colour="orange",linetype="dashed",size=0.5) + 
  geom_hline(yintercept=cumsum(percent)[nlevel],colour="orange",linetype="dashed",size=0.5)

#' ### PCA plot
pc.dat <- bind_cols(PC1=pc$x[,1],
                    PC2=pc$x[,2],
                    as.data.frame(colData(dds)))

p <- ggplot(pc.dat,aes(x=PC1,y=PC2,shape=treatment, col=time)) + 
  geom_point(size=2) + 
  ggtitle("Principal Component Analysis",subtitle="variance stabilized counts")

ggplotly(p) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")))

# add a frame around different sample groups in pca
hull_group <- pc.dat %>%
  dplyr::mutate(SampleID = sample) %>%
  dplyr::group_by(time) %>%
  dplyr::slice(chull(PC1, PC2))

p2 <- ggplot() +
  geom_point(data = pc.dat, mapping = aes(x = PC1, y = PC2, color = treatment, 
                                          shape=time), size = 2) +
  ggplot2::geom_polygon(data = hull_group, aes(x = PC1, y = PC2, fill = treatment, 
                                               group = time), alpha = 0.2)

# ggplotly(p2, tooltip = "label")
ggplotly(p2) %>% 
  layout(xaxis=list(title=paste("PC1 (",percent[1],"%)",sep="")),
         yaxis=list(title=paste("PC2 (",percent[2],"%)",sep="")))

#' ## Sequencing depth
#' Number of genes expressed per condition at different cutoffs:
conds <- factor(paste(dds$Time))
dev.null <- rangeSamplesSummary(counts=vst,conditions=conds,nrep=4)

#' ## Heatmap
#' 
#' Filter for noise
#' 
sels <- rangeFeatureSelect(counts=vst,
                           conditions=conds,
                           nrep=2) %>%
  suppressWarnings()
vst.cutoff <- 2

#' * Heatmap of "all" genes
#' 
hm <- heatmap.2(t(scale(t(vst[sels[[vst.cutoff+1]],]))),
                distfun=pearson.dist,
                hclustfun=function(X){hclust(X,method="ward.D2")},
                labRow = NA,trace = "none",
                labCol = conds,
                col=hpal)

plot(as.hclust(hm$colDendrogram))

#' ## Clustering of samples
#' Done to assess the previous dendrogram's reproducibility
hm.pvclust <- pvclust(data = t(scale(t(vst[sels[[vst.cutoff+1]],]))),
                      method.hclust = "ward.D2", 
                      nboot = 1000, parallel = TRUE)

#' plot the clustering with bp and au
plot(hm.pvclust, labels = conds)
pvrect(hm.pvclust)
#' bootstrapping results as a table
print(hm.pvclust, digits=3)

#' 
#' ```{tech rep, echo=FALSE, eval=FALSE}
#' # First create a new variable in your sample object called BioID that identifies uniquely technical replicates, so one value for all tech rep of the same bio rep
#' samples$BioID <- CHANGEME
#' # Merging technical replicates
#' txi$counts <- sapply(split.data.frame(t(txi$counts),samples$BioID),colSums)
#' txi$length <- sapply(split.data.frame(t(txi$length),samples$BioID),colMaxs)
#' # Counts are now in alphabetic order, check and reorder if necessary
#' stopifnot(colnames(txi$counts) == samples$BioID)
#' samples <- samples[match(colnames(txi$counts),samples$BioID),]
#' # Recreate the dds
#' dds <- DESeqDataSetFromTximport(
#'   txi=txi,
#'   colData = samples,
#'   design = ~ Tissue)
#'```
#'
#' # Session Info
#' <details><summary>Session Info</summary>
#' ```{r session info}
#' sessionInfo()
#' ```
