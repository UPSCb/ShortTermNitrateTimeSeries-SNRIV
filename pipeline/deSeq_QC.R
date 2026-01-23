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

#' samples
samples <- read_table(here("doc/sampleTremula.tsv"))

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
tx2gene <- suppressMessages(read_delim(here("../single_cell_analysis_poplar/reference/annotatio/tx2gene.tsv.gz"),
                                       delim="\t", col_names=c("TXID","GENE")))

#' Read the expression at the gene level
#' ```{r CHANGEME4,eval=FALSE,echo=FALSE}
#' If the species has only one transcript per gene, or if you are conducting QA 
#' in transcript level replace with the following:
txi <- suppressMessages(tximport(files = samples$file,type = "salmon",
                                 tx2gene=tx2gene))

counts <- txi$counts
colnames(counts) <- samples$sample

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
#' This needs to be adapted to your study design. Add or drop variables as needed.
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
#' # Session Info
#' <details><summary>Session Info</summary>
#' ```{r session info}
#' sessionInfo()
#' ```
