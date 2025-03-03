suppressMessages({
  source(here("UPSCb-common/Rtoolbox/src/plotEnrichedTreemap.R"))
  source(here("UPSCb-common/src/R/featureSelection.R"))
  source(here("UPSCb-common/src/R/volcanoPlot.R"))
  source(here("UPSCb-common/src/R/topGoUtilities.R"))
})

suppressPackageStartupMessages({
  library(here)
  library(matrixStats)
  library(Mfuzz)
  library(readr)
  library(Biobase)
})

# Clustering of dds from Case 4.1: 
# identify groups of genes that share a pattern of expression change across
# sample groups (levels). 

# `degPatterns` uses a hierarchical clustering approach based on pair-wise
# correlations between genes, then cuts the hierarchical tree to generate groups
# of genes with similar expression profiles. The tool cuts tree in a way 
# to optimize the diversity of the clusters, such that the 
# variability inter-cluster > the variability intra-cluster.

# Before clustering, first subset our rlog transformed normalized counts to 
# retain only the differentially expressed genes (padj < 0.01). 
# may take some time depending on the number of genes
# vst recommende for n>30
# Obtain vst (faster)/rlog values for those significant genes
# blind = FALSE, means that differences between time and treatment 
# (the variables in the design) will not contribute to the expected 
# variance-mean trend of the experiment. The experimental design is not used
# directly in the transformation, only in estimating the global amount of 
# variability in the counts. For a fully unsupervised transformation, 
# one can set blind = TRUE (which is the default).

# vst transformed counts for the significant genes are input to `degPatterns` 
# to show gene clusters across sample groups. additional arguments required are:
# metadata`: the metadata dataframe that corresponds to samples
# time`: character column name in metadata that is used as variable that changes
# col`: character column name in metadata to separate samples

# head(clusters$df)
# group1 <- clusters$df %>% dplyr::filter(cluster == 1)

# To reduce the size of the object, and to increase the speed of our functions, 
# remove the rows that have no or nearly no information about the amount of gene
# expression. So do pre-filtering to keep only rows that have a count of 
# at least 10 for a minimal number of samples. The count of 10 is a reasonable 
# choice for bulk RNA-seq. A recommendation for the minimal number of samples 
# is to specify the smallest group size, e.g. here there are atleast 3 samples 
# in each group. If there are not discrete groups, one can use the minimal 
# number of samples where non-zero counts would be considered interesting. 
# nrow(dds2)
# smallestGroupSize <- 3
# keep <- rowSums(counts(dds2) >= 10) >= smallestGroupSize
# dds <- dds2[keep,]
# nrow(dds)

# Extract results for LRT
dds2 <- readRDS("~/shruti/ShortTermNitrateTimeSeries-SNRIV/data/dds2.rds")

# also from case 4.5.a.: reduced = ~ treatment + time
dds2 <- readRDS("~/shruti/ShortTermNitrateTimeSeries-SNRIV/data/dds.rds")

vsd <- varianceStabilizingTransformation(dds2, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)

# clusters <- degPatterns(cluster_vst, metadata = samples, time = "treatment",
#                         col="treatmentTime")

samples <- read_table(here("doc/sampleTremula.csv"))
colnames(vst) <- dds2$sample
# conds <- apply(samples[match(colnames(vst),samples$sample),
#                        c("treatmentTime")],
#                1,paste,collapse="-")
conds <- apply(samples[match(colnames(vst),samples$sample),
                       c("treatment", "time")],
               1,paste,collapse="-")
eset <- ExpressionSet(sapply(split.data.frame(t(vst),conds),colMeans))

#' Remove genes with too little a variation
#' 
#' First we look at the SD distribution
plot(density(rowSds(vst)))
plot(density(rowSds(vst)),xlim=c(-0.2,0.5))

#' A cutoff at 0.1 seems adequate
eset <- filter.std(eset,min.std=0.1)

#' Standardise the values
eset <- standardise(eset) 

#' ## Clustering
#' parameter estimation
m1 <- mestimate(eset)

#' cluster
cl <- mfuzz(eset,m1,c=24)

#' ##Plot
colnames(eset)
dir.create(here("data/Mfuzz"),showWarnings = FALSE)
pdf(file=here("data/Mfuzz/clusters50.pdf"),width = 16,height=24)
mfuzz.plot2(eset,cl,x11 = FALSE,mfrow=c(3,10),time.labels = colnames(eset),
            centre = TRUE,las=2)
dev.off()

#' ##Membership
str(cl)
barplot(cl$size)

# genes for cluster 20
names(cl$cluster)[cl$cluster == 20]

# enrichment of all clusters
threshold <- 0.1
min_replicates <- 2

replicate_means <- t(apply(vst, 1, function(row) {
  tapply(row, c(samples$time,samples$treatment), mean)
}))

replicate_counts <- rowSums(replicate_means >= threshold)
selected_genes <- replicate_counts >= min_replicates
background <- rownames(vst)[selected_genes]

# goannot <- prepAnnot(mapping = "/mnt/reference/Populus-tremula/v2.2/annotation/blast2go/Potra22_blast2go_GO_export.txt.gz")
goannot <- prepAnnot(mapping = "/mnt/picea/storage/reference/Populus-tremula/v2.2/gopher/gene_to_go.tsv")
# goannot <- prepAnnot(mapping = "/mnt/reference/Arabidopsis-thaliana/TAIR10/gopher/tair10_gene_to_go.tsv")

suppressMessages(enr.list <- lapply(1:length(cl$size),function(i,clc){lapply(names(clc)[clc==i],topGO,background=background,annotation = goannot)},cl$cluster))

# TODO
# export enr.list as in DifferentialExpression.R to files
# I chose the one from Fai as follows:

extractEnrichmentResults <- function(enrichment,
                                     go.namespace=c("BP","CC","MF"),
                                     count=100,plot=TRUE){
  
  # sanity
  if(is.null(unlist(enrichment)) | length(unlist(enrichment)) == 0){
    message("No GO enrichment for",names(enrichment))
  } else {
    if(plot){
      gocatname <- c(BP="Biological Process",
                     CC="Cellular Component",
                     MF="Molecular Function")
      lapply(names(enrichment),function(n){
        lapply(names(enrichment[[n]]),function(gocat){
          dat <- enrichment[[n]][[gocat]]
          if(is.null(dat)){
            message("No GO enrichment for ",n," in category ",gocatname[gocat])
          } else {
            dat$GeneRatio <- dat$Significant/dat$Annotated
            dat$adjustedPvalue <- as.numeric(dat$FDR)
            dat$Count <- as.numeric(dat$Significant)
            dat <- dat[order(dat$GeneRatio),]
            if(nrow(dat) > count){ dat <- dat[1:count,] }
            dat$Term <- factor(dat$Term, levels = unique(dat$Term))
            ggplot(dat, aes(x =Term, y = GeneRatio, color = adjustedPvalue, size = Count)) + 
              geom_point() +
              scale_color_gradient(low = "red", high = "blue") +
              theme_bw() + 
              ylab("GeneRatio") + 
              xlab("") + 
              ggtitle(paste0("GO enrichment: ",n," ",gocatname[gocat])) +
              coord_flip()
          }
        })
      })
    }
  }
}

suppressWarnings(extractEnrichmentResults(enr.list, count = 30))

# cluster membership
# find the highest membership score per gene
max.membership <- sapply(1:nrow(cl$membership),function(i,m,p){
  m[i,p[i]]
},cl$membership,cl$cluster)

# plot it
plot(density(max.membership))

# find the genes in every cluster that is above a given membership value (0.5 in the example below)
cluster.genes <- lapply(1:length(cl$size),function(i){
  dat <- cl$membership[cl$cluster == i,i]
  plot(density(dat),main=paste("cluster",i))
  abline(v=0.75,lty=2,lwd=2,col="gray")
  names(cl$cluster)[cl$cluster == i][cl$membership[cl$cluster == i,i] >= 0.5]
})

# assess how many genes are still shared between clusters
sum(duplicated(unlist(cl.specific.enr.list)))

# run the enrichment
cl.specific.enr.list <- lapply(cluster.genes,gopher,background=background,task="go",url="athaliana")

# run all enrichment at once
# cl.specific.enr.list <- lapply(cluster.genes,gopher,background=background,task=c("go","kegg","mapman"),url="athaliana")
suppressMessages(
  cl.specific.enr.list <- lapply(cluster.genes,topGO,background=background,
                               annotation = goannot, alpha=0.05, p.adjust="BH")
)

saveRDS(cluster.genes,file=here("data/clusterMember.rds"))
saveRDS(eset,file=here("data/eset.rds"))
saveRDS(cl,file=here("data/cluster.rds"))

dev.null <- lapply(1:length(cl.specific.enr.list),function(i){
  r <- cl.specific.enr.list[[i]]
  write_tsv(r$go,path=file.path(file.path(here("data/analysis/Mfuzz",
                                               paste0("cluster_",i,"_GO-enrichment.tsv")))))
  write_tsv(r$go[,c("id","padj")],path=file.path(file.path(here("data/analysis/Mfuzz",
                                                                paste0("cluster_",i,"_GO-enrichment_for-REVIGO.tsv")))))
  # write_tsv(r$kegg,path=file.path(file.path(here("data/analysis/Mfuzz",
  #                                              paste0("cluster_",i,"_KEGG-enrichment.tsv")))))
  # write_tsv(r$kegg[,c("id","padj")],path=file.path(file.path(here("data/analysis/Mfuzz",
  # ....                                                              
  
})

# Overlap analysis - check how much overlap there is between clusters and
# how similar clusters are # useful to refine the number of clusters
O <- overlap(cl)
Ptmp <- overlap.plot(cl,over=O,thres=0.05)
overlap.plot(cl,over=O,thres=0.1)