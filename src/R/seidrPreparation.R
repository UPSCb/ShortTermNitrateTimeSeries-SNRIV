#' ---
#' title: "Network data preparation"
#' author: "Happy friends and Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(tidyverse)
})

#' Helper
source(here("UPSCb-common/src/R/featureSelection.R"))

#' # Data
#' ```{r CHANGEME2,eval=FALSE,echo=FALSE}
#'  CHANGEME is the variance stabilised data, where the transformation was done taking the
#'  model into account (_i.e._ `blind=FALSE`)
#' ```
load(here("data/analysis/DE/vst-noNA.rda"))

vst <- as_tibble(vst_subset, rownames="ID") %>% 
  column_to_rownames("ID")

#' # Filter
#' Different ways of filtering are considered. An cutoff on the absolute expression level
#' or on the variance (standard score), respectively. The latter might be more adequate to
#' avoid filtering lowly expressed genes that show sufficient, metadata relevant, variance.
samples <- read_tsv(here("doc/samples_full_rank.txt"))

#Check order
samples <- samples %>% 
  mutate(VstNames = paste(sample, treatment, time, replicate,sep="_"))


sels <- rangeFeatureSelect(counts=as.matrix(vst),
                           conditions=factor(samples$treatmentTime),
                           nrep=3)

e.cutoff <- 1

s.sels <- rangeFeatureSelect(counts=as.matrix(vst),
                             conditions=factor(samples$treatmentTime),
                             nrep=3,scale=TRUE)
s.cutoff <- 1

table(s.sels[[s.cutoff]],sels[[e.cutoff  + 1]])
sum(featureSelect(as.matrix(vst),
                  factor(samples$treatmentTime),
                  exp=0.1,3))
sum(featureSelect(as.matrix(vst),
                  factor(samples$treatmentTime),
                  exp=0.05,3,scale=TRUE))

#' ```{r CHANGEME3,eval=FALSE,echo=FALSE}
#'  CHANGEME is the vst cutoff devised from the plot above. The goal is to remove / reduce
#'  the signal to noise. Typically, this means trimming the data after the first sharp decrease 
#'  on the y axis, most visible in the non logarithmic version of the plot
#' ```
vst.cutoff <- 0.5

#' # Export
dir.create(here("data/seidr"), recursive=TRUE)

#' * gene by column, without names matrix
write.table(t(vst[sels[[vst.cutoff+1]],]),
            file=here("data/seidr/headless.tsv"),
            col.names=FALSE,
            row.names=FALSE,
            sep="\t",quote=FALSE)

#' * gene names, one row
write.table(t(sub("\\.[0-9]+$","",rownames(vst)[sels[[vst.cutoff+1]]])),
            file=here("data/seidr/genes.tsv"),
            col.names=FALSE,
            row.names=FALSE,
            sep="\t",quote=FALSE)

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
