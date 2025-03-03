# Interaction terms can be added to design formula, in order to test, example, 
# if log2 fold change attributable to a given condition is different based on
# another factor, for eg if the condition effect differs across genotype.

# Case1: https://support.bioconductor.org/p/62684/
# Case1a: find genes showing difference in any time point between Treat and control. 
# Answer: a reduced formula of ~ time + treat means that genes which show a consistent
# difference from time 0 onward will not have a small p value. 
# This is what we think you want, because differences between groups at time 0
# should be controlled for, e.g. random differences between the individuals 
# chosen for each treatment group which are observable before the treatment 
# takes effect.

# In the time series example in the workflow, we include the effect of treatment
# as time=0 to control for these differences, so we use the second reduced design

design(dds) <- ~ time + treat + time:treat
dds <- DESeq(dds, test="LRT", reduced = ~ time + treat)

# Case1b: # include genes which are DE at time 0

# Answer: ...to generate p-values which test for any difference due to treatment 
# including time=0 (and which allows for differences in treatment at different times).
# In other words, if we remove treatment as well as its interaction with time 
# in the reduced model, tests if there is any effect of treatment at any time 
# point include time=0.
full = formula(~ Time + Treat + Time:Treat)
reduced = formula(~ Time)
dds <- DESeq(dds, test="LRT", reduced = ~ time)

# to test the treatment effect only at time=0, use a Wald test:
res <- results(dds, name="treat_yes_vs_no", test="Wald")

# Case2a: https://support.bioconductor.org/p/67600/#67612 
# 1. For KO condition, comparison between each cell type
# 2. For each cell type, compare between KO and WT
# Answer: The easiest way for you to make all these comparisons is to join 
# the two factors and then use contrast to compare in the end:
dds$group <- factor(paste0(dds$celltype, dds$condition))
design(dds) <- ~ group
dds <- DESeq(dds)
resultsNames(dds)
# e.g. for condition KO cell type 2 vs cell type 1
results(dds, contrast=c("group","2KO","1KO")) 
# e.g. for cell type 1 KO vs WT
results(dds, contrast=c("group","1KO","1WT")) 

# Case2b: https://support.bioconductor.org/p/9136775/
# the model matrix is not full rank, so the model cannot be fit as specified. 
# Levels or combinations of levels without any samples have resulted in column(s) 
# of zeros in the model matrix.
# only one timepoint was included for the mock condition
# to compare are the timepoints of the infected conditions contrasted to 
# single mock condition.
# Answer: combine the two variables into one and make pairwise comparisons 

# Case3: Example of LRT for time series: 
# In a likelihood ratio test, the p values and test statistic the stat column) 
# are values for the test that removes all of the variables which are present in
# the full design and not in the reduced design. This tests the null hypothesis 
# that all the coefficients from these variables and levels of these factors are
# equal to zero. LRT p values therefore represent a test of all variables and 
# all levels of factors which are among these variables. However, the results 
# table only has space for one column of log fold change, so a single variable 
# and a single comparison is shown (among the potentially multiple log fold changes
# which were tested in LRT). This is indicated at the top of the results table 
# with the text, e.g., log2 fold change (MLE): condition C vs A, followed by, 
# LRT p-value: batch + condition’ vs ‘~ batch’. This indicates that the p value 
# is for LRT of all the variables and all the levels, while the
# log fold change is a single comparison from among those variables and levels. 

# Case4: Example of time series: https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#time-series-experiments
# 
# Case5: https://support.bioconductor.org/p/65676/#66860
# 5a: 60 vs 0 minutes in mutant
# Answer: strainmut.minute60 is an interaction- names of both variables
# strain and minute. It is a test for if mut vs WT fold change is
# different at minute 60 than at 0. 
# lfc of 60 vs 0 minutes for WT strain: results(dds, name="minute_60_vs_0")
# lfc of 60 vs 0 minutes for mut strain would be sum of WT term above and 
# interaction term which is an additional effect beyond the effect for reference
# (WT): results(dds, contrast=list(c("minute_60_vs_0","strainmut.minute60"))

# 5b: "Imagine a gene that behaves in the same way in all timepoints 
# (between the two groups) but time point 0 to 15. Will this be detected by deseq?"
# Answer: This is the description of an interaction term in the model. 
# find this by testing the interaction term strainmut.minute15:
# results(dds, name="strainmut.minute15")

# 5c: take all time points into account- find a gene DE in all time points 
# between two groups?
# Answer: find significant mut vs wt fold change for each and every time point, 
# not simple. test the fold change for mut vs WT at minute 0:
# results(dds, name="strain_mut_vs_wt")
# ...and then test the fold change for mut vs wt at every other time point:
# results(dds, contrast=list(c("strain_mut_vs_wt","strainmut.minute15")))
# And you would then look for genes which have a small adjusted p-value in 
# all of the results tables.
# first comparison takes into account the difference at t0- recommended.
# The second doesn't take into account the difference at t0- not recommended.
# strainmut.minute15" is the difference between Mut vs WT at minute 15, 
# controlling for baseline. If you add "strain_mut_vs_wt" to this, you get LFC
# for Mutant vs WT at minute 15, not controlling for baseline. So the second one
# is the observed difference at minute 15 between the two groups (because you 
# added in the change that was present at time=0).

# get genes, in just one specific condition, change over time? 
# For reference condition, pull out results for e.g. minute_15_vs_0 with name 
# For non-reference condition, add the main effect minute_15_vs_0 and interaction
# results(dds, contrast=list(c("minute_15_vs_0","strainmut.minute15")))

# Case6: https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/08b_time_course_analyses.html
# We use the LRT to explore whether there are any significant differences across 
# a series of timepoints and further evaluate differences observed between sample classes.
# an experiment looking at the effect of treatment over time on mice of 
# two different genotypes. We could use a design formula for our ‘full model’
# that would include the major sources of variation in our data: genotype, 
# treatment, time, and our main condition of interest, 
# which is the difference in the effect of treatment over time (treatment:time).
full_model <- ~ genotype + treatment + time + treatment:time
# For LRT test, also provide a reduced model, that is the full model without
# treatment:time term:
reduced_model <- ~ genotype + treatment + time

dds <- DESeqDataSetFromMatrix(countData = raw_counts, colData = metadata, 
                              design = ~ genotype + treatment + time + treatment:time)
dds_lrt_time <- DESeq(dds, test="LRT", reduced = ~ genotype + treatment + time)

# Case7: Time course exp by M. Love of DeSeq
# https://bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#time-course-experiments