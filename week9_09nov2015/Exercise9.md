---
title: "Exercise 9"
output: html_document
---


In this exercise, we will explore two popular pipelines for differential expression of RNA-seq data, (already) given the counts.  See Exercise 8 for some details on the counting, but more will follow in coming lectures.  The dataset used here is the so-called 'pasilla' data, which compares the knockout of pasilla (a splicing factor) to a wild-type control.

In this exercise, some code will be given to get started.  But, there are good resources on how to go through these (now) fairly standard analysis pipelines and these are referred to below for more details.

First, we have the samples already organized into a table of metadata and this is used to set the count filenames.


```r
library("edgeR")
samples <- read.table("samples.txt", header=TRUE,
                      row.names=5, stringsAsFactors=FALSE)
samples
```

```
##                   rep condition libtype shortname           countfile
## CG8144_RNAi-1.bam  T1         T      SE     T1.SE CG8144_RNAi-1.count
## CG8144_RNAi-3.bam  T3         T      PE     T3.PE CG8144_RNAi-3.count
## CG8144_RNAi-4.bam  T4         T      PE     T4.PE CG8144_RNAi-4.count
## Untreated-1.bam    C1         C      SE     C1.SE   Untreated-1.count
## Untreated-6.bam    C6         C      SE     C6.SE   Untreated-6.count
## Untreated-3.bam    C3         C      PE     C3.PE   Untreated-3.count
## Untreated-4.bam    C4         C      PE     C4.PE   Untreated-4.count
```


Here, we read in the count files and consolidate it into a single table as opposed to 7 individual files -- the readDGE() function saves having to do this manually:


```r
counts <- readDGE(samples$countfile)$counts
```

Here, we can trim the column names to get a 'nice' grouping variable to use in the design matrix:


```r
grp <- gsub("\\-.[0-9]*","",colnames(counts))
grp
```

```
## [1] "CG8144_RNAi" "CG8144_RNAi" "CG8144_RNAi" "Untreated"   "Untreated"  
## [6] "Untreated"   "Untreated"
```


Below is a pipeline using likelihood-ratio tests that is adapted from Quick Start section, early in the [edgeR user's guide](http://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)



```r
group <- factor(c(1,1,2,2))
y <- DGEList(counts=x,group=group)
y <- calcNormFactors(y)
design <- model.matrix(~group)

y <- estimateDisp(y,design)
#To perform likelihood ratio tests:
fit <- glmFit(y,design)
lrt <- glmLRT(fit,coef=2)
topTags(lrt)
```

This will lead you to a list of putative/detected differentially expressed genes.

#### Exercise 1: Fix the code above to work on the count table that was already read in.  Add in some spot checks, including an MDS plot from plotMDS(), a look at the dispersion-mean plot using plotBCV() and a look at the overall M vs A plot using plotSmear().

As an alternative, one can run through a standard voom pipeline, such as:


```r
v <- voom(d$counts, design=model.matrix(~grp), plot = TRUE)
vf = lmFit(v, design = model.matrix(~grp))  # note similarity to limma
                                            # analysis performed earlier
vf = eBayes(vf)
```

For more details, see Section 15.3 Differential expression of the [limma user's guide](http://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf).

#### Exercise 2. Again, the above code for voom is not run here.  Fix this code to work with the same count dataset and then make some comparisons of the genes called DE.  For example, make a plot of the plot of estimated fold changes for the two methods, P-values, or a venn diagram of the called sets of DE genes at a set threshold.
