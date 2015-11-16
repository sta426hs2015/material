---
title: "Exercise 10"
output: html_document
---

In the first part of this exercise, we will revisit the analysis from last week, to account for a covariate.  In the second part, we will use (preprocessed) exon-level counts to look for changes in splicing.


```r
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

#### Exercise 1. Take the data from last week and produce an MDS plot again, but this time colour the points according to the covariate in the samples table: libtype (SE = single end, PE = paired end); perhaps also label the points on the MDS plot using the 'shortname' column to make them easy to distinguish.  Comment on the organization of the samples.

#### Exercise 2. Put a factor variable for the 'libtype' covariate in the design matrix and redo the edgeR or limma analysis from last week (i.e., include also the biological factor in the design matrix).  Compare the set of genes called DE from last week's exercise (i.e., without accounting for the covariate) to this analysis.  Identify and plot (normalized) expression levels of a gene that is affected solely by library type.

Next, we will explore "differential splicing", using the same pasilla dataset (Note: this was used in the pasilla manuscript).


```r
library(pasilla)
sdir <- file.path(system.file(package="pasilla"),"extdata/")
anno <- file.path(sdir, "Dmel.BDGP5.25.62.DEXSeq.chr.gff")

samplesX = data.frame(cond = rep( c("trt","untrt"), c(3,4) ),
                      type = c("SE","PE","PE","SE","SE","PE","PE"),
                      countfile = dir(sdir,pattern="fb.txt"),
                      stringsAsFactors = TRUE)
samplesX
```

```
##    cond type        countfile
## 1   trt   SE   treated1fb.txt
## 2   trt   PE   treated2fb.txt
## 3   trt   PE   treated3fb.txt
## 4 untrt   SE untreated1fb.txt
## 5 untrt   SE untreated2fb.txt
## 6 untrt   PE untreated3fb.txt
## 7 untrt   PE untreated4fb.txt
```

Below is some unevaluated code that represents a standard DEXSeq pipeline:



```r
dxd <- DEXSeqDataSetFromHTSeq(
           countfiles=file.path( sdir, filename ),
           sampleData = samples,
           design = ~ sample + exon + type:exon + condition:exon)
dxd <- estimateSizeFactors( dxd )
dxd <- estimateDispersions( dxd )
dxd <- testForDEU( dxd )
dxr <- DEXSeqResults( dxd )
```

As usual, refer to the vignette for the DEXSeq or the documentation for further details.

#### Exercise 3 (optional). Fix the above code to run a standard DEXSeq analysis and plot one of the top differentially spliced genes -- for example, see the plotDEXSeq() function.
