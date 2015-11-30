---
title: "Exercise 12"
output: html_document
---

This exercise is designed to introduce you to two types of DNA methylation data: from the Illumina 450k array (bisulphite sequencing of DNA + microarray readout), from sequencing of bisulphite-treated DNA.  Most of what is below is snippets from a previously-established pipeline, not all of the details of which are given.  But, you can look these up in the documentation for the packages we use.

The data consists of DNA methylation measurements (2 platforms) for replicates of [LNCaP](http://www.lgcstandards-atcc.org/products/all/CRL-1740.aspx?geo_country=ch) cells (metastatic prostate cancer cells from lymph node) and PrEC (prostate epithelial cells).  We will compare these two platforms and run an algorithm to find differentially methylated regions.

First, we load some new packages (install these on your system if needed):


```r
library("minfi")
```

```
## Loading required package: BiocGenerics
## Loading required package: parallel
## 
## Attaching package: 'BiocGenerics'
## 
## Die folgenden Objekte sind maskiert von 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## 
## Das folgende Objekt ist maskiert 'package:stats':
## 
##     xtabs
## 
## Die folgenden Objekte sind maskiert von 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, do.call, duplicated, eval, evalq, Filter, Find, get,
##     intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rep.int, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unlist, unsplit
## 
## Loading required package: Biobase
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
## 
## Loading required package: lattice
## Loading required package: GenomicRanges
## Loading required package: S4Vectors
## Loading required package: stats4
## Creating a generic function for 'nchar' from package 'base' in package 'S4Vectors'
## Loading required package: IRanges
## Loading required package: GenomeInfoDb
## Loading required package: Biostrings
## Loading required package: XVector
## Loading required package: bumphunter
## Loading required package: foreach
## foreach: simple, scalable parallel programming from Revolution Analytics
## Use Revolution R for scalability, fault tolerance and more.
## http://www.revolutionanalytics.com
## Loading required package: iterators
## Loading required package: locfit
## locfit 1.5-9.1 	 2013-03-22
## Setting options('download.file.method.GEOquery'='auto')
```

```r
library("charm")
```

```
## Loading required package: SQN
## Loading required package: mclust
##     __  ___________    __  _____________
##    /  |/  / ____/ /   / / / / ___/_  __/
##   / /|_/ / /   / /   / / / /\__ \ / /   
##  / /  / / /___/ /___/ /_/ /___/ // /    
## /_/  /_/\____/_____/\____//____//_/    version 5.1
## Type 'citation("mclust")' for citing this R package in publications.
## Loading required package: nor1mix
## Loading required package: fields
## Loading required package: spam
## Loading required package: grid
## Spam version 1.3-0 (2015-10-24) is loaded.
## Type 'help( Spam)' or 'demo( spam)' for a short introduction 
## and overview of this package.
## Help for individual functions is also obtained by adding the
## suffix '.spam' to the function name, e.g. 'help( chol.spam)'.
## 
## Attaching package: 'spam'
## 
## Das folgende Objekt ist maskiert 'package:stats4':
## 
##     mle
## 
## Die folgenden Objekte sind maskiert von 'package:base':
## 
##     backsolve, forwardsolve
## 
## Loading required package: maps
## 
##  # ATTENTION: maps v3.0 has an updated 'world' map.        #
##  # Many country borders and names have changed since 1990. #
##  # Type '?world' or 'news(package="maps")'. See README_v3. #
## 
## 
## 
## Attaching package: 'maps'
## 
## Das folgende Objekt ist maskiert 'package:mclust':
## 
##     map
## 
## Loading required package: RColorBrewer
## Loading required package: genefilter
## 
## Attaching package: 'genefilter'
## 
## Das folgende Objekt ist maskiert 'package:base':
## 
##     anyNA
## 
## Welcome to charm version 2.14.0
## 
## Attaching package: 'charm'
## 
## Das folgende Objekt ist maskiert 'package:minfi':
## 
##     qcReport
## 
## Das folgende Objekt ist maskiert 'package:bumphunter':
## 
##     clusterMaker
```

```r
library("GenomicRanges")
```

After you have unzipped the exercise files, you can follow the steps below to: i) read in a sheet of metadata; ii) read in the raw probe-level data from IDAT files (the format used by Illumina); iii) preprocess the data.


```r
f <- read.450k.sheet(".", "450k_SampleSheet", recursive = FALSE)
```

```
## [read.450k.sheet] Found the following CSV files:
## [1] "./450k_SampleSheet.csv"
```

```
## Warning in readLines(file): unvollständige letzte Zeile in './
## 450k_SampleSheet.csv' gefunden
```

```r
# hack to get software to recognize
f$Basename <- paste0(f$Slide,"_",f$Array) 

raw <- read.450k.exp(base="6264509024/", 
                     targets = f, verbose=TRUE, recursive = TRUE)
```

```
## [read.450k] Reading 6264509024_R01C02_Grn.idat 
## [read.450k] Reading 6264509024_R02C02_Grn.idat 
## [read.450k] Reading 6264509024_R03C02_Grn.idat 
## [read.450k] Reading 6264509024_R04C02_Grn.idat 
## [read.450k] Reading 6264509024_R05C02_Grn.idat 
## [read.450k] Reading 6264509024_R06C02_Grn.idat 
## [read.450k] Reading 6264509024_R01C02_Red.idat 
## [read.450k] Reading 6264509024_R02C02_Red.idat 
## [read.450k] Reading 6264509024_R03C02_Red.idat 
## [read.450k] Reading 6264509024_R04C02_Red.idat 
## [read.450k] Reading 6264509024_R05C02_Red.idat 
## [read.450k] Reading 6264509024_R06C02_Red.idat 
## [read.450k] Read idat files in  9.859 seconds
## [read.450k] Creating data matrices ... done in 0.954 seconds
## [read.450k] Instantiating final object ... done in 1.926 seconds
```

```r
ppi <- preprocessIllumina(raw, bg.correct = TRUE, normalize = "controls")
```

```
## Loading required package: IlluminaHumanMethylation450kmanifest
```

Next, we pull out the *beta* values (consensus methylation levels) in to a matrix:


```r
b <- getBeta(ppi)
colnames(b) <- f$Sample_Name  # give descriptive column names
head(b)
```

```
##               LNCaP_1    LNCaP_2    LNCaP_3     PrEC_1     PrEC_2
## cg00050873 0.48922038 0.52516234 0.44031977 0.83689336 0.82305069
## cg00212031 0.97269231 0.99762244 0.99946314 0.00000000 0.01416765
## cg00213748 0.74040632 0.68767123 0.61498708 0.92773109 0.96774194
## cg00214611 0.48949580 0.60976409 0.45129163 0.02387068 0.01940594
## cg00455876 0.77377220 0.72318257 0.77623643 0.80812400 0.95655254
## cg01707559 0.07329267 0.07586994 0.05261299 0.02813853 0.01836801
##                PrEC_3
## cg00050873 0.84062681
## cg00212031 0.01217985
## cg00213748 0.89467593
## cg00214611 0.01533785
## cg00455876 0.86848677
## cg01707559 0.03179948
```

To make the problem a little smaller for the exercise, we subset here to just the probes for chromosome 22:


```r
anno <- read.csv("450_probes_subset.csv", stringsAsFactors=FALSE, header=TRUE)
probeids <- intersect(rownames(b), anno$Name)
```

With a couple steps of matching identifiers, we can align the table of data and the table of annotation information:


```r
m <- match(probeids, anno$Name)
anno <- anno[m, ]

m <- match(probeids, rownames(b))
b <- b[m, ]

# spot checks
all( rownames(b) == anno$Name )
```

```
## [1] TRUE
```

```r
nrow(anno) == nrow(b)
```

```
## [1] TRUE
```


Next, we filter out some probes because some of the methods (to find differentially methylated regions) below do not handle them well:


```r
k <- !is.na(anno$MAPINFO) & rowSums( is.na(b) )==0
anno <- anno[k,]
b <- b[k,]
```

#### Question 1.   Make histograms of "beta" values (3 LNCaP replicates, 3 PrEC replicates).  Are there (global) differences between cancer and  normal?  What if you restrict to probes only in CpG islands?  (You can get this information for chr22 in the annotation file).


Next, we use a *GRanges* object as a container for the location information (you will see why below):


```r
# use a 'GRanges' object to save the location info
gr450k <- GRanges(seqnames=paste("chr",anno$CHR,sep=""), 
                  ranges=IRanges(start=anno$MAPINFO,width=1))
gr450k
```

```
## GRanges object with 6995 ranges and 0 metadata columns:
##          seqnames               ranges strand
##             <Rle>            <IRanges>  <Rle>
##      [1]    chr22 [18632618, 18632618]      *
##      [2]    chr22 [43253521, 43253521]      *
##      [3]    chr22 [24302043, 24302043]      *
##      [4]    chr22 [30901640, 30901640]      *
##      [5]    chr22 [39883192, 39883192]      *
##      ...      ...                  ...    ...
##   [6991]    chr22 [27837205, 27837205]      *
##   [6992]    chr22 [28073464, 28073464]      *
##   [6993]    chr22 [31690508, 31690508]      *
##   [6994]    chr22 [36895405, 36895405]      *
##   [6995]    chr22 [48731367, 48731367]      *
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```


The RRBS data ("reduced representaton" bisulphite sequencing; essentially a protocol that selects for regions near CpG islands) to compare against was collected by a consortium and they make their (processed) data available in BED format.  In the code below, I've selected a subset of columns to get only the information we want:


```r
rf <- dir(,".bedRrbs")
names(rf) <- gsub("_chr22.bedRrbs","",rf)

rd <- lapply(rf, function(u) {
  d <- read.table(u, header = FALSE, stringsAsFactors = FALSE)
  d[,c(1:3,5:6,11)]
})

# have a look at how RRBS data comes
lapply(rd,head)
```

```
## $LNCAP1
##      V1       V2       V3  V5 V6 V11
## 1 chr22 16106635 16106636 239  +  15
## 2 chr22 16106731 16106732  74  -  28
## 3 chr22 16106742 16106743  74  -  39
## 4 chr22 16107857 16107858   2  +   0
## 5 chr22 16157561 16157562  51  +  73
## 6 chr22 16157575 16157576  51  +  84
## 
## $LNCAP2
##      V1       V2       V3  V5 V6 V11
## 1 chr22 16106635 16106636 212  +  14
## 2 chr22 16106731 16106732  99  -  31
## 3 chr22 16106742 16106743  99  -  27
## 4 chr22 16157561 16157562  46  +  65
## 5 chr22 16157575 16157576  46  +  80
## 6 chr22 16157576 16157577   5  -  80
## 
## $PREC1
##      V1       V2       V3  V5 V6 V11
## 1 chr22 16106635 16106636 132  +  45
## 2 chr22 16106731 16106732  25  -  96
## 3 chr22 16106742 16106743  25  -  76
## 4 chr22 16157561 16157562  33  +  73
## 5 chr22 16157575 16157576  33  +  64
## 6 chr22 16157576 16157577   3  -  67
## 
## $PREC2
##      V1       V2       V3 V5 V6 V11
## 1 chr22 16106635 16106636 13  +  69
## 2 chr22 16106731 16106732  1  - 100
## 3 chr22 16106742 16106743  1  - 100
## 4 chr22 16157561 16157562 38  +  68
## 5 chr22 16157575 16157576 38  +  45
## 6 chr22 16157576 16157577  4  - 100
```

In this table, *V11* is the percent methylation and *V5* is the total number of reads at that position.

Below is a clunky function that was written several years ago to collapse the information across both strands to a single position (in principle, the positive and negative positions offset by a base are measuring the same thing):


```r
mapNeg2PosRRBS <- function(p) {
  xs <- split(p, p$V6)
  pos <- xs[["+"]]
  neg <- xs[["-"]]

  key <- paste(pos$V1, pos$V3-1, sep = ".")
  keyn <- paste(neg$V1, neg$V2, sep = ".")

  # finding the pairs
  m.pos <- match(key, keyn)
  n <- is.na(m.pos)
 
  # positive strand
  x <- data.frame(chr=pos$V1, position=pos$V2+1, n=pos$V5, 
                  nC=round( pos$V5*pos$V11/100 ))

  # negative strand that match with positive strand
  mn <- m.pos[!n]
  x$n[!n] <- x$n[!n] + neg$V5[mn]
  x$nC[!n] <- x$nC[!n] + round( neg$V5[mn]*neg$V11[mn]/100 )
 
  # negative strand only
  m.neg <- match(keyn, key)
  nas.neg <- is.na(m.neg)
  y <- data.frame(chr=neg$V1, position=neg$V2, n=neg$V5,
                  nC=round( neg$V5*neg$V11/100 ))[nas.neg, ]

  res <- rbind(x,y)
  GRanges(seqnames=res$chr, IRanges(start=res$position, width=1),
          n=res$n, nC=res$nC)
}
```

This function can be applied to every element of the list (the tables that we read in above) and returns a list of GRanges objects:


```r
rdc <- lapply(rd, mapNeg2PosRRBS)
lapply(rdc,head) # spot check
```

```
## $LNCAP1
## GRanges object with 6 ranges and 2 metadata columns:
##       seqnames               ranges strand |         n        nC
##          <Rle>            <IRanges>  <Rle> | <integer> <numeric>
##   [1]    chr22 [16106636, 16106636]      * |       239        36
##   [2]    chr22 [16107858, 16107858]      * |         2         0
##   [3]    chr22 [16157562, 16157562]      * |        51        37
##   [4]    chr22 [16157576, 16157576]      * |        51        43
##   [5]    chr22 [16157584, 16157584]      * |        51        49
##   [6]    chr22 [16157595, 16157595]      * |        51        48
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
## 
## $LNCAP2
## GRanges object with 6 ranges and 2 metadata columns:
##       seqnames               ranges strand |         n        nC
##          <Rle>            <IRanges>  <Rle> | <integer> <numeric>
##   [1]    chr22 [16106636, 16106636]      * |       212        30
##   [2]    chr22 [16157562, 16157562]      * |        46        30
##   [3]    chr22 [16157576, 16157576]      * |        46        37
##   [4]    chr22 [16157584, 16157584]      * |        46        43
##   [5]    chr22 [16157595, 16157595]      * |        46        41
##   [6]    chr22 [17045625, 17045625]      * |         1         0
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
## 
## $PREC1
## GRanges object with 6 ranges and 2 metadata columns:
##       seqnames               ranges strand |         n        nC
##          <Rle>            <IRanges>  <Rle> | <integer> <numeric>
##   [1]    chr22 [16106636, 16106636]      * |       132        59
##   [2]    chr22 [16157562, 16157562]      * |        33        24
##   [3]    chr22 [16157576, 16157576]      * |        33        21
##   [4]    chr22 [16157584, 16157584]      * |        33        19
##   [5]    chr22 [16157595, 16157595]      * |        33        21
##   [6]    chr22 [17047244, 17047244]      * |        28        20
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
## 
## $PREC2
## GRanges object with 6 ranges and 2 metadata columns:
##       seqnames               ranges strand |         n        nC
##          <Rle>            <IRanges>  <Rle> | <integer> <numeric>
##   [1]    chr22 [16106636, 16106636]      * |        13         9
##   [2]    chr22 [16157562, 16157562]      * |        38        26
##   [3]    chr22 [16157576, 16157576]      * |        38        17
##   [4]    chr22 [16157584, 16157584]      * |        38        16
##   [5]    chr22 [16157595, 16157595]      * |        38        25
##   [6]    chr22 [17047244, 17047244]      * |        55        47
##   -------
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```


#### Exercise 2. Plot the concordance of 450k and RRBS data (for the subset where they both measure the same site).  Does the concordance between 450k and RRBS improve when you require a higher depth in the RRBS data? n.b. a useful function to match up coordinates is findOverlaps().



Next, we run one of the pipelines discussed in lectures: finding differentially methylated regions (DMRs).


```r
# make design matrices for DMR finding 
# ('d1' with group variable, 'd0' just intercept)
grp <- gsub("_[1-3]","",colnames(b))
d1 <- model.matrix(~grp)
d0 <- d1[,1,drop=FALSE]

b1 <- b*.98+.01
range(b)
```

```
## [1] 0 1
```

```r
range(b1)
```

```
## [1] 0.01 0.99
```

```r
o <- order(anno$CHR, anno$MAPINFO)

# find differentially methylated regions
library(charm)
pns <- clusterMaker(anno$CHR, anno$MAPINFO, maxGap=500)
dmrs <- dmrFind(p=b1[o,], mod=d1, mod0=d0, coef=2, pns=pns[o],
                chr=anno$CHR[o], pos=anno$MAPINFO[o], svs=0, 
                use.limma=TRUE, use="swald", Q=0.97)
```

```
## 
## Regression (limma)
## Obtaining estimates for  grpPrEC 
## Smoothing
## ===
```

```
## Warning in lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth, : Estimated
## rdf < 1.0; not estimating variance
```

```
## Warning in lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth, : Estimated
## rdf < 1.0; not estimating variance
```

```
## ==
```

```
## Warning in lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth, : Estimated
## rdf < 1.0; not estimating variance
```

```
## ==
```

```
## Warning in lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth, : Estimated
## rdf < 1.0; not estimating variance
```

```
## ==============
```

```
## Warning in lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth, : Estimated
## rdf < 1.0; not estimating variance
```

```
## ============================
```

```
## Warning in lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth, : Estimated
## rdf < 1.0; not estimating variance
```

```
## ======================
```

```
## Warning in lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth, : Estimated
## rdf < 1.0; not estimating variance
```

```
## ===
```

```
## Warning in lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth, : Estimated
## rdf < 1.0; not estimating variance
```

```
## Warning in lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth, : Estimated
## rdf < 1.0; not estimating variance
```

```
## Warning in lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth, : Estimated
## rdf < 1.0; not estimating variance
```

```
## Warning in lfproc(x, y, weights = weights, cens = cens, base = base, geth = geth, : Estimated
## rdf < 1.0; not estimating variance
```

```
## =
## ==============================================================================================
## ..
```

```
## Covariate recognized as categorical.
```

```
## 
## Found 52 potential DMRs
```

You may want to look up what all of these settings in dmrFind are doing.  In the end, the method has discovered a list of DMRs and these are available in the 'dmrs' element of the output list:


```r
head(dmrs$dmrs)
```

```
##     chr    start      end     value     area  pns indexStart indexEnd nprobes        avg
## 321  22 45809244 45809952 -36.54717 548.2075 1984       5403     5417      15 -0.7857086
## 195  22 24890690 24891220 -32.98841 428.8494  705       1784     1796      13 -0.7530544
## 27   22 22901145 22902237  29.80695 387.4903  577       1437     1449      13  0.4918608
## 303  22 42896105 42896846 -26.58877 212.7102 1772       4815     4822       8 -0.6895311
## 182  22 21386798 21387058 -26.52531 212.2025  468       1164     1171       8 -0.6100558
## 354  22 50720468 50721626 -29.91718 179.5031 2538       6646     6651       6 -0.7536264
##            max  area.raw
## 321 -0.9720718 11.785628
## 195 -0.9441602  9.789707
## 27   0.8346143  6.394191
## 303 -0.9071111  5.516249
## 182 -0.9047373  4.880447
## 354 -0.9642824  4.521758
```

#### Exercise 3. For one of the DMRs found using 450k array data (as above), take a look in the corresponding RRBS data for the same region (assuming measurements are made) to verify that there is validating evidence of the DMR.


