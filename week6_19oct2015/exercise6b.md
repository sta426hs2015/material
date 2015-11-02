---
title: "Exercise 6b"
output: html_document
---


Setup: if you want to repeat this exercise on your own computer, you will/may need to install the following software packages: 

blat: https://genome.ucsc.edu/FAQ/FAQblat.html#blat3 

STAR: https://code.google.com/p/rna-star/ 

samtools: http://sourceforge.net/projects/samtools/files/samtools/ 

The files we will work with correspond to: 

**hs_ch19_subset.fa**  - human genome sequence from chr19 from 2,000,000-4,000,000, treat this as a "small genome" to speed things up 

**galaxy_small.fastq** - subset of RNA-seq reads from a real dataset 

**chr19_rescaled.gff** - table of annotation (helps the aligner with cDNA reads)


#### Part 1: This exercise is just to get you familiar with a tool called **BLAT**

However, raw reads (RNA-seq, DNA-seq, etc.) typically come as FASTQ files, but BLAT pre-dates high-throughput sequencing, before quality scores were used regularly and does not use quality information.  So, we first need to extract a FASTA file from the FASTQ file.

There are many tools that could do this, but it may be useful later to know how to manipulate sequence files in R, so we can do something like the following:


```r
library("ShortRead")
library("Biostrings")

rfq <- readFastq("galaxy_small.fastq")
class(rfq)

seqs <- sread(rfq)
names(seqs) <- id(rfq)
seqs  # a DNAStringSet object - ?"DNAStringSet-class"

## Use writeXStringSet() to make a FASTA file of reads -- read
# ?writeXStringSet on how to use it
```

You can run BLAT from the Terminal:

```
blat hs_ch19_subset.fa file_of_reads.fa output.psl
```

This creates a [PSL output file](http://www.ensembl.org/info/website/upload/psl.html).  Look up the details if you are unsure.

Note also that if you want to manually check some reads/mapping, you can use the online version of blat: http://genome.ucsc.edu/cgi-bin/hgBlat (make sure you pick the right reference genome).

Using R or unix command line tools, how many reads do you have to start with ?  How many are successfully aligned?  


#### Part 2: This exercise is to get you familiar with a commonly-used splice-aware aligner called **STAR**

For more details, see: [STAR github repo](https://github.com/alexdobin/STAR).

A standard run of STAR includes the following steps:

```
# ---------------------------------------------------------------
# Step 1: build an index (only need to do this once)
# ---------------------------------------------------------------
mkdir star  # make directory to store index
STAR --runMode genomeGenerate --genomeDir ./star  --genomeFastaFiles hs_ch19_subset.fa

# ---------------------------------------------------------------
# Step 2: map reads against reference genome, giving the known annotation
# ---------------------------------------------------------------
STAR --outFileNamePrefix STAR_mapped_ --genomeDir ./star --readFilesIn galaxy_small.fastq --sjdbGTFfile chr19_rescaled.gtf 

# ---------------------------------------------------------------
# Step 3: check output report
# ---------------------------------------------------------------
Look at [..]Log.final.out file to get some mapping statistics

# ---------------------------------------------------------------
# Step 4: do the dance for the alignments: create binary version, sort, index 
# ---------------------------------------------------------------
samtools view -bS STAR_mapped_Aligned.out.sam -o STAR_mapped.bam
samtools sort STAR_mapped.bam STAR_mapped_s
samtools index STAR_mapped_s.bam

# Perhaps remove some of the extra files
rm -i STAR_mapped_Aligned.out.sam STAR_mapped.bam
```


Also, a quick print out to the screen of the alignments can be achieved with a command like this (details of the columns can be found in the [SAM specification](http://samtools.github.io/hts-specs/SAMv1.pdf): 
```
samtools view STAR_mapped_s.bam | more 
```

How many alignments were recorded in the SAM/BAM file?  And, how many unique reads have alignments?


#### Part 3: Assuming the creation of the BAM worked fine, it is generally useful to visualize the mapped reads in a genome browser, such as [IGV](https://www.broadinstitute.org/igv/).  For IGV, you can follow these steps:

Genomes .. Load Genome from File .. select hs_ch19_subset.fa (for "genome") 

File .. Load from File .. [select sorted, indexed BAM file] (for mapped reads) 

File .. Load from File .. [select chr19_rescaled.gtf]  (for annotation) 

Create a screenshot of a gene that attracted reads.


#### Part 4: Reading the mapped reads in R

The Bioconductor project has significant infrastructure for processing various file formats, including BAM files.  For example:

 Run commands similar to:

```r
library("GenomicAlignments")
rga <- readGAlignments("STAR_mapped_s.bam")
cs <- strsplit( cigar(rga), "M|N" )
```

From the CIGAR strings, you can understand how many reads have gaps. How many reads have a single gap?  How many have 0 or 2 gaps?
