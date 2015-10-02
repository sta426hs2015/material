
# Part 1 - Brainstorm: Statistics

## Distributions
### e.g., Gaussian, Poisson, ..

## Statistical Models
## Methods for Estimation
## Methods for Hypothesis Testing

# Part 2 - Brainstorm: Technologies in Biology

## microarray, sequencing, etc.

## In-class exercise 1 (3 groups): microarray, Illumina seq, other types of seq

### Goal: 
#### produce a 1 paragraph summary
#### links to a few (<5) good resources
#### submit a pull request to brainstorm_answers.md

# Part 3 - Brainstorm: Applications/protocols in genomics OR subfields

### e.g., gene expression

# Part 4 - Brainstorm: Linking Technologies to Applications to Statistics (mainstream)

## e.g., microarray -> gene expression -> normally distributed (log intensities)

## In-class exercise 2 (3 groups): 
### Goal: make the link between technologies, what is measured, statistical models used
### For example:
##### RNA-seq: sequencing -> gene expression -> ?
##### BS-seq: bisulphite + sequencing -> DNA methylation -> ?
##### ChIP-seq: sequencing -> protein/DNA interactions -> ?


# Part 5 - Brainstorm: Methods/algorithms/data structures that are associated more to computer science

# Part 6 - Pick another "technology" (from those above or from [1]) to briefly describe

## In-class exercise 3 (individually): 
### Goal: 
#### write ~2 sentences about what the method does
#### again, make the link (technology -> application -> statistics)
#### submit a pull request to brainstorm_answers.md

[1] [https://liorpachter.wordpress.com/seq/](https://liorpachter.wordpress.com/seq/)

### Verena Steffen: Chem-Seq [(source)](http://www.nature.com/nbt/journal/v32/n1/pdf/nbt.2776.pdf)
Chem-Seq is a method to investigate small molecule/DNA interaction. It is used to gain insights of how a small molecule (chem. compound) influences cellular functions, for example, it could reveal the action mechanism of a DNA associated (influencing replication, transcription, etc.) drug.
The method is similar to ChIP-seq, which measures DNA/protein interactions, but instead of antibodies, chemical compounds are used. 

*(technology -> application -> statistics)*

**Chem-Seq -> small molecule/DNA interactions -> Poisson / Negative Binomial (NB) or zero-inflated NB**

### Natacha Bodenhausen: metagenomics [(source)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000352)

Reference: White JR, Nagarajan N, Pop M (2009) Statistical Methods for Detecting Differentially Abundant Features in Clinical Metagenomic Samples. PLoS Comput Biol 5(4): e1000352. doi: 10.1371/journal.pcbi.1000352 

Metagenomics is a method where all the DNA of a particular environment is sequenced. In comparative metagenomics, the goal is to identify features (taxa or biological pathway) which are more abundant in one populations of samples (for example, healthy individual) compared to another population of samples (for example, sick individuals).

*(technology -> application -> statistics)*

**DNA Sequencing -> Comparative metagenomics -> nonparametric t test (permutation method as described in Storey and Tibshirani)** 

### Marisa Silva: PD-Seq [(source)](http://www.pnas.org/content/110/24/E2153) 
Reference: Daniel Arango et al., “Molecular Basis for the Action of a Dietary Flavonoid Revealed by the Comprehensive Identification of Apigenin Human Targets,” Proceedings of the National Academy of Sciences 110, no. 24 (June 11, 2013): E2153–E2162, doi:10.1073/pnas.1303726110.

PD-Seq is a method that allow to understand protein-protein interaction. Phage display has been used to study protein-DNA, protein-protein and protein-peptide interactions. This technic combined with deep sequencing allows the analysis of high complex different interactions.

*(technology -> application -> statistics)* 

**PD-Seq-> small molecule-protein interaction ->  Binomial, Fisher's exact**


### Stefan Oberlin: GRO-Seq (Global run-on sequencing) [(source)](http://www.sciencemag.org/content/322/5909/1845.short) 
Reference: Leighton J. Core, Joshua J. Waterfall, John T. Lis: Nascent RNA Sequencing Reveals Widespread Pausing and Divergent Initiation at Human Promoters, Science (2008) 
 
In order to monitor global transcriptional activity, run-on technology was coupled to RNAseq. This technique allowed to generalise the phenomenon of stalled polymerases at the transcriptional start site (TSS) as well as their widespread bidirectional activity. 
 
*(technology -> application -> statistics)* 
 
**GRO-Seq-> transcriptional activity, eg. comparing activity at the TSS to gene body  ->  Fisher's exact**


