
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

