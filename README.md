# sQTLseekeR2.int

`sQTLseekeR2.int` is a slightly modified version of the `sQTLseekeR2` R package to detect splicing QTLs (sQTLs).

In contrast to `sQTLseekeR2`, `sQTLseekeR2.int` asesses a model with a second main effect (a 2-level factor such as gender, disease status, etc.), 
in addition to the genotype, and the interaction between the two of them. Hence, it allows to identify condition-biased sQTLs.

## Installation 

From R:

```r
install.packages("devtools")
devtools::install_github("guigolab/sQTLseekeR2@interaction")
```

**R 3.3 or higher** is required.

## Analysis steps

* Index the genotype file. `index.genotype` compresses and indexes the genotype file to optimize the access to particular regions.

* Preprocess the transcript expression data. `prepare.trans.exp` removes transcripts with low expression and genes expressing only
one transcript, with low splicing dispersion, with few different splicing patterns or low expression.
Then relative transcript expression is computed.

* Test for association. `sqtl.seeker` assesses the model: splicing ratios ~ genotype + condition + genotype:condition,
and computes a P-value per term, for each *cis* genetic variant-gene pair.

* Multiple testing correction available via [eigenMT](https://github.com/joed3/eigenMT) in sqtlseeker2.int-nf (see below)

## Running sQTLseekeR2.int on computing clusters (recommended)

`sQTLseekeR2.int` can be easily used on a cluster thanks to [Nextflow](https://www.nextflow.io). See [sqtlseeker2.int-nf](https://github.com/guigolab/sqtlseeker2-nf/tree/interaction) for details.

