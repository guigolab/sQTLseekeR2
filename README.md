# sQTLseekeR2

[![R-CMD-check](https://github.com/guigolab/sQTLseekeR2/workflows/R-CMD-check/badge.svg)](https://github.com/guigolab/sQTLseekeR2/actions)

`sQTLseekeR2` is an enhanced version of [sQTLseekeR](https://github.com/guigolab/sQTLseekeR), an R package to detect splicing QTLs (sQTLs).

## Installation 

From R:

```r
install.packages("devtools")
devtools::install_github("guigolab/sQTLseekeR2")
```

**R 3.3 or higher** is required.

## Enhancements 

* Covariate correction
* Faster pre-processing of transcript expression
* Improved P-value calculation with higher precision
* Permutation-based schema for multiple testing correction 
* LD-based variant clustering 

## Analysis steps

* Index the genotype file. `index.genotype` compresses and indexes the genotype file to optimize the access to particular regions.

* Preprocess the transcript expression data. `prepare.trans.exp` removes transcripts with low expression and genes expressing only
one transcript, with low splicing dispersion, with few different splicing patterns or low expression.
Then relative transcript expression is computed.

* Test for association between splicing ratios and genetic variants in *cis* (nominal pass). `sqtl.seeker` computes a nominal P-value for
each variant-gene pair, testing for the association between the genotype and the transcript relative expression.

* Obtain an empirical P-value for each phenotype (permutation pass, optional). `sqtl.seeker.p` implements a permutation scheme that 
empirically characterizes, for each gene, the distribution of nominal P-values expected under the null hypothesis of no association. 
This null distribution is then modeled using a beta distribution as in [FastQTL](https://academic.oup.com/bioinformatics/article/32/10/1479/1742545).

* Control for multiple testing. `sqtls` computes FDR (Benjamini-Hochberg's or Storey's) across all nominal tests. 
Alternatively, if `sqtl.seeker.p` has been run, `sqtls.p` performs FDR on empirical P-values and then,
to recover all significant variant-gene pairs, implements a procedure identical to the one depicted [here](https://media.nature.com/original/nature-assets/nature/journal/v550/n7675/extref/nature24277-s1.pdf).

For additional information on the analysis steps you can have a look at [sQTLseekeR](https://github.com/guigolab/sQTLseekeR).

## Running sQTLseekeR2 on computing clusters (recommended)

`sQTLseekeR2` can be easily used on a cluster thanks to [Nextflow](https://www.nextflow.io). See [sqtlseeker2-nf](https://github.com/guigolab/sqtlseeker2-nf) for details.

## Cite sQTLseekeR2

If you find `sQTLseekeR2` useful in your research please cite the related publication:

Garrido-Martín, D., Borsari, B., Calvo, M., Reverter, F., Guigó, R. Identification and analysis of splicing quantitative trait loci across multiple tissues in the human genome. *Nat Commun* 12, 727 (2021). [https://doi.org/10.1038/s41467-020-20578-2](https://doi.org/10.1038/s41467-020-20578-2)


