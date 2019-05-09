# sQTLseekeR2

`sQTLseekeR2` is an enhanced version of [sQTLseekeR](https://github.com/guigolab/sQTLseekeR), an R package to detect splicing QTLs (sQTLs).

### Installation 

From R, source the [install.R](install.R) file.

```r
source("install.R")
```

**R 3.3 or higher** is required.

### Enhancements 

* Covariate correction
* Faster pre-processing of transcript expression
* Improved P-value calculation with higher precision
* Permutation-based schema for multiple testing correction 
* LD-based variant clustering 

For additional information on the analysis steps you can have a look at [sQTLseekeR](https://github.com/guigolab/sQTLseekeR).

### Running sQTLseekeR2 on computing clusters

`sQTLseekeR2` can be easily used on a cluster thanks to [Nextflow](https://www.nextflow.io). See [sqtlseeker2-nf](https://github.com/dgarrimar/sqtlseeker2-nf) for details.

