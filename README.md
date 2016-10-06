sQTLseekeR
==========

sQTLseekeR is a R package to detect splicing QTLs (sQTLs), which are variants associated with change in the splicing pattern of a gene. Here, splicing patterns are modeled by the relative expression of the transcripts of a gene.

For more information about the method and performance see article :
Monlong, J. et al. Identification of genetic variants associated with alternative splicing using sQTLseekeR. Nat. Commun.
5:4698 doi: [10.1038/ncomms5698](http://www.nature.com/ncomms/2014/140820/ncomms5698/full/ncomms5698.html) (2014).

### Installation

First some [Bioconductor](http://bioconductor.org/) packages are required. They can be installed with :

```
source("http://bioconductor.org/biocLite.R")
biocLite(c("vegan", "Rsamtools", "qvalue"))
```

Then, install [`devtools` package](https://github.com/hadley/devtools).

```
install.packages("devtools")
```

Finally, install either the latest development version:

```
devtools::install_github("guigolab/sQTLseekeR")
```

Or a specific release:

```
devtools::install_github("guigolab/sQTLseekeR", ref="2.1")
```

All these installation commands are also written in `install.R` script (*source* it to install). These commands work for R 3.2 or higher. For R 3.1, some packages need their older version (see `installOnR.3.1.R`).

**R 3.1 or higher** is required.

### Analysis steps

For an example analysis with the dummy data in [Data](https://github.com/guigolab/sQTLseekeR/tree/master/Data), see the [R-markdown overview](docs/sQTLseekeR-Overview.md).

The first step is to prepare the input data. `sQTLseekeR` requires three inputs:
* transcript expression. Column *trId* and *geneId*, corresponding to the transcript and gene ID are required. Then each column represents a sample and is filled with the expression values. Relative expression will be used hence both **read counts or RPKMs** works as the expression measure. However, it is not recommended to use transformed (log, quantile, or any non-linear transformation) counts or RPKMs because Hellinger distance is suited for relative expression.
* gene location information. In a BED-like format, the location of each gene is explicitly defined in this file.
* genotype information. The genotype of each sample is coded as follow: 0 for ref/ref; 1 for ref/mutated; 2 for mutated/mutated; -1 for missing value. Furthermore the first four columns should gather information about the SNP: *chr*, *start*, *end* and *snpId*. Finally **this file needs to be ordered** per *chr* and *start* position.

When all input files are correctly formatted `sQTLseekeR` prepares the data through functions `prepare.trans.exp` and `index.genotype`.
* `prepare.trans.exp` will :
  * remove transcripts with low expression.
  * remove genes with less than two expressed transcript.
  * remove genes with low splicing dispersion.
  * remove genes with not enough different splicing patterns.
  * flag samples with low gene expression.
* `index.genotype` compresses and indexes the genotype file to optimize further accession of particular regions. Note that the input file should not be compressed.

Once the input files are ready, `sqtl.seeker` function will compute the P-values for each pair of gene/SNP testing the association between the genotype and transcript relative expression. Here is a quick description of the parameters that would most likely be tweaked:
* `genic.window` the window(bp) around the gene in which the SNPs are tested. Default is 5000 (i.e. 5kb).
* `svQTL` should svQTLs test be performed in addition to sQTLs (default is FALSE). svQTLs are used to identify potential false positive among the significant sQTLs. svQTLs represents situation where the variance in transcript relative expression is different between genotype groups. In this particular situation identification of sQTLs is less robust as we assume homogeneity of the variance between groups, hence it might be safer to remove svQTLs from the list of reported sQTLs. However computation of svQTLs cannot rely on an asymptotic approximation, hence the heavy permutations will considerably increase the running time.
* `nb.perm.max` the maximum number of permutation/simulation to compute the P-value. The higher this number, the lower the P-values can potentially get but the longer the computation.

The output is a data.frame with for each gene/SNP:

+ *geneId*: the gene name
+ *snpId*: the SNP name
+ *F*: the F score
+ *nb.groups*: the number of genotype groups (2 or 3)
+ *md*: the maximum difference in relative expression between genotype groups
+ *tr.first*/*tr.second*: the transcript IDs of the two transcripts that change the most (and symetrically).
+ *pv*: the P-value
+ *nb.perms*: the number of permutation used for the P-value computation
+ *F.svQTL*/*pv.svQTL*/*nb.perms.svQTL*: idem for svQTLs, if `svQTL=TRUE` was used.

Finally, function `sqtls` is used to retrieve significant associations. The user can manually define a false discovery rate(FDR) or perform further filtering afterwards. Of note, there is a separate FDR threshold for svQTL removal (if svQTLs were computed), which is usually preferred to stay low (e.g. around 0.01).

An example of an analysis can be found in folder `scripts`.

### Running sQTLseekeR on computing clusters

`sQTLseekeR` can be used on a cluster using package `BatchJobs`. An example of an analysis using `BatchJobs` can be found in folder `scripts`.

`BatchJobs` is a potent package but basic functions are enough in our situation. Here is a quick practical summary of `BatchJobs` commands used in the script:
* `makeRegistry` creates a registry used to manipulate jobs for a particular analysis step.
* `batchMap` adds jobs to a registry. Simply, the user gives a function and a list of parameters. One job per parameter will be created to compute the output of the function using this specific parameter.
* `submitJobs` submits the jobs to the cluster. This is where the queue, maximum computation time, number of core can be specified. Moreover, if needed, a subset of the jobs can be sent to the cluster. Functions `findNotDone` and `findErrors` are particularly useful to find which jobs didn't finish or were lost in the limbo of the cluster management process.
* `showStatus` outputs the status of the computations.
* `loadResult` retrieves the output of one specific job, while `reduceResultsList` retrieves output for all jobs into a list format.

Another important point about `BatchJobs` is its configuration for the computing cluster. An example of the configuration files can be found in the `scripts` folder:
* If present in the working directory, `.BatchJobs.R` is loaded when the `BatchJobs` package is loaded. It defines which template to use and `BatchJobs` functions. In practice, it loads another R script file (here `makeClusterFunctionsAdaptive.R`) with the functions to use. In `.BatchJobs.R` users would only need to change the email address where to send the log messages to.
* In `makeClusterFunctionsAdaptive.R`, users just need to check/replace `qsub`/`qdel`/`qstat` calls with the correct bash commands (sometimes `msub`/`canceljob`/`showq`). This file should also be in the working directory when `BatchJobs` is loaded.
* Finally `cluster.tmpl` is a template form of a job bash script that would be send to the cluster. There the correct syntax for the resources or parameters of the cluster are defined. This file should also be in the working directory when `BatchJobs` is loaded.


### FAQ

#### `sqtl.seeker` outputs `NULL`. What am I doing wrong ?

An output of `NULL` means there were no gene/SNPs to analyze. It is likely due to inconsistent input files. To debug check that :

+ Gene IDs in the transcript expression and gene location are similar.
+ Genomic coordinates in the gene location and genotype information are consistent. E.g. *1* vs *chr1*.
+ Sample names in the transcript expression and genotype are similar.

If the input files are fine, an output of `NULL` might be caused by inappropriate transcript expression (e.g. genes with low expression) or genotypes (e.g. many missing values).

#### What are these svQTLs ?

svQTLs are SNPs associated with splicing variability of a gene. Here the relative transcript expression wight be globally similar between genotype groups but much more variable in specific one. Although the biological interpretation is not straightforward, we use them to flag potential false sQTLs. Indeed, the test for differential transcript relative expression assumes similar variability between genotype groups. Hence if a conservative approach to find sQTLs would be to retrieve significant sQTLs that are not svQTLs.

#### What about trans-sQTLs ?

By default, sQTLseekeR tests association between a gene and SNPs within or close by (defined by `genic.window` parameter). Testing association between all SNPs and all genes is not feasible : it would require too much computation and we don't see a good multiple-testing correction for this design that would fit our permuted/simulated P-values.

However, **the user can specify exactly which regions should be tested for each gene**. Instead of using the gene location for `gene.loc` parameter, the user can feed the locations of the regions to test. Similarly to the gene location, *chr*, *start*, *end* and *geneId* columns define the region and the gene to test. Several regions can be defined for a same gene. Of note, in this design, `genic.window=0` could be used to ensure that no flanking regions are added to the regions to test.

#### In the paper you associated a splicing event from the *tr.first*/*tr.second*. How ?

We can either use [AStalavista](http://genome.crg.es/astalavista/) or a build-in function in the package (`classify.events`).
For more infos, have a look at the [R-markdown tutorial on this](docs/SplicingEventClassification.md).

#### What is this "maximum difference" (*md*) column in the output ?

In addition to the F score and P-value, the maximum difference(MD) in relative expression between genotype groups is reported.
This is to be used as a measure of the size of the effect.
For example, if 'md' is 0.2, it means that there is one transcript whose relative expression shifted by 20% between two genotype groups.
Hence the higher the more dramatic the change in transcript usage.

#### For the same gene I get a several sQTLs but the transcripts (*tr.first*/*tr.second*) are different. Why ?

The primary result is the association between SNP(s) and a gene.
The two transcripts IDs in the output are here to provide some quick information about what changes the most but it's far from the complete picture and should be taken with a grain of salt.

If the same gene has several sQTLs with different *tr.first*/*tr.second* it could be the same changes affecting more than two transcripts.
Especially if the SNPs are close together or in LD, it is likely the same association but because of technical noise the most affected transcripts are not always the same two.

In practice, the first thing is to look at the transcript relative expression in the different genotype groups.
There you will see if it's mostly a pair of transcript affected or more.
Then, it's good to look at the structure of the transcripts.
Some transcripts might be very similar.
Or maybe they differ only in parts not affected by the SNP.
In this situation, it's not surprising for more than two transcripts to be affected.

