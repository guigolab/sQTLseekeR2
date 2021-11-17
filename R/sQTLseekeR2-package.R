##' sQTLseekeR2: R package for splicing QTL mapping
##'
##' sQTLseekeR2 tests for association between genetic variants and alternative
##' splicing (relative transcript abundances of a gene measured from RNA-seq), 
##' using a multivariate non-parametric approach. It is an enhanced version of sQTLseekeR
##' (see Details). 
##' 
##' Main enhancements in sQTLseekeR2: 
##' \itemize{
##' \item{Covariate correction}
##' \item{Faster pre-processing of transcript expression}
##' \item{Improved P-value calculation with higher precision}
##' \item{Permutation-based schema for multiple testing correction}
##' \item{LD-based variant clustering}
##' }
##' @author Diego Garrido-Martín, Jean Monlong, Miquel Calvo, Ferran Reverter, Roderic Guigó
##' @references Garrido-Martín, D. et al. Identification and analysis of splicing quantitative 
##' trait loci across multiple tissues in the human genome. Nat Commun 12, 727 (2021). 
##' \url{https://doi.org/10.1038/s41467-020-20578-2}
##' @docType package
##' @name sQTLseekeR2-package
##' @seealso Useful links:
##' \itemize{
##' \item{\url{https://github.com/guigolab/sQTLseekeR2}}
##' \item{\url{https://github.com/guigolab/sqtlseeker2-nf}}
##' }
##' @keywords internal
"_PACKAGE"

NULL