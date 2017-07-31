##' Retrieves sQTLs after two-level multiple testing correction and svQTL removal (if requested). 
##'
##' We consider two levels of multiple testing:
##' \itemize{
##' \item{1. Multiple genetic variants are tested per gene.}
##' \item{2. Multiple genes are tested genome-wide.}}
##' sQTLseekeR uses a permutation approach, implemented and detailed in \code{\link{sqtl.seeker.p}} to correct for
##' the former and false discovery rate (FDR) estimation to control for the latter. Both Benjamini-Hochberg (\code{\link{p.adjust}} method) 
##' and Storey (\code{\link{qvalue}} method) approaches for FDR can be applied.
##'
##' If \code{svQTL.removal = TRUE} and svQTLs were tested in 'sqtl.seeker', gene/SNPs with
##' significant svQTL association (after multiple testing correction and similar FDR threshold)
##' are removed from the final set of significant sQTLs.
##' @title sQTL retrieval from permuted phenotypes
##' @param res.nominal.df a data.frame, output of 'sqtl.seeker' with the nominal P-values
##' for each gene/SNP test.
##' @param res.permuted.df a data.frame, output of 'sqtl.seeker.p' with the empirical P-values
##' and the MLE estimates for the beta distribution parameters for each gene.
##' @param FDR the False Discovery Rate to call an association significant. Default is 0.05.
##' @param method the FDR approach. Either \code{"BH"} (Benjamini-Hochberg) or \code{"qvalue"} (Storey). The latter 
##' tends to be less stringent.  Default is "BH".
##' @param md.min the minimum MD (Maximum Difference) in relative expression required. Maximum
##' difference in relative expression (MD) gives an idea of the effect size of the association.
##' Default is 0.05.
##' @param svQTL.removal if TRUE (and column 'pv.svQTL' is present in 'res.df') significant
##' sQTL which are also significant svQTLs are not reported. Default is FALSE.
##' @param FDR.svQTL the False Discovery Rate to call a svQTL, that may be removed from the final set of sQTLs. 
##' Note that svQTL FDR is computed on the pooled nominal P-values from the svQTL test.
##' @return a data.frame with the columns of sqtl.seeker and sqtl.seeker.p outputs with the significant sQTLs.
##' @author Diego Garrido-Martín
##' @export
sqtls.p <- function(res.nominal.df, res.permuted.df, FDR = 0.05, method = "BH", md.min = 0.05, svQTL.removal = FALSE, FDR.svQTL = 0.05){
  
  pv <- md <- NULL ## Uglily suppress R checks for ggplot2
  traceBack <- function(nominals.gene, permuted.df, p_t){
    permuted.gene <- subset(permuted.df, geneId == nominals.gene$geneId[1])
    p_tn <- qbeta(p = p_t, shape1 = permuted.gene$shape1, shape2 = permuted.gene$shape2)
    permuted.gene$p_tn <- p_tn
    return(merge(permuted.gene, subset(nominals.gene, pv <= p_tn)))
  }
  if (any(colnames(res.nominal.df) == "pv.svQTL") & svQTL.removal) {
    if (method == "BH"){
      res.nominal.df$fdr.svQTL <- stats::p.adjust(res.nominal.df$pv.svQTL, method = "BH")
    }else if (method == "qvalue"){
      res.nominal.df$fdr.svQTL <- qvalue::qvalue(res.nominal.df$pv.svQTL)$qvalues
    }else{
      stop("Available methods for FDR are 'BH' and 'qvalue'.")
    }
  }
  if (method == "BH"){
    res.permuted.df$fdr <- stats::p.adjust(res.permuted.df$pv.emp.beta, method = "BH")
  }else if (method == "qvalue"){
    res.permuted.df$fdr <- qvalue::qvalue(res.permuted.df$pv.emp.beta)$qvalues
  }else{
    stop("Available methods for FDR are 'BH' and 'qvalue'.")
  }
  
  set0 <- subset(res.permuted.df, fdr <= FDR)
  set1 <- subset(res.permuted.df, fdr > FDR)
  p_t <- (sort(set1$pv.emp.beta)[1] - sort(-1.0 * set0$pv.emp.beta)[1]) / 2
  fdr_pt <- (sort(set1$fdr)[1] - sort(-1.0 * set0$fdr)[1]) / 2
  
  res.permuted.df <- set0
  res.nominal.df <- subset(res.nominal.df, geneId %in% res.permuted.df$geneId)
  
  err <- abs(fdr_pt - FDR)
  if (err > 0.5 * FDR){
    stop(paste0(sprintf("|closest observed FDR - FDR threshold set at %0.3f| > %0.3f.\n", FDR, 0.5 * FDR),
                "   Empirical P-value threshold cannot be estimated. Use sqtls function in this package instead."))
  } else if (err > 0.2 * FDR) {
      warning(sprintf("|closest observed FDR - FDR threshold set at %0.3f| > %0.3f.", FDR, 0.2 * FDR))
  } else {
      message(sprintf("|closest observed FDR - FDR threshold set at %0.3f| = %0.2e.", FDR, err))
  }
  message(sprintf("Global empirical P-value threshold = %0.2e", p_t))
  
  res.df <- as.data.frame(dplyr::do(dplyr::group_by(res.nominal.df, geneId), traceBack(., res.permuted.df, p_t)))
  
  if (svQTL.removal){
    res.df <- subset(res.df, fdr.svQTL > FDR.svQTL & md >= md.min)
  }else{
    res.df <- subset(res.df, md >= md.min)
  }
  return(res.df)
}
