##' Retrieves sQTLs after two-level multiple testing correction and svQTL removal (if requested). The distribution of the P-values and 
##' a semi-volcano plot can be also displayed.
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
##' @param out.pdf the name of the pdf file to create. If NULL (default), no pdf output
##' will be created. If non-NULL, the distribution of the P-values and a semi-volcano plot
##' (P-value vs MD) will be shown.
##' @author Diego Garrido-Martín
##' @export
sqtls.p <- function(res.nominal.df, res.permuted.df, FDR = 0.05, method = "BH", md.min = 0.05, out.pdf = NULL, 
                    svQTL.removal = FALSE, FDR.svQTL = 0.05){
 
  pv <- md <- NULL ## Uglily suppress R checks for ggplot2
  traceBack <- function(nominals.gene, permuted.df, p_t){
    permuted.gene <- subset(permuted.df, geneId == nominals.gene$geneId[1])
    p_tn <- qbeta(p = p_t, shape1 = permuted.gene$shape1, shape2 = permuted.gene$shape2)
    permuted.gene$p_tn <- p_tn
    return(merge(permuted.gene, subset(nominals.gene, pv <= p_tn)))
  }
  
  if (any(colnames(res.nominal.df) == "pv.svQTL") & svQTL.removal) {
    if (method == "BH"){
      res.nominal.df$fdr.svQTL <- stats::p.adjust(res.nominal.df$pv.svQTL)
    }else if (method == "qvalue"){
      res.nominal.df$fdr.svQTL <- qvalue::qvalue(res.nominal.df$pv.svQTL)$qvalues
    }else{
      stop("Available methods for FDR are 'BH' and 'qvalue'.")
    }
  }
  
  if (method == "BH"){
    res.permuted.df$fdr <- stats::p.adjust(res.permuted.df$pv.emp.beta)
  }else if (method == "qvalue"){
    res.permuted.df$fdr <- qvalue::qvalue(res.permuted.df$pv.emp.beta)$qvalues
  }else{
    stop("Available methods for FDR are 'BH' and 'qvalue'.")
  }
  
  res.permuted.df <- subset(res.permuted.df, fdr <= FDR)
  res.nominal.df <- subset(res.nominal.df, geneId %in% res.permuted.df$geneId)
  err <- min(abs(res.permuted.df$fdr - FDR))
  if (err > 0.5 * FDR){
    stop(paste0(sprintf("|closest observed FDR - FDR threshold set at %0.3f| > %0.3f.\n", FDR, 0.5 * FDR),
                "   Empirical P-value threshold cannot be estimated. Use sqtls function in this package instead."))
  } else if (err > 0.2 * FDR) {
      warning(sprintf("|closest observed FDR - FDR threshold set at %0.3f| > %0.3f.", FDR, 0.2 * FDR))
  }
  p_t <- res.permuted.df[which.min(abs(res.permuted.df$fdr - FDR)), "pv.emp.beta"] 
  message(sprintf("Global empirical P-value threshold = %0.2e",p_t))
  res.df <- as.data.frame(dplyr::do(dplyr::group_by(res.nominal.df, geneId), traceBack(., res.permuted.df, p_t)))
  
  if (!svQTL.removal){
    res.df <- subset(res.df, fdr.svQTL > FDR.svQTL & md >= md.min)
    if(!is.null(out.pdf)){
      grDevices::pdf(out.pdf, 8, 6)
      suppressWarnings(print(ggplot2::ggplot(res.df, ggplot2::aes(x = pv)) +
                               ggplot2::geom_histogram() + ggplot2::theme_bw() + 
                               ggplot2::xlab("P-value") + ggplot2::ylab("number of gene/SNP pairs")
      ))
      grDevices::dev.off()
    }
  }else{
    res.df <- subset(res.df, md >= md.min)
    if(!is.null(out.pdf)){
      grDevices::pdf(out.pdf, 8, 6)
      suppressWarnings(print(ggplot2::ggplot(res.df, ggplot2::aes(x = pv)) +
                               ggplot2::geom_histogram() + ggplot2::theme_bw() + 
                               ggplot2::xlab("P-value") + ggplot2::ylab("number of gene/SNP pairs") +
                               ggplot2::ggtitle("After svQTL removal")
      ))
      suppressWarnings(print(ggplot2::ggplot(res.df, ggplot2::aes(y = -log10(pv), x = md)) +
                               ggplot2::geom_bin2d(bins = 200) + ggplot2::theme_bw() +
                               ggplot2::ylab(expression('-log'[10]*' P-value')) + 
                               ggplot2::xlab("MD (Maximum difference in relative expression)")
      ))
      grDevices::dev.off()
    }
  } 
  return(res.df)
}
