##' Retrieves sQTLs after multiple testing correction and svQTL removal (if requested).
##' The distribution of the P-values and a semi-volcano plot can also be shown.
##'
##' Multiple testing correction can be performed using either Benjamini-Hochberg (\code{\link{p.adjust}} method) 
##' or Storey (\code{\link{qvalue}} method) approach for FDR. 
##'
##' If \code{svQTL.removal = TRUE} and svQTLs were tested in 'sqtl.seeker', gene/SNPs with
##' significant svQTL association (after multiple testing correction and similar FDR threshold)
##' are removed from the final set of significant sQTLs.
##' @title sQTL retrieval
##' @param res.df a data.frame, output of 'sqtl.seeker' with the P-values
##' for each gene/SNP test.
##' @param FDR the False Discovery Rate to call an association significant. Default is 0.05.
##' @param method FDR approach: either \code{"BH"} (Benjamini-Hochberg) or \code{"qvalue"} (Storey).
##' Default is "BH".
##' @param md.min the minimum MD (Maximum Difference) in relative expression desired. Maximum
##' difference in relative expression (MD) gives an idea of the effect size of the association.
##' Default is 0.05.
##' @param out.pdf the name of the pdf file to create. If NULL (default), no pdf output
##' will be created. If non-NULL, the distribution of the P-values and a semi-volcano plot
##' (P-value vs MD) will be shown.
##' @param svQTL.removal if TRUE (and column 'pv.svQTL' is present in 'res.df') significant
##' sQTL which are also significant svQTLs are not reported.
##' @param FDR.svQTL the False Discovery Rate to call a svQTL, that may be removed from the final set of sQTLs.
##' @return a subset of the input data.frame with only significant sQTLs and FDR estimates.
##' @author Jean Monlong, Diego Garrido-Martín
##' @export
sqtls <- function(res.df, FDR = 0.05, method = "BH", md.min = 0.05, out.pdf = NULL, svQTL.removal = TRUE, FDR.svQTL = 0.05){
  pv <- md <- NULL ## Uglily suppress R checks for ggplot2
  
  if (method == "BH"){
    res.df$fdr <- stats::p.adjust(res.df$pv)
  }else if (method == "qvalue"){
    res.df$fdr <- qvalue::qvalue(res.df$pv)$qvalues
  }else{
    stop("Available methods for FDR are 'BH' and 'qvalue'.")
  }
  
  if(!is.null(out.pdf)){
    grDevices::pdf(out.pdf, 8, 6)
    suppressWarnings(print(ggplot2::ggplot(res.df, ggplot2::aes(x = pv)) +
          ggplot2::geom_histogram() + ggplot2::theme_bw() + 
          ggplot2::xlab("P-value") + ggplot2::ylab("number of gene/SNP pairs")
          ))
  }
  
  if(any(colnames(res.df)=="pv.svQTL")){
    if (method == "BH"){
      res.df$fdr.svQTL <- stats::p.adjust(res.df$pv.svQTL)
    }else if (method == "qvalue"){
      res.df$fdr.svQTL <- qvalue::qvalue(res.df$pv.svQTL)$qvalues
    }
    if(svQTL.removal){
      res.df <- res.df[which(res.df$fdr.svQTL >= FDR.svQTL), ]
      if(!is.null(out.pdf)){
        suppressWarnings(print(ggplot2::ggplot(res.df, ggplot2::aes(x = pv)) +
              ggplot2::geom_histogram() + ggplot2::theme_bw() +
              ggplot2::xlab("P-value") + ggplot2::ylab("number of gene/SNP pairs") +
              ggplot2::ggtitle("After svQTL removal")
              ))
      }
    }
  }

  res.df = res.df[which(res.df$fdr <= FDR & res.df$md >= md.min), ]
  rownames(res.df) = NULL

  if(!is.null(out.pdf)){
    if(nrow(res.df)>0){
      suppressWarnings(print(ggplot2::ggplot(res.df, ggplot2::aes(y = -log10(pv), x = md)) +
            ggplot2::geom_bin2d(bins = 200) + ggplot2::theme_bw() +
            ggplot2::ylab(expression('-log'[10]*' P-value')) + ggplot2::xlab("MD (Maximum difference in relative expression)")
            ))
    }
    grDevices::dev.off()
  }
  return(res.df)
}
