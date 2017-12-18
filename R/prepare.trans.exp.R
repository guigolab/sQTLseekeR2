##' Prepare the transcript expression matrix. Transcripts with low
##' expression are removed. Genes expressing only one transcript
##' are removed because they are not informative for splicing analysis. 
##' The relative transcript expression is retrieved and genes with
##' low variability in splicing are removed. Also genes with a small number (<25) of
##' different splicing patterns are removed as they don't fit the
##' permutation process.
##'
##' After removal of inappropriate transcripts/genes, transcript
##' relative expression is computed by dividing its expression
##' by the total expression of all transcripts for the gene.
##' @title Transcript expression preparation
##' @param te.df a data.frame with the transcript expression. The first
##' two columns are the transcript and gene ids, named 'trId' and
##' 'geneId'. The rest of the columns contain the expression values for each sample.
##' @param min.transcript.exp the minimum transcript expression. Transcripts
##' with lower expression in all samples are removed.
##' @param min.gene.exp the minimum gene expression. Samples with lower
##' gene expression are removed from the analysis of this gene. If after this step
##' the proportion of retained samples for a given gene is smaller than 
##' \code{min.prop} (see below), the gene is removed from the analysis.
##' @param min.prop the minimum proportion of samples with gene expression above
##' \code{min.gene.exp} for the gene to be analyzed.
##' @param min.dispersion the minimum dispersion of the transcript
##' relative expression. Genes with lower dispersion are removed.
##' @param verbose If TRUE the names of filtered genes and transcripts will be displayed. 
##' Default is FALSE.
##' @return a data.frame with the relative transcript expression for the
##' genes to study.
##' @author Diego Garrido-Martín, Jean Monlong
##' @export
prepare.trans.exp <- function(te.df, min.transcript.exp = 0.01, min.gene.exp = 0.01, min.prop = 0.8, min.dispersion = 0.1, verbose = FALSE){
  if(!all(c("geneId", "trId") %in% colnames(te.df))){
    stop("Missing column in 'te.df' : 'geneId' and 'trId' are required.")
  }
  if(min.transcript.exp < 0 || min.gene.exp < 0){
    stop("transcript/gene minimum expression should be greater or equal to 0.")
  }
  if(min.prop < 0 || min.prop > 1){
    stop("'min.prop' value should be in [0, 1].")
  }
  
  ## Convert into character, just in case
  te.df$geneId <- as.character(te.df$geneId)
  te.df$trId <- as.character(te.df$trId)
  ##

  samples <- setdiff(colnames(te.df), c("chr", "start", "end", "geneId", "trId"))
  if(length(samples) < 5){
    stop("Not enough samples; at least 5 samples required (although at least 20 is recommended).")
  }
  if(length(samples) < 20){
    warning("Low sample size : it's recommended to have at least 20 samples.")
  }
  trans.to.keep <- apply(te.df[, samples], 1, function(r) any(r >= min.transcript.exp))
  if(all(!trans.to.keep)){
    stop("No transcripts with expression above threshold")
  }
  if(verbose && any(!trans.to.keep)){
    message("Filtered transcripts due to low expression : ", paste(te.df$trId[which(!trans.to.keep)], collapse = " "))
  }
  te.df <- te.df[which(trans.to.keep), ]
  nb.trans <- table(te.df$geneId)
  trans2 <- names(which(nb.trans > 1))
  if(verbose && any(nb.trans <= 1)){
    message("Filtered single-transcript genes : ", paste(setdiff(unique(te.df$geneId), trans2), collapse = " "))
  }
  te.df <- te.df[which(te.df$geneId %in% trans2), ]
  
  relativize.filter.dispersion <- function(df) {
    df[, samples] <- apply(df[, samples], 2, relativize, min.gene.exp = min.gene.exp)
    disp <- te.dispersion(df[, samples])
    non.na.p <- sum(!is.na(df[1, samples])) / length(samples) 
    if (disp >= min.dispersion && non.na.p >= min.prop  && nbDiffPt(df[, samples]) >=                  
        min(25, length(samples) * 0.8)){
      return(df)
    } 
    else {return(data.frame())}
  }
  
  te.df <- plyr::ldply(lapply(unique(te.df$geneId), function(gene.i){
    df <- te.df[which(te.df$geneId == gene.i), ]
    relativize.filter.dispersion(df)
  }), identity)

  if(verbose && length(unique(te.df$geneId)) != length(trans2)){
    message("Filtered low exp/disp genes : ", paste(setdiff(trans2,unique(te.df$geneId)), collapse = " "))
  }

  if(nrow(te.df) == 0){
    stop("No genes found with suitable transcript expression.")
  }

  return(te.df)
}

##' Compute the relative expression of transcripts.
##' @title Relative expression computation
##' @param x a vector of transcript expression.
##' @param min.gene.exp the minimum gene expression in the sample. If the
##' gene expression is too low it is safer to remove the sample from the
##' analysis.
##' @return a vector with the relative expression.
##' @author Jean Monlong
##' @keywords internal
relativize <- function(x, min.gene.exp = 0.01){
  x <- as.numeric(x)
  if (!any(is.na(x)) && sum(x) > min.gene.exp) {
    x/sum(x)
  } else {
    return(rep(NA, length(x)))
  }
}

##' Compute the dispersion from a gene's transcript expression matrix. 
##' Dispersion is computed as the mean Hellinger distance to the centroid.
##' @title Splicing dispersion computation
##' @param tr a data.frame with the splicing ratios (transcript x sample). 
##' @return a value for the dispersion.
##' @author Diego Garrido-Martín
##' @keywords internal
te.dispersion <- function (tr) {
  c <- as.numeric(apply(tr, 1, function(x)(mean(x, na.rm = T))))         
  d <- mean(apply(tr, 2, function(x)(hellingerDist.p(x, c))), na.rm = T)  
  return(d)
}

##' Compute the number of non-identical splicing ratios. Too few splicing ratios are
##' problematic for the permutation process. If many samples have exactly the same
##' splicing ratios, the permutation might create many identical permuted F scores.
##' @title Number of non-identical splicing ratio points
##' @param sr a data.frame with the splicing ratios (transcript x sample).
##' @return the number of non-identical splicing ratios.
##' @author Jean Monlong
##' @keywords internal
nbDiffPt <- function(sr){
  length(unique(apply(sr, 2, paste, collapse = "-")))
}
