##' Compute the maximum difference in transcript relative expression
##' between genotype groups.
##' @title MaxDiff splicing ratios computation
##' @param sr.o a matrix with the transcript relative expression (samples x transcripts).
##' @param groups.o a factor with the genotype for each sample.
##' @return A list with:
##' \item{md}{the maximum difference in splicing ratios between genotype groups.}
##' \item{tr.first, tr.second}{the two transcripts that change the most.}
##' @author Jean Monlong, Diego Garrido-Martín
##' @keywords internal
md.trans <- function(sr.o, groups.o){
  mTrans <- apply(sr.o, 2, function(sr.r) tapply(sr.r, groups.o, mean, na.rm = TRUE))
  lr <- nrow(mTrans)
  ind1 <- rep(1:(lr-1), (lr-1):1)
  ind2 <- NULL
  for(ii in 2:lr){
    ind2 <- c(ind2, ii:lr)
  }
  MDtrans <- apply(mTrans, 2, function(r) diff(rbind(r[ind1], r[ind2])))
  if(!is.matrix(MDtrans)){
    MDtrans <- matrix(MDtrans, 1)
  }
  gpMD <- apply(MDtrans, 1, function(e) max(abs(e)))
  gpMD.max <- which.max(gpMD)
  tr.first <- names(which.max(abs(MDtrans[gpMD.max, ])))
  tr.second <- names(which.max(-sign(MDtrans[gpMD.max, tr.first]) * MDtrans[gpMD.max, ]))
  return(list(md = max(gpMD, na.rm = TRUE), tr.first = tr.first, tr.second = tr.second))
}
