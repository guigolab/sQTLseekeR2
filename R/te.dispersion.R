##' Compute the dispersion from a gene's transcript expression matrix. 
##' Dispersion is computed as the mean Hellinger distance to the centroid.
##' @title Dispersion computation
##' @param tr a transcript expression matrix 
##' @return a value for the dispersion
##' @author Diego Garrido-Martín
##' @keywords internal

te.dispersion = function (tr) {
  hellingerDist.p <- function(x1, x2) {
    a <- (sqrt(x1) - sqrt(x2))
    b <- sqrt(sum(a * a))
    return(b)
  }
  c <- as.numeric(apply(tr, 1, function(x)(mean(x, na.rm=T))))         
  d <- mean(apply(tr, 2, function(x)(hellingerDist.p(x, c))), na.rm=T)  
  return(d)
}