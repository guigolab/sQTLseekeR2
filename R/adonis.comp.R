##' Compute either observed F score or permuted F scores.
##' @title Compute F score.
##' @param dis the distance object of the transcript relative expression.
##' @param groups a factor with the group information.
##' @param permutations the number of permutations.
##' @param f.perms should the permuted F scores be returned instead of the
##' real F score. Default is FALSE.
##' @param svQTL should svQTL test be performed instead of sQTL. Default is FALSE.
##' @param approx should the asymptotic distribution be used instead of permutations.
##' Default is TRUE.
##' @return a vector with the F score or the permuted F scores.
##' @author Jean Monlong, Diego Garrido-Martín
##' @keywords internal
adonis.comp <- function(tre.mt, groups, permutations = 99, f.perms = FALSE, svQTL = FALSE, approx = TRUE){
  if(svQTL){
    bd <- vegan::betadisper(dist(tre.mt), groups, type = "centroid")
    bd.perm <- permutest.betadisper(bd, control = permute::how(nperm = permutations)) 
    if(f.perms){
      return(bd.perm$f.perms)
    } else {
      return(bd.perm$F)
    }
  } else {
    if (approx) {
      return(approx.dist(tre.mt, permutations, groups))
    } else {
      res <- vegan::adonis(dist(tre.mt) ~ groups, permutations = permutations)
      return(as.numeric(res$f.perms[, 1]))
    }
  }
}

approx.dist <- function(Y, nb.mont, groups) {
  
  nb.gp <- nlevels(groups)
  n <- nrow(Y)
  fit <- lm(Y ~ groups)
  R <- fit$residuals
  Df <- nb.gp - 1
  Df.e <- fit$df.residual
  e <- eigen(cov(R)*(n-1)/Df.e, symmetric = T, only.values = T)$values
  eigenStats <- c(length(e), sum(e > 0), sum(e < 0))
  
  if (eigenStats[3] > 0) 
    e <- abs(e)
  
  randomChisqN <- matrix(stats::rchisq(nb.mont * 
                                         eigenStats[1], df = nb.gp - 1), nrow = eigenStats[1], 
                         ncol = nb.mont)
  randomChisqD <- matrix(stats::rchisq(nb.mont * 
                                         eigenStats[1], df = n - nb.gp), nrow = eigenStats[1], 
                         ncol = nb.mont)
  asymptNume <- e %*% randomChisqN
  asymptDeno <- e %*% randomChisqD
  asymptF <- asymptNume/asymptDeno * (n - nb.gp)/(nb.gp - 1)
  return(asymptF)
}