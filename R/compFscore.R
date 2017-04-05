##' Compute the F score, max diff ratio difference and transcripts that change the most.
##' Additionally, it can compute a P-value to assess the significance of the F score.
##' @title F score computation
##' @param geno.df a data.frame of one row with the genotype information for each sample.
##' @param tre.dist a distance object from the transcript relative expression.
##' @param tre.df a data.frame with the transcript relative expression. 
##' @param svQTL should svQTL test be performed in addition to sQTL. Default is FALSE.
##' @param qform should significance for the F score (sQTL test) be computed using 
##' the \code{\link[CompQuadForm]{davies}} method in the \code{CompQuadForm} package. 
##' Default is TRUE.
##' @return a data.frame with columns:
##' \item{F}{the F score.}
##' \item{nb.groups}{the number of groups created by the genotypes.}
##' \item{md}{the maximum difference in splicing ratios between genotype groups.}
##' \item{tr.first, tr.second}{the two transcripts that change the most.}
##' \item{pv}{if \code{qform = TRUE} a P-value for the F score is computed}
##' @author Diego Garrido-Martín, Jean Monlong
##' @keywords internal
##' @import CompQuadForm
compFscore <- function(geno.df, tre.dist, tre.df, svQTL = FALSE, qform = TRUE){
  
  if(class(tre.dist) != "dist"){
    stop("'tre.dist' must be a distance object.")
  }
  if(nrow(geno.df) > 1){
    stop(geno.df$snpId[1], " SNP is duplicated in the genotype file.")
  }
  # # if(!any(colnames(geno.df) %in% labels(tre.dist))){
  # #   stop("No common samples between genotype and transcript ratios files.")
  # # } # Checked before
  
  geno.snp <- as.numeric(geno.df[,labels(tre.dist)])
  names(geno.snp) <- labels(tre.dist)
  if(any(geno.snp == -1)){
    non.na <- geno.snp > -1
    geno.snp <- geno.snp[non.na]
    tre.dist <- stats::as.dist(as.matrix(tre.dist)[non.na, non.na])
  }
  groups.snp.f <- factor(as.numeric(geno.snp))
  mdt <- md.trans(tre.df, groups.snp.f, labels(tre.dist))
  if(qform){
    G <- gower(tre.dist)
    X <- stats::model.matrix(~., data = data.frame(genotype = groups.snp.f), contrasts.arg = list("genotype" = "contr.sum"))     
    p <- ncol(X) - 1
    n <- nrow(X)
    H <- tcrossprod(tcrossprod(X, solve(crossprod(X))), X)
    numer <- crossprod(c(H), c(G))
    trG <- sum(diag(G))
    denom <- trG - numer                                                          
    f.snp <- as.numeric(numer/denom)                                                  
    lambda <- eigen(G, only.values = T)$values                                    
    lambda <- abs(lambda[abs(lambda) > 1e-12 ])                                     
    start.acc <- 1e-14   
    item.acc <- start.acc
    pv.snp <- pcqf(q = f.snp, lambda = lambda, k = p, p = p, n = n, acc = start.acc)     
    while (length(pv.snp) > 1) {                                                      
      item.acc <- item.acc * 10
      pv.snp <- pcqf(q = f.snp, lambda = lambda, k = p, p = p, n = n, acc = item.acc)  
    }
    if(pv.snp < item.acc) {
      pv.snp <- item.acc
    }
    res.df <- data.frame(F = f.snp * (n - p - 1) / p,  
                         nb.groups = nlevels(groups.snp.f),
                         md = mdt$md,
                         tr.first = mdt$tr.first,
                         tr.second = mdt$tr.second,
                         pv = pv.snp,
                         stringsAsFactors = FALSE) # Here note that p = df(factor) - 1
  }else{
    F.snp <- adonis.comp(tre.dist, groups.snp.f, permutations = 2, svQTL = FALSE)
    res.df <- data.frame(F = F.snp,
                         nb.groups = nlevels(groups.snp.f),
                         md = mdt$md,
                         tr.first = mdt$tr.first,
                         tr.second = mdt$tr.second,
                         stringsAsFactors = FALSE)
  }
  if(svQTL){
    res.df$F.svQTL <- adonis.comp(tre.dist, groups.snp.f, permutations = 2, svQTL = TRUE)
  }
  if (any(colnames(geno.df) == "LD")){
    res.df$LD <- geno.df$LD
  }
  return(res.df)
}

gower <- function (d.mat) {
  d.mat <- as.matrix(d.mat)
  n <- nrow(d.mat)
  A <- -0.5 * d.mat^2
  As <- A - rep(colMeans(A), rep.int(n, n))
  return(t(As) - rep(rowMeans(As), rep.int(n, n)))
}

pcqf <- function(q, lambda, k, p, n = length(lambda), lim = 5e4, acc = start.acc) {
  gamma <- c(lambda, -q * lambda)                                               
  nu <- c(rep(k, length(lambda)), rep(n - p - 1, length(lambda)))               
  pv <- CompQuadForm::davies(0, lambda = gamma, h = nu, lim = lim, acc = acc)
  if (pv$ifault != 0) {                                                        
    return(pv)                                                                  
  }                                                                             
  if (pv$Qq < 0) {                                                              
    return(pv)
  }                                                                             
  if (pv$ifault == 0) {
    return(pv$Qq)                                    
  }
} 