##' Groups SNPs according to LD.
##' 
##' Merges the SNP IDs of those variants that fulfill the following criteria:
##' \itemize{
##' \item{LD (r2) >= than a given user-defined threshold.}
##' \item{Relative differences in sQTL F score <= 0.05}
##' \item{(if svQTL = TRUE) Relative differences in sQTL F score <= 0.05}}
##' 
##' @title Linkage disequilibrium (LD) filter
##' @param genotype.gene a data.frame with genotype info produced by 'read.bedix'.
##' @param tre.dist a distance object from the transcript relative expression.
##' @param com.samples a character vector with the common sample names between the genotype
##' and transcript relative expression files.
##' @param th r2 threshold over which SNPs will be merged. Default is 1.
##' @param svQTL \code{LD.filter} only merges the SNP IDs if the LD r2 is higher than
##' \code{th} and the F scores for the eQTL test differ no more than a 5\%. If svQTL = TRUE, 
##' it also requires that the F scores for the svQTL test differ no more than the 20\%. Default is TRUE.
##' @return a data.frame with columns:
##' \item{F}{the F score.}
##' \item{nb.groups}{the number of groups created by the genotypes.}
##' \item{md}{the maximum difference in splicing ratios between genotype groups.}
##' \item{tr.first, tr.second}{the two transcripts that change the most.}
##' \item{pv}{if \code{qform = TRUE} a P-value for the F score is computed}
##' \item{LD}{IDs of the SNPs that are in LD >= th, NA otherwise}
##' @author Diego Garrido-Martín
##' @keywords internal
LD.filter <- function(genotype.gene, com.samples, tre.dist, th = 1, tol = 0.05, svQTL = FALSE){
if (th > 1){
  stop ("Treshold for LD must be <= 1.")
}
ids <- genotype.gene$snpId
g <- genotype.gene[, com.samples]
colnames(g) <- rownames(g) <- NULL
nG <- apply(g, 1, function(x)(length(table(x)))) # Get nb of groups (2 or 3) for each SNP ID
M <- t(g)
M[M == -1] <- NA

res3 <- computeLD(M = M[, nG == 3], ids = ids[nG == 3], tre.dist = tre.dist, th = th, tol = tol, svQTL = svQTL)
res2 <- computeLD(M = M[, nG == 2], ids = ids[nG == 2], tre.dist = tre.dist, th = th, tol = tol, svQTL = svQTL)

if(is.null(res3) & !is.null(res2)){ # Some checks
  res <- res2
}else if (!is.null(res3) & is.null(res2)){
  res <- res3 
}else{
  res <- rbind(res3,res2)
}

if (is.null(res)){
  return(genotype.gene)
} else{
  res <- as.data.frame(res)
  colnames(res) <- "LD"
  res$LD <- as.character(unlist(res$LD))
  
  linked <- c(rownames(res), unlist(strsplit(res$LD, ", ", fixed = TRUE)))
  indep_names <- ids[!ids%in%linked]
  indep <- rep(NA, length(indep_names))
  names(indep) <- indep_names
  res <- rbind(res, data.frame(LD = indep))
  
  rownames(genotype.gene) <- ids
  genotype.gene <- merge(genotype.gene, res, by = "row.names") 
  genotype.gene <- with(genotype.gene, genotype.gene[order(chr, start), ])
  genotype.gene$Row.names <- NULL
  
  return(genotype.gene)
}
}

computeLD <- function(M, ids, tre.dist, th = 1, tol = 0.05, svQTL = FALSE){
  M <- as.matrix(M)
  if(ncol(M) > 1){
    s <- ncol(M)
    R <- cor(M, use = "pairwise.complete.obs")
    R <- R^2
    blocks <- lapply(as.data.frame(R),function(x) which(x >= th))
    names(blocks) <- 1:s 
    if (th < 1){
      blocks <- F.filter(blocks = blocks, tre.dist = tre.dist, M = M, tol = tol, svQTL = FALSE)
      if(svQTL){
        blocks <- F.filter(blocks = blocks, tre.dist = tre.dist, M = M, tol = tol, svQTL = TRUE)
      }
    }
    store <- c()
    bsc <- prune_reorder(blocks, R)
    blocks <- bsc[[1]]
    bsizes <- bsc[[2]]
    res <- list()
    while(length(bsizes) > 0){
      i <- names(bsizes[1])
      group <- blocks[[i]]
      res[[i]] <- group[group != as.numeric(i)]
      store <- unique(c(store, group)) 
      blocks[as.character(group)] <- NULL
      blocks <- lapply(blocks, function(x) x[!x%in%store])
      bsc <- prune_reorder(blocks, R)
      blocks <- bsc[[1]]
      bsizes <- bsc[[2]]
    }
    if(length(unlist(res)) > 0){
      res <- lapply(res, function(x) paste(ids[x], collapse=", "))
      names(res) <- ids[as.numeric(names(res))]
      res <- res[res != ""]
      return(as.matrix(res))
    }else{
      return(NULL)
    }
  }else{
    return(NULL)
  }
}

prune_reorder <- function (edges, R){
  nedges <- unlist(lapply(edges, length))
  nedges <- nedges[nedges > 1]
  if(length(nedges) == 0){
    return(list(edges,nedges))
  }
  edges <- edges[names(nedges)]
  
  redges <- c()
  for (e in names(edges)){
    subs <- edges[[e]]
    redges[e] <- mean(R[as.numeric(e), subs])
  }
  dd <- data.frame(redges,nedges)
  ordered <- rownames(dd[with(dd, order(-redges, -nedges)), ])
  nedges <- nedges[ordered]
  return(list(edges, nedges))
}

F.calc <- function(tre.dist, snp, svQTL = FALSE){
  
  if (any(is.na(snp))) {
    non.na <- !is.na(snp)
    snp <- snp[non.na]
    tre.dist <- stats::as.dist(as.matrix(tre.dist)[non.na, non.na])
  }
  snp.f <- factor(snp)
  if(!svQTL){
    G <- gower(tre.dist)
    dfnum <- nlevels(snp.f) - 1
    dfden <- attr(tre.dist,"Size") - dfnum - 1
    X <- stats::model.matrix(~., data = data.frame(genotype = snp.f), contrasts.arg = list("genotype" = "contr.sum"))     
    H <- tcrossprod(tcrossprod(X, solve(crossprod(X))), X)
    numer <- crossprod(c(H), c(G))
    trG <- sum(diag(G))
    denom <- trG - numer                                                          
    f.snp <- as.numeric((numer*dfden)/(denom*dfnum))  
    return(f.snp)
  }else{
    bd <- vegan::betadisper(tre.dist, snp.f, type = "centroid")
    bd.perm <- permutest.betadisper(bd, control = permute::how(nperm = 2)) 
    return(bd.perm$F)
  }
  
}

F.filter <- function(blocks, tre.dist, M, tol = 0.05, svQTL = FALSE){
  if(svQTL){
    tol <- 0.20
  }
  Fs <- apply(M, 2, function(x) F.calc(tre.dist = tre.dist, snp = x, svQTL = svQTL))
  names(Fs) <- 1:length(blocks)
  d <- list()
  for (e in names(blocks)){
    subs <- blocks[[e]]
    d[[e]] <- abs(Fs[as.character(subs)]-Fs[e])/Fs[e]
  }
  h <- list()
  for (k in names(blocks)){
    h[k] <- list(blocks[[k]] [ which(d[[k]] <= tol) ]) 
  }
  return(h)
}

