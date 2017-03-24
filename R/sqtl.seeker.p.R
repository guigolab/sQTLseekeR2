##' \code{sqtl.seeker.p} description
##' 
##' Details
##' 
##' @title permuted sQTL seeker 
##' @param tre.df a data.frame with transcript relative expression
##' produced by 'prepare.trans.exp'. Same as in \code{sqtl.seeker}.
##' @param genotype.f the name of the genotype file. This file needs to
##' be ordered by position, compressed and indexed using 'index.genotype' or externally using tabix (samtools).
##' Must have column 'snpId'. Same as in \code{sqtl.seeker}.
##' @param gene.loc a data.frame with the genes location. Columns 'chr', 'start',
##' 'end' and 'geneId' are required. Same as in \code{sqtl.seeker}.
##' @param genic.window the window(bp) around the gene in which the SNPs are tested. Default is 5000 (i.e. 5kb).
##' Same as in \code{sqtl.seeker}.
##' @param min.nb.ext.scores the minimum number of permuted  nominal P-values lower than
##' the lowest observed nominal P-value to allow the computation to stop. Default is 100. 
##' @param nb.perm.max the maximum number of permutations. Default is 1000. 
##' @param verbose Default is TRUE. Mainly for debugging.
##' @return a data.frame with columns:
##' \item{geneId}{the gene name}
##' \item{variants.cis}{the number of variants tested in cis.}
##' \item{LD}{a linkage disequilibrium estimate for the genomic window (median r2).}
##' \item{best.snp}{ID of the SNP with the smallest observed nominal P-value.}
##' \item{best.nominal.pv}{P-value corresponding to the best SNP.}
##' \item{shape1}{Beta distribution parameter shape1.}
##' \item{shape2}{Beta distribution parameter shape2.}
##' \item{nb.perms}{the number of permutations used for the empirical P-value computation.}
##' \item{pv.emp}{empirical P-value based on permutations.}
##' \item{pv.emp.beta}{empirical P-value based on the beta approximation.}
##' \item{runtime}{approximated computation time per gene.}
##' @author Diego Garrido-Martín
##' @export
##' @import fitdistrplus
sqtl.seeker.p <- function(tre.df, genotype.f, gene.loc, genic.window=5e3, min.nb.ext.scores=1e2, nb.perm.max=1e3, verbose=TRUE){
  
  . <- nb.groups <- snpId <- NULL ## Uglily appease R checks (dplyr)
  
  analyze.gene.f <- function(tre.gene){
    if(verbose) message(tre.gene$geneId[1])
    if(sum(duplicated(gene.loc$geneId)) > 1){
      stop(tre.gene$geneId[1], " Repeated gene in gene location file.")
    }
    ## Load genotype
    gr.gene <- with(gene.loc[which(gene.loc$geneId == tre.gene$geneId[1]), ], GenomicRanges::GRanges(chr, IRanges::IRanges(start, end)))
    if(genic.window > 0){
      gr.gene <- GenomicRanges::resize(gr.gene, GenomicRanges::width(gr.gene) + 2 * genic.window, fix = "center")
    }
    if(length(gr.gene) > 0){
      ## Remove samples with non expressed genes
      tre.gene <- tre.gene[, !is.na(tre.gene[1, ])]
      ## Focus on common samples
      genotype.headers <- as.character(utils::read.table(genotype.f, as.is = TRUE, nrows = 1))
      com.samples <- intersect(colnames(tre.gene), genotype.headers)
      if(length(com.samples) == 0){
        stop("No common samples between genotype and transcript ratios files.")
      }
      tre.dist <- hellingerDist(tre.gene[, com.samples])
      res.df <- data.frame()
      genotype.gene <- read.bedix(genotype.f, gr.gene) 
      if (verbose & is.null(genotype.gene)) {
        message("\tNo SNPs in the genomic range.")
      }
      if (!is.null(genotype.gene)) {
        snps.to.keep <- check.genotype(genotype.gene[, com.samples], tre.gene[, com.samples])
        if (verbose) {
          snps.to.keep.t <- table(snps.to.keep)
          message("\t", paste(names(snps.to.keep.t), snps.to.keep.t, sep = ": ", collapse = ", "))
        }
        if (any(snps.to.keep == "PASS")) {
          genotype.gene <- genotype.gene[snps.to.keep == "PASS", ]
          observed <- dplyr::do(dplyr::group_by(genotype.gene, snpId), compute.nominal.pv(., tre.dist, permute = FALSE))
          min.pv.obs <- min(observed$pv.snp) 
          best.snp <- observed$snpId[which.min(observed$pv.snp)]
          res.df <- compute.empirical.pv(genotype.gene = genotype.gene, tre.dist = tre.dist, min.nb.ext.scores = min.nb.ext.scores,
                                         nb.perm.max = nb.perm.max, min.pv.obs = min.pv.obs, best.snp = best.snp, verbose = verbose)
          return(data.frame(done = TRUE, res.df))
        }
      }
    } else {
      if (verbose) {
        warning("Issue with the gene location.")
      }
    }
    return(data.frame(done = FALSE))
  }

  ret.df <- lapply(unique(tre.df$geneId), function(gene.i){
    df <- tre.df[which(tre.df$geneId == gene.i), ]
    data.frame(geneId = gene.i, analyze.gene.f(df))
  })
  done <- which(unlist(lapply(ret.df, ncol)) > 2)
  if(length(done) > 0){
    ret.df <- ret.df[done]
    ret.df <- do.call(rbind, ret.df)
    ret.df$done <- NULL
    return(ret.df)
  } else {
    return(NULL)
  }
}

##' Computes nominal P-value for the association between the genotype of a SNP and the splicing ratios of a gene.
##' @title Compute nominal P-value
##' @param geno.df a data.frame of a single row corresponding to the genotype of one SNP 
##' observed in several individuals (columns).
##' @param tre.dist a interdistance matrix produced by \code{\link{hellingerDist}}.
##' @param permute should the interdistance matrix be permuted. Default is FALSE.
##' @return A data.frame containing a P-value for the association.
##' @author Diego Garrido-Martín 
##' @keywords internal
compute.nominal.pv <- function(geno.df, tre.dist, permute = FALSE){
  if (nrow(geno.df) > 1) {
    stop(geno.df$snpId[1], " SNP is duplicated in the genotype file.")
  } 
  geno.snp <- as.numeric(geno.df[, labels(tre.dist)])
  names(geno.snp) <- labels(tre.dist)
  D <- as.matrix(tre.dist)
  if (any(geno.snp == -1)) {
    non.na <- geno.snp > -1
    geno.snp <- geno.snp[non.na]
    D <- D[non.na, non.na]
  }
  groups.snp.f <- factor(as.numeric(geno.snp))
  X <- stats::model.matrix(~., data = data.frame(genotype = groups.snp.f), contrasts.arg = list("genotype" = "contr.sum"))     
  p <- ncol(X) - 1  # nb.gp - 1
  n <- nrow(X)      # nb of non-NA individuals
  H <- tcrossprod(tcrossprod(X, solve(crossprod(X))), X)
  if (permute){
    perm <- sample(1:n)
    D <- D[perm, perm]
  }
  G <- gower(D)
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
  
  return(data.frame(pv.snp))
} 

##' Computes empirical P-value for a gene.
##' @title Compute empirical P-value
##' @param genotype.gene a data.frame of genotypes produced by \code{\link{read.bedix}}.
##' @param tre.dist a interdistance matrix produced by \code{\link{hellingerDist}}.
##' @param best.snp SNP with the smallest observed nominal P-value, computed by \code{\link{compute.nominal.pv}}.
##' @param min.pv.obs smallest observed nominal P-value.
##' @param min.nb.ext.scores the minimum number of permuted  nominal P-values lower than
##' the smallest observed nominal P-value to allow the computation to stop. Default is 100. 
##' @param nb.perm.max the maximum number of permutations. Default is 1000. 
##' @param verbose Default is TRUE. Mainly for debugging.
##' @return a data.frame with columns:
##' \item{variants.cis}{the number of variants tested in cis.}
##' \item{LD}{a linkage disequilibrium estimate for the genomic window (median r2).}
##' \item{best.snp}{ID of the SNP with the smallest observed nominal P-value.}
##' \item{best.nominal.pv}{P-value corresponding to the best SNP.}
##' \item{shape1}{Beta distribution parameter shape1.}
##' \item{shape2}{Beta distribution parameter shape2.}
##' \item{nb.perms}{the number of permutations used for the empirical P-value computation.}
##' \item{pv.emp}{empirical P-value based on permutations.}
##' \item{pv.emp.beta}{empirical P-value based on the beta approximation.}
##' \item{runtime}{approximated computation time per gene.}
##' @author Diego Garrido-Martín 
##' @keywords internal
compute.empirical.pv <- function(genotype.gene, tre.dist, best.snp, min.pv.obs, min.nb.ext.scores = 1e2, nb.perm.max = 1e3, comp.ld = TRUE, verbose = FALSE){
  if(comp.ld) {
    compute.ld <- function(df){
      M <- as.matrix(df)
      ns <- dim(M)[1]
      if(ns == 1) {
        return(NA)
      }
      R <- matrix(NA, ncol = ns, nrow = ns)               
      for (i in 1:(ns - 1)){
        for (j in (1 + i):ns){
          R[i, j] <- cor(M[i, ], M[j, ])^2
        }
      }
      return(median(R, na.rm = T)) 
    }
    ld <- compute.ld(genotype.gene[, labels(tre.dist)])
  }
  if (verbose) {
    message ("\tAdaptative permutation scheme")
  }
  t0 <- Sys.time()
  store.perm = c()
  i <- 1 
  ext <- 0
  pv <- 1
  while ( ext < min.nb.ext.scores & i <= nb.perm.max ) {
    # Note that i starts in 1. Thus "<=" instead of "<" in the while condition
    if (verbose & i%%100 == 0) message (sprintf("\t\tpermutation %s",i))
    min.pv.perm <- min(dplyr::do(dplyr::group_by(genotype.gene, snpId),compute.nominal.pv(., tre.dist, permute = TRUE))$pv.snp)
    store.perm <- c(store.perm, min.pv.perm)
    ext <- sum(store.perm <= min.pv.obs)
    pv <- (ext + 1)/(i + 1)
    i <- i + 1
  }
  fit <- fitdistrplus::mledist(store.perm, "beta")  # Estimate beta parameters on permutations 
  shape1 <- fit$estimate[1]                                               
  shape2 <- fit$estimate[2]
  names(shape1) <- names(shape2) <- NULL
  pv.beta <- pbeta(min.pv.obs, shape1 = shape1, shape2 = shape2)  
  t1 <- Sys.time()
  t.run <- as.numeric(difftime(t1, t0, units = "mins"))
  res.df <- data.frame(variants.cis = dim(genotype.gene)[1], LD = ld, best.snp = best.snp, best.nominal.pv = min.pv.obs, 
                       shape1 = shape1, shape2 = shape2, nb.perms = length(store.perm), pv.emp = pv, pv.emp.beta = pv.beta,
                       runtime = t.run) 
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

