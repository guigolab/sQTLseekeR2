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
##' @param covariates a data.frame with covariate information per sample (samples x covariates).
##' Rownames should be the sample ids. Covariates can be either \code{numeric} or \code{factor}. 
##' When provided, they are regressed out before testing the genotype effect. Default is \code{NULL}.
##' @param genic.window the window(bp) around the gene in which the SNPs are tested. Default is 5000 (i.e. 5kb).
##' Same as in \code{sqtl.seeker}.
##' @param min.nb.ext.scores the minimum number of permuted  nominal P-values lower than
##' the lowest observed nominal P-value to allow the computation to stop. Default is 100. 
##' @param nb.perm.min the minimum number of permutations. Default is 100.
##' @param nb.perm.max the maximum number of permutations. Default is 1000. 
##' @param min.nb.ind.geno SNPs with less samples than \code{min.nb.ind.geno} in any genotype group
##' are filtered out. Default is 10.
##' @param verbose Default is TRUE. 
##' @return A data.frame with columns:
##' \item{geneId}{the gene name.}
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
sqtl.seeker.p <- function(tre.df, genotype.f, gene.loc, covariates = NULL, genic.window = 5000, nb.perm.min = 100, nb.perm.max = 1000, min.nb.ext.scores = 100, min.nb.ind.geno = 10, verbose = FALSE){
  
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
      if(!is.null(covariates)){
        if(!all(rownames(covariates) %in% colnames(tre.gene))){
          stop("All samples should have covariate information (either a value or NA).")                  
        }
        cov.na <- apply(covariates, 1, function(x){any(is.na(x))})
        covariates <- covariates[!cov.na, ]                                     
        if(sum(cov.na) > 0){
          warning(sprintf("%s samples with NA values for at least one covariate have been removed.", sum(cov.na)))     
        }
        if(sum(cov.na) > round(0.05 * nrow(covariates))){
          stop("More than 5% of the samples contain NA values for at least one covariate")  
        }
      }
      ## Remove samples with non expressed genes
      tre.gene <- tre.gene[, !is.na(tre.gene[1, ])]
      ## Focus on common samples
      genotype.headers <- as.character(utils::read.table(genotype.f, as.is = TRUE, nrows = 1))
      if(!is.null(covariates)){
        com.samples <- Reduce(intersect, list(colnames(tre.gene), genotype.headers, rownames(covariates)))   
        if(length(com.samples) == 0){
          stop("No common samples between genotype, covariate and transcript files.")
        }
      } else{
        com.samples <- intersect(colnames(tre.gene), genotype.headers)
        if(length(com.samples) == 0){
          stop("No common samples between genotype and transcript files.")
        }
      }
      tre.gene <- tre.gene[, c("trId", "geneId", com.samples)]                     
      tre.tc <- t(sqrt(tre.gene[,com.samples]))
      tre.tc <- scale(tre.tc, center = T, scale = F)
      colnames(tre.tc) <- tre.gene$tr
      
      # Here regress covariates and keep residual if applies
      if(!is.null(covariates)){
        fit <- lm(tre.tc ~ ., data = covariates)
        tre.tc <- fit$residual
      }
      
      res.df <- data.frame()
      genotype.gene <- read.bedix(genotype.f, gr.gene) 
      if (verbose & is.null(genotype.gene)) {
        message("\tNo SNPs in the genomic range.")
      }
      if (!is.null(genotype.gene)) {
        snps.to.keep <- check.genotype(genotype.gene[, com.samples], tre.gene[, com.samples], min.nb.ind.geno = min.nb.ind.geno)
        if (verbose) {
          snps.to.keep.t <- table(snps.to.keep)
          message("\t", paste(names(snps.to.keep.t), snps.to.keep.t, sep = ": ", collapse = ", "))
        }
        if (any(snps.to.keep == "PASS")) {
          genotype.gene <- genotype.gene[snps.to.keep == "PASS", ]
          observed <- dplyr::do(dplyr::group_by(genotype.gene, snpId), compute.nominal.pv(., tre.mt = tre.tc))
          min.pv.obs <- min(observed$pv.snp) 
          best.snp <- observed$snpId[which.min(observed$pv.snp)]
          res.df <- compute.empirical.pv(genotype.gene = genotype.gene, tre.mt = tre.tc, min.nb.ext.scores = min.nb.ext.scores, nb.perm.min = nb.perm.min,
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

##' Computes nominal P-value for the association between the genotype of a locus (e.g. SNP) and the splicing ratios of a gene.
##' @title Compute nominal P-value
##' @param geno.df a data.frame of one row with the genotype information for each sample.
##' @param tre.mt a matrix of splicing ratios (samples x transcripts).
##' @param permute should the rows of the splicing ratio matrix be permuted. Default is FALSE.
##' @param seed if \code{permute} is TRUE, value provided to \code{\link{set.seed}}. Default is 1.
##' @return A data.frame containing a P-value for the association.
##' @author Diego Garrido-Martín 
##' @keywords internal
compute.nominal.pv <- function(geno.df, tre.mt, permute = FALSE, seed = 1){
  
  if (nrow(geno.df) > 1) {
    stop(geno.df$snpId[1], " SNP is duplicated in the genotype file.")
  } 
  geno.snp <- as.numeric(geno.df[, rownames(tre.mt)])
  names(geno.snp) <- rownames(tre.mt)
  if (any(geno.snp == -1)) {
    non.na <- geno.snp > -1
    geno.snp <- geno.snp[non.na]
    tre.mt <- tre.mt[non.na, ]
  }
  groups.snp.f <- factor(as.numeric(geno.snp))
  n <- nrow(tre.mt)

  if (permute){
    set.seed(seed)
    perm <- sample(1:n)
    tre.mt <- tre.mt[perm, ]
  }
  
  fit <- lm(tre.mt ~ groups.snp.f)
  G <- tcrossprod(tre.mt)
  X <- stats::model.matrix(fit) # Note contrast type
  H <- tcrossprod(tcrossprod(X, solve(crossprod(X))), X)
  numer <- crossprod(c(H), c(G))
  trG <- sum(diag(G))
  denom <- trG - numer
  f.tilde <- as.numeric(numer/denom)  
  R <- fit$residuals
  df.e <- fit$df.residual
  df.i <- nlevels(groups.snp.f) - 1
  e <- eigen(cov(R)*(n-1)/df.e, symmetric = T, only.values = T)$values
  lambda <- abs(e[abs(e) > 1e-12])
  
  item.acc <- 1e-14
  pv.snp <- pcqf(q = f.tilde, lambda = lambda, df.i = df.i, df.e = df.e, acc = item.acc)
  while (length(pv.snp) > 1) {
    item.acc <- item.acc * 10
    pv.snp <- pcqf(q = f.tilde, lambda = lambda, df.i = df.i, df.e = df.e, acc = item.acc)
  }
  if (pv.snp < item.acc) {
    pv.snp <- item.acc
  }
  return(data.frame(pv.snp))
} 


##' Computes empirical P-value for a gene.
##' @title Compute empirical P-value
##' @param genotype.gene a data.frame of genotypes produced by \code{\link{read.bedix}}.
##' @param tre.mt a matrix of splicing ratios (samples x transcripts).
##' @param best.snp SNP with the smallest observed nominal P-value, computed by \code{\link{compute.nominal.pv}}.
##' @param min.pv.obs smallest observed nominal P-value.
##' @param nb.perm.min the minimum number of permutations. Default is 100.
##' @param nb.perm.max the maximum number of permutations. Default is 1000. 
##' @param min.nb.ext.scores the minimum number of permuted  nominal P-values lower than
##' the smallest observed nominal P-value to allow the computation to stop. Default is 100. 
##' @param comp.ld should linkage disequilibrium estimates be computed (median r2). Default is TRUE.
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
compute.empirical.pv <- function(genotype.gene, tre.mt, best.snp, min.pv.obs, nb.perm.min = 100, nb.perm.max = 1000, min.nb.ext.scores = 100, comp.ld = TRUE, verbose = FALSE){
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
    ld <- compute.ld(genotype.gene[, rownames(tre.mt)])
  }
  variants.cis <- dim(genotype.gene)[1]
  genotype.gene <- LD.filter(genotype.gene = genotype.gene, tre.mt = tre.mt, th = 1, tol = 0, tol.svqtl = NULL)             
  genotype.gene$LD <- NULL
  if (verbose) {
    message(sprintf("\tPASS not in perfect LD: %s", nrow(genotype.gene)))
    message ("\tAdaptative permutation scheme")
  }
  t0 <- Sys.time()
  store.perm = c()
  i <- 1 
  ext <- 0
  pv <- 1
  while ( i <= nb.perm.min || (ext < min.nb.ext.scores && i <= nb.perm.max) ) {
    # Note that i starts in 1. Thus "<=" instead of "<" in the while condition
    if (verbose & i%%100 == 0) message (sprintf("\t\tpermutation %s",i))
    min.pv.perm <- min(dplyr::do(dplyr::group_by(genotype.gene, snpId),compute.nominal.pv(., tre.mt, permute = TRUE, seed = i))$pv.snp)
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
  res.df <- data.frame(variants.cis = variants.cis, LD = ld, best.snp = best.snp, best.nominal.pv = min.pv.obs, 
                       shape1 = shape1, shape2 = shape2, nb.perm = length(store.perm), pv.emp.perm = pv, pv.emp.beta = pv.beta,
                       runtime = t.run) 
  return(res.df)
} 
