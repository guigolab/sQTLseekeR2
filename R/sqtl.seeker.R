##' \code{sqtl.seeker} is the main function of \code{sQTLseekeR2} package. From transcript relative expression,
##' prepared using \code{prepare.trans.exp}, information about the gene location and the path to an ordered genotype
##' file, indexed by the function \code{index.genotype}, association between each SNP and the transcript relative
##' expression is tested. Eventually, svQTL, i.e. SNPs affecting splicing variability can also be tested to pinpoint
##' potentially false sQTL (see Details).
##'
##' A set of filters is automatically set to remove SNPs which are unpractical or not informative (see also \code{\link{check.genotype}} function). Precisely, these
##' filters remove SNP with:
##' \itemize{
##' \item{more than 3 unknown genotypes}
##' \item{less than min.nb.ind.geno samples in any genotype group}
##' \item{less than 5 different splicing pattern (needed for permutation efficiency) in any genotype group}}
##'
##' Testing differences in transcript relative expression between genotype groups assumes homogeneity of the variances
##' in these groups. Testing this assumption is more complex and computationnally intensive but if needed the user
##' can choose to test for svQTL (splicing variability QTL), i.e. gene/SNPs where this assumption is violated,
##' by using the \code{svQTL = TRUE}. This test is run in parallel to the sQTL tests, but the computation time
##' will be considerably higher. For this reason, another parameter can be tweaked, \code{nb.perm.max.svQTL},
##' to reduce the number of permutations for the svQTL tests if needed for feasibility reasons.
##'
##' The permutation process is optimized by computing one permuted distribution per gene and using a number
##' of permutations dependant on how extreme are the observed scores when compared to the permuted ones. To decrease
##' even more the computation time, an approximation of the null F distribution was given by Anderson &
##' Robinson (2003), as a misture of Chi-square distributions whose parameters are derived from the eigen
##' values of the distance matrix. 
##'
##' In addition to the F score and P-value, the maximum difference(MD) in relative expression between genotype
##' groups is reported. This is to be used as a measure of the size of the effect. For example, if 'md' is 0.2
##' there is one transcript whose relative expression shifted by 20 percent between two genotype groups.
##' @title sQTL seeker
##' @param tre.df a data.frame with transcript relative expression
##' produced by 'prepare.trans.exp'.
##' @param genotype.f the name of the genotype file. This file needs to
##' be ordered by position, compressed and indexed using 'index.genotype' or externally using tabix (samtools).
##' Must have column 'snpId'.
##' @param gene.loc a data.frame with the genes location. Columns 'chr', 'start',
##' 'end' and 'geneId' are required.
##' @param covariates a data.frame with covariate information per sample (samples x covariates).
##' Rownames should be the sample ids. Covariates can be either \code{numeric} or \code{factor}. 
##' When provided, they are regressed out before testing the genotype effect. Default is \code{NULL}.
##' @param genic.window the window(bp) around the gene in which the SNPs are tested. Default is 5000 (i.e. 5kb).
##' @param min.nb.ext.scores the minimum number of permuted score higher than
##' the highest true score to allow the computation to stop. Default is 1000.
##' @param nb.perm.max the maximum number of permutations. Default is 1e6.
##' @param nb.perm.max.svQTL the maximum number of permutations for the svQTL computation. Default is 1e4.
##' @param svQTL should svQTLs test be performed in addition to sQTLs. Default is \code{FALSE}. Warning:
##' computation of svQTLs cannot rely on asymptotic approximation, hence the heavy permutations will
##' considerably increase the running time.
##' @param asympt should the asymptotic null distribution be used to assess significance instead of permutations. 
##' The \code{\link[CompQuadForm]{davies}} method in the \code{CompQuadForm} package is employed to compute P-values.
##' Default is \code{TRUE}.
##' @param ld.filter Linkage disequilibrium threshold (r2) over which variants should be merged,
##' so that only one SNP per LD block is tested. Only variants over the treshold that have highly similar 
##' pseudo F scores will be merged. Default is \code{NULL}.  
##' @param min.nb.ind.geno SNPs with less samples than \code{min.nb.ind.geno} in any genotype group
##' are filtered out. Default is 10.
##' @param verbose Should the gene IDs be outputed when analyzed. Default is \code{TRUE}. Mainly for debugging.
##' @return a data.frame with columns
##' \item{geneId}{the gene name.}
##' \item{snpId}{the SNP name.}
##' \item{F}{the F score.}
##' \item{nb.groups}{the number of genotype groups (2 or 3).}
##' \item{md}{the maximum difference in relative expression between genotype groups (see Details).}
##' \item{tr.first/tr.second}{the transcript IDs of the two transcripts that change the most (and symetrically).}
##' \item{info}{comma separated list with the individuals per genotype group: -1,0,1,2.}
##' \item{pv}{the P-value.}
##' \item{nb.perms}{the number of permutations used for the P-value computation.}
##' \item{F.svQTL/pv.svQTL/nb.perms.svQTL}{idem for svQTLs, if \code{svQTL} is \code{TRUE}'.}
##' \item{LD}{if ld.filter is not NULL, the variants in high LD (r2 >= ld.filter) with the tested variant that also have a similar Fscore.}
##' @author Jean Monlong, Diego Garrido-Martín
##' @export
sqtl.seeker <- function(tre.df, genotype.f, gene.loc, covariates = NULL, 
                        genic.window = 5000, min.nb.ext.scores = 1000, 
                        nb.perm.max = 1e6, nb.perm.max.svQTL = 1e4, 
                        svQTL = FALSE, asympt = TRUE, ld.filter = NULL, 
                        min.nb.ind.geno = 10, verbose = FALSE)
{
    . <- nb.groups <- snpId <- NULL ## Uglily appease R checks (dplyr)
    analyze.gene.f <- function(tre.gene){
        if(verbose) message(tre.gene$geneId[1])
        if(sum(duplicated(gene.loc$geneId)) > 1){
            stop(tre.gene$geneId[1], " Repeated gene in gene location file.")
        }
        gr.gene <- with(gene.loc[which(gene.loc$geneId == tre.gene$geneId[1]), ],
                        GenomicRanges::GRanges(chr, IRanges::IRanges(start, end)))
        if(genic.window > 0){
            gr.gene <- GenomicRanges::resize(gr.gene, GenomicRanges::width(gr.gene) +
                                               2 * genic.window, fix = "center")
        }
        if(length(gr.gene) > 0){
            if(!is.null(covariates)){
                if(!all(rownames(covariates) %in% colnames(tre.gene))){
                    stop("All samples should have covariate information (either a value or NA).")                  
                }
                cov.na <- apply(covariates, 1, function(x){any(is.na(x))})
                covariates <- covariates[!cov.na,  ,drop = FALSE]                                     
                if(sum(cov.na) > 0){
                    warning(sprintf("%s samples with NA values for at least one covariate have been removed.", 
                                    sum(cov.na)))     
                }
                if(sum(cov.na) / nrow(covariates) > 0.05){
                    stop("More than 5% of the samples contain NA values for at least one covariate")  
                }
            }
            tre.gene <- tre.gene[, !is.na(tre.gene[1, ])]
            genotype.headers <- as.character(utils::read.table(genotype.f,
                                                               as.is = TRUE, nrows = 1))
            if(!is.null(covariates)){
                com.samples <- Reduce(intersect, list(colnames(tre.gene),
                                                      genotype.headers, rownames(covariates)))   
                if (length(com.samples) == 0) {
                    stop("No common samples between genotype, covariate and transcript files.")
                }
            } else{
                com.samples <- intersect(colnames(tre.gene), genotype.headers)
                if (length(com.samples) == 0) {
                    stop("No common samples between genotype and transcript files.")
                }
            }
            tre.gene <- tre.gene[, c("trId", "geneId", com.samples)]
            allzero <- apply(tre.gene[, com.samples], 1, function(x){ sum(x) == 0 })
            if (any(allzero)){
                tr.out <- tre.gene$trId[allzero]
                tre.gene <- tre.gene[!allzero,]
                warning(sprintf("Transcript(s) %s is(are) removed due to zero expression in all common samples.", 
                                paste(tr.out, collapse = ", ")) )
            }
            tre.tc <- t(sqrt(tre.gene[, com.samples]))
            colnames(tre.tc) <- tre.gene$tr
            if(!is.null(covariates)){
                covariates <- covariates[com.samples, , drop = FALSE]
                multiclass <- apply(covariates, 2, function(x){length(table(x)) > 1})
                covariates <- covariates[, multiclass, drop = FALSE]
                if (verbose & any(!multiclass)){
                  message("\t", "Covariates removed due to only one value: ", 
                          paste(names(multiclass)[!multiclass], collapse = ", "))
                }
                fit <- stats::lm(tre.tc ~ ., data = covariates)
                if (ncol(covariates) > 1){
                    vifs <- car::vif(stats::lm(tre.tc[, 1] ~ ., data = covariates))
                    if (verbose){
                        message("\t", "Covariates VIF - ", 
                                paste(names(vifs), round(vifs, 2), sep = ": ", collapse = ", "))
                    }
                    if (any(vifs > 5)){
                        warning("Check multicollinearity. VIF > 5 for some covariates:", "\n", 
                                paste(names(vifs), round(vifs, 2), sep = ": ", collapse = ", "))
                    }  
                }
                tre.tc <- fit$residual
            }
            res.df <- data.frame()
            if(GenomicRanges::width(gr.gene) > 20000 && is.null(ld.filter)){
                pos.breaks <- unique(round(seq(GenomicRanges::start(gr.gene), 
                                               GenomicRanges::end(gr.gene), 
                                               length.out = floor(GenomicRanges::width(gr.gene)/10000) + 1)))
                gr.gene.spl <- rep(gr.gene, length(pos.breaks) - 1)
                GenomicRanges::start(gr.gene.spl) <- pos.breaks[-length(pos.breaks)]
                pos.breaks[length(pos.breaks)] <- pos.breaks[length(pos.breaks)] + 1
                GenomicRanges::end(gr.gene.spl) <- pos.breaks[-1] - 1
            } else {
                gr.gene.spl <- gr.gene
            }
            res.df <- lapply(1:length(gr.gene.spl), function(ii){
                res.range <- data.frame()
                if(verbose) message("  Sub-range ", ii)
                genotype.gene <- read.bedix(genotype.f, gr.gene.spl[ii])
                if(verbose && is.null(genotype.gene)) message("\tNo SNPs in the genomic range.") 
                if(!is.null(genotype.gene)){
                    snps.to.keep <- check.genotype(genotype.gene[, com.samples], 
                                                   tre.gene[, com.samples], 
                                                   min.nb.ind.geno = min.nb.ind.geno)
                    if(verbose){
                        snps.to.keep.t <- table(snps.to.keep)
                        message("\t", paste(names(snps.to.keep.t), snps.to.keep.t, 
                                            sep = ": ", collapse=", "))
                    }
                    if(any(snps.to.keep == "PASS")){
                        genotype.gene <- genotype.gene[snps.to.keep == "PASS", ]
                        if(!is.null(ld.filter)){
                            if(verbose) message("\tLD filtering")
                            genotype.gene <- LD.filter(genotype.gene = genotype.gene,
                                                       tre.mt = tre.tc, th = ld.filter,
                                                       svQTL = svQTL)
                        }
                        res.range <- dplyr::do(dplyr::group_by(genotype.gene, snpId),
                                               compFscore(., tre.tc, svQTL = svQTL, 
                                                          asympt = asympt, 
                                                          res = !is.null(covariates)))
                    }
                }
                return(res.range)
            })
            range.done <- which(unlist(lapply(res.df, nrow)) > 0)
            if(length(range.done) > 0){
                res.df <- res.df[range.done]
                res.df <- do.call(rbind, res.df)
                if (!is.null(ld.filter)){
                    ld <- res.df[,c("snpId","LD")] 
                    res.df$LD <- NULL
                }
                if(!asympt){
                    res.df <- dplyr::do(dplyr::group_by(res.df, nb.groups),
                                        compPvalue(., tre.tc, min.nb.ext.scores = min.nb.ext.scores,
                                                   nb.perm.max = nb.perm.max)) 
                }
                if(svQTL){
                    res.df <- dplyr::do(dplyr::group_by(res.df, nb.groups),
                                        compPvalue(., tre.tc, svQTL = TRUE,
                                                   min.nb.ext.scores = min.nb.ext.scores,
                                                   nb.perm.max = nb.perm.max.svQTL))
                }
                if (!is.null(ld.filter)){
                    res.df <- merge(res.df, ld, by = "snpId")
                }
                res.df <- dplyr::arrange(res.df, pv)
                return(data.frame(done = TRUE, res.df))
            }
        } else {
            if(verbose) warning("Issue with the gene location.")
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

##' Labels SNPs according to their suitability for sQTL mapping.  
##' @title Check genotype
##' @param geno.df a data.frame of genotypes produced by \code{\link{read.bedix}}.
##' @param tre.df a data.frame with transcript relative expression values
##' produced by \code{\link{prepare.trans.exp}} corresponding to one gene.
##' @return A character vector of SNP suitabilities. Possible labels:
##' \item{"Missing genotype" (more than 3 missing genotype values)}{ }
##' \item{"Only one genotype group"}{ }
##' \item{"Not all the groups with >10 samples"}{ }
##' \item{"Not all the groups with >5 different splicing patterns"}{ }
##' \item{"PASS"}{ }
##' @author Jean Monlong, Diego Garrido-Martín 
##' @keywords internal
check.genotype <- function(geno.df, tre.df, min.nb.ind.geno = 10)
{
    apply(geno.df, 1, function(geno.snp){
        if(sum(as.numeric(geno.snp) == -1) > 2){
            return("Missing genotype")
        }
        geno.snp.t <- table(geno.snp[geno.snp > -1])
        if (length(geno.snp.t) < 2) {
            return("Only one genotype group")                    
        }
        if (sum(geno.snp.t >= min.nb.ind.geno) < length(geno.snp.t)) {
            return(sprintf("Not all the groups with >%s samples", min.nb.ind.geno))        
        }
        nb.diff.pts <- sapply(names(geno.snp.t)[geno.snp.t > 1], function(geno.i){
            nbDiffPt(tre.df[, which(geno.snp == geno.i)])
        })
        if(sum(nb.diff.pts >= 5) < length(geno.snp.t)){
            return("Not all the groups with >5 different splicing patterns")
        }
        return("PASS")
    })
}
