##' Group SNPs according to their LD.
##' 
##' Clusters the SNP IDs of those variants that fulfill the following criteria:
##' \itemize{
##' \item{LD (r2) >= \code{th}.}
##' \item{Relative differences in sQTL F score <= \code{tol}}
##' \item{Relative differences in svQTL F score <= \code{tol.svqtl}}}
##' @title Linkage disequilibrium (LD) filter
##' @param genotype.gene a data.frame with genotype info produced by 'read.bedix'.
##' @param tre.mt a matrix of transcript relative expression (samples x transcripts).
##' @param th r2 threshold over which SNPs will be clustered. Default is 1.
##' @param tol maximum relative difference in sQTL F score allowed to group the SNPs. Default is 0.05.
##' @param tol.svqtl same as \code{tol} for svQTL F score. Default is 0.25. 
##' If it is NULL, the svQTL test is not performed.
##' @return A data.frame with genotype info identical to \code{genotype.gene} with 
##' an additional column named LD, containing the IDs of the SNPs that are 
##' in LD >= th and NA otherwise.
##' @author Diego Garrido-Martín
##' @keywords internal
LD.filter <- function(genotype.gene, tre.mt, th = 1, tol = 0.05, tol.svqtl = 0.25)
{
    if (th < 0 || th > 1 || !is.numeric(th)){
        stop ("'th' must be a numeric value in [0,1].")
    }
    if (!is.numeric(tol) || tol < 0 || tol > 1){
        stop ("'tol' should be a numeric value in [0,1].") 
    }
    if (!is.null(tol.svqtl)){
        if(!is.numeric(tol.svqtl) || tol.svqtl < 0 || tol.svqtl > 1){
            stop ("'tol.svqtl' should be either NULL or a numeric value in [0,1].") 
        }
    }
    ids <- genotype.gene$snpId
    g <- genotype.gene[, rownames(tre.mt)]
    colnames(g) <- rownames(g) <- NULL
    nG <- apply(g, 1, function(x)(length(table(x[x > -1])))) # Get nb of groups 
    M <- t(g)
    M[M == -1] <- NA
    res3 <- computeLD(M = M[, nG == 3], ids = ids[nG == 3], 
                      tre.mt = tre.mt, th = th, tol = tol, tol.svqtl = tol.svqtl)
    res2 <- computeLD(M = M[, nG == 2], ids = ids[nG == 2], 
                      tre.mt = tre.mt, th = th, tol = tol, tol.svqtl = tol.svqtl)
    if(is.null(res3) & !is.null(res2)){ # Some checks
        res <- res2
    } else if (!is.null(res3) & is.null(res2)){
        res <- res3 
    } else {
        res <- rbind(res3,res2)
    }
    if (is.null(res)){
        genotype.gene$LD <- rep(NA, nrow(genotype.gene))
        return(genotype.gene)
    } else {
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

computeLD <- function(M, ids, tre.mt, th = 1, tol = 0.05, tol.svqtl = 0.25)
{
    M <- as.matrix(M)
    if(ncol(M) > 1){
        s <- ncol(M)
        R <- cor(M, use = "pairwise.complete.obs")
        R <- R^2
        blocks <- lapply(as.data.frame(R),function(x) which(x >= th))
        names(blocks) <- 1:s 
        blocks <- F.filter(blocks = blocks, tre.mt = tre.mt, 
                           M = M, tol = tol, tol.svqtl = NULL)
        if(!is.null(tol.svqtl)){
            blocks <- F.filter(blocks = blocks, tre.mt = tre.mt, 
                               M = M, tol = tol, tol.svqtl = tol.svqtl)
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
        } else {
            return(NULL)
        }
    } else {
        return(NULL)
    }
}

prune_reorder <- function (edges, R)
{
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

F.filter <- function(blocks, tre.mt, M, tol = 0.05, tol.svqtl = 0.25, svQTL = FALSE)
{
    if(!is.null(tol.svqtl)){
        tol <- tol.svqtl
        svQTL <- TRUE
    }
    Fs <- apply(M, 2, function(x) F.calc(tre.mt = tre.mt, snp = x, svQTL = svQTL))
    names(Fs) <- 1:length(blocks)
    d <- list()
    for (e in names(blocks)){
        subs <- blocks[[e]]
        d[[e]] <- abs(Fs[as.character(subs)] - Fs[e])/Fs[e]
    }
    h <- list()
    for (k in names(blocks)){
        h[k] <- list(blocks[[k]] [ which(d[[k]] <= tol) ]) 
    }
    return(h)
}

F.calc <- function(tre.mt, snp, svQTL = FALSE)
{
    if (any(is.na(snp))) {
        non.na <- !is.na(snp)
        snp <- snp[non.na]
        tre.mt <- tre.mt[non.na, ]
    }
    snp.f <- factor(snp)
    if(svQTL){
        bd <- vegan::betadisper(dist(tre.mt), snp.f, type = "centroid")
        bd.perm <- permutest.betadisper(bd, control = permute::how(nperm = 2)) 
        return(bd.perm$F)
    } else {
        dfnum <- nlevels(snp.f) - 1
        dfden <- nrow(tre.mt) - dfnum - 1
        tre.mt <- scale(tre.mt, center = T, scale = F)
        G <- tcrossprod(tre.mt)
        X <- stats::model.matrix(~., data = data.frame(genotype = snp.f), 
                                 contrasts.arg = list("genotype" = "contr.sum")) # Note contrast type    
        H <- tcrossprod(tcrossprod(X, solve(crossprod(X))), X)
        numer <- crossprod(c(H), c(G))
        trG <- sum(diag(G))
        denom <- trG - numer
        f.snp <- as.numeric((numer*dfden)/(denom*dfnum))  
        return(f.snp)
    }
}


