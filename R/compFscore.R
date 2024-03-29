##' Compute the F score, max diff ratio difference and transcripts that change the most.
##' Additionally, it can compute a P-value to assess the significance of the F score.
##' @title F score computation
##' @param geno.df a data.frame of one row with the genotype information for each sample.
##' @param tre.mt a matrix with the transcript relative expression (samples x transcripts). 
##' @param svQTL should svQTL test be performed in addition to sQTL. Default is \code{FALSE}.
##' @param asympt should significance for the F score (sQTL test) be computed using 
##' the \code{\link[CompQuadForm]{davies}} method in the \code{CompQuadForm} package. 
##' Default is \code{TRUE}.
##' @param res is \code{tre.mt} the residual of the regression of additional covariates. Default is \code{FALSE}
##' @return A data.frame with columns:
##' \item{F}{the F score.}
##' \item{nb.groups}{the number of genotype groups.}
##' \item{md}{the maximum difference in splicing ratios between genotype groups.}
##' \item{tr.first, tr.second}{the two transcripts that change the most.}
##' \item{info}{comma separated list with the number of individuals per genotype group: -1,0,1,2.}
##' \item{pv}{if \code{asympt} is \code{TRUE} a P-value for the F score is computed.}
##' @author Diego Garrido-Martín, Jean Monlong
##' @keywords internal
compFscore <- function(geno.df, tre.mt, svQTL = FALSE, asympt = TRUE, res = FALSE)
{
    if(nrow(geno.df) > 1){
        stop(geno.df$snpId[1], " SNP is duplicated in the genotype file.")
    }
    geno.snp <- as.numeric(geno.df[, rownames(tre.mt)])
    names(geno.snp) <- rownames(tre.mt)
    info.snp <- c()
    tb.snp <- table(geno.snp)
    for (gt in c("-1", "0", "1", "2")){
        info.snp[gt] <- ifelse(is.na(tb.snp[gt]), 0, tb.snp[gt])
    }
    if (any(geno.snp == -1)) {
        non.na <- geno.snp > -1
        geno.snp <- geno.snp[non.na]
        tre.mt <- tre.mt[non.na, ]
    }
    info.snp <- paste(info.snp, collapse =",")
    groups.snp.f <- factor(as.numeric(geno.snp))
    if(res) {
        tre.mt2md <- tre.mt
    } else{ # w/o covariates, undo sqrt transf and keep original MD
        tre.mt2md <- tre.mt^2 
    }   
    mdt <- md.trans(tre.mt2md, groups.snp.f)
    n <- nrow(tre.mt)
    nb.gp <- nlevels(groups.snp.f)
    dfnum <- nb.gp - 1
    dfden <- n - dfnum - 1
    tre.mt <- scale(tre.mt, center = TRUE, scale = FALSE)
    G <- tcrossprod(tre.mt)
    X <- stats::model.matrix(~., data = data.frame(genotype = groups.snp.f),
                             contrasts.arg = list("genotype" = "contr.sum")) # Note contrast type    
    H <- tcrossprod(tcrossprod(X, solve(crossprod(X))), X)
    numer <- crossprod(c(H), c(G))
    trG <- sum(diag(G))
    denom <- trG - numer
    f.tilde <- as.numeric(numer/denom)  
    if (asympt) {
        fit <- stats::lm(tre.mt ~ groups.snp.f)
        R <- fit$residuals
        e <- eigen(stats::cov(R)*(n-1)/dfden, symmetric = TRUE, only.values = TRUE)$values
        lambda <- abs(e[abs(e) > 1e-12])
        item.acc <- 1e-14
        pv.snp <- pcqf(q = f.tilde, lambda = lambda, df.i = dfnum, 
                       df.e = dfden, acc = item.acc)
        while (length(pv.snp) > 1) {
            item.acc <- item.acc * 10
            pv.snp <- pcqf(q = f.tilde, lambda = lambda, df.i = dfnum, 
                           df.e = dfden, acc = item.acc)
        }
        if (pv.snp < item.acc) pv.snp <- item.acc
        res.df <- data.frame(F = f.tilde*dfden/dfnum, nb.groups = nb.gp,
                             md = mdt$md, tr.first = mdt$tr.first, tr.second = mdt$tr.second,
                             info = info.snp, pv = pv.snp, stringsAsFactors = FALSE)
    } else {
        res.df <- data.frame(F = f.tilde*dfden/dfnum, nb.groups = nb.gp, 
                         md = mdt$md, tr.first = mdt$tr.first, tr.second = mdt$tr.second, 
                         info = info.snp, stringsAsFactors = FALSE)
    }
    if (svQTL) {
        bd <- vegan::betadisper(stats::dist(tre.mt), groups.snp.f, type = "centroid")
        bd.perm <- permutest.betadisper(bd, control = permute::how(nperm = 2)) 
        res.df$F.svQTL <- bd.perm$F
    }
    if (any(colnames(geno.df) == "LD")) {
        res.df$LD <- geno.df$LD
    }
    return(res.df)
}

pcqf <- function (q, lambda, df.i, df.e, lim = 50000, acc = 1e-14)
{
    gamma <- c(lambda, -q * lambda)
    nu <- c(rep(df.i, length(lambda)), rep(df.e, length(lambda)))
    pv <- suppressWarnings(CompQuadForm::davies(0, lambda = gamma, h = nu, lim = lim, acc = acc))
    if (pv$ifault != 0) {
        return(pv)
    }
    if (pv$Qq < 0 || pv$Qq > 1) {
        return(pv)
    }
    if (pv$ifault == 0) {
        return(pv$Qq)
    }
}