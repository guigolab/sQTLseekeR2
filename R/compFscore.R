##' Compute the F score, max diff ratio difference and transcripts that change the most.
##' Additionally, it can compute a P-value to assess the significance of the F score.
##' @title F score computation
##' @param geno.df a data.frame of one row with the genotype information for each sample.
##' @param tre.mt a matrix with the transcript relative expression (samples x transcripts). 
##' @param svQTL should svQTL test be performed in addition to sQTL. Default is \code{FALSE}.
##' @param asympt should significance for the F score (sQTL test) be computed using 
##' the \code{\link[CompQuadForm]{farebrother}} method in the \code{CompQuadForm} package. 
##' Default is \code{TRUE}.
##' @param res is \code{tre.mt} the residual of the regression of additional covariates. Default is \code{FALSE}
##' @param condition covariate (factor with two levels) to be tested together with the genotype. The interaction term will be also assessed. 
##' Default is \code{NULL}. This will enable the asymptotic mode. Ignored if \code{res} is \code{FALSE}. 
##' @return A data.frame with columns:
##' \item{F}{the F score.}
##' \item{nb.groups}{the number of genotype groups.}
##' \item{md}{the maximum difference in splicing ratios between genotype groups.}
##' \item{tr.first, tr.second}{the two transcripts that change the most.}
##' \item{info}{comma separated list with the number of individuals per genotype group: -1,0,1,2.}
##' \item{pv}{if \code{asympt} is \code{TRUE} a P-value for the F score is computed.}
##' @author Diego Garrido-Martín, Jean Monlong
##' @keywords internal
compFscore <- function(geno.df, tre.mt, svQTL = FALSE, asympt = TRUE, res = FALSE,
                       condition = NULL)
{
    if(nrow(geno.df) > 1){
        stop(geno.df$snpId[1], " SNP is duplicated in the genotype file.")
    }
    if(!is.null(condition) && !res){
        condition <- NULL
        warning("'condition' will be set to NULL (res is FALSE)") # Double-check
    }
    geno.snp <- as.numeric(geno.df[, rownames(tre.mt)])
    names(geno.snp) <- rownames(tre.mt)
    info.snp <- c()
    tb.snp <- table(geno.snp)
    for (gt in c("-1", "0", "1", "2")){
        info.snp[gt] <- ifelse(is.na(tb.snp[gt]), 0, tb.snp[gt])
    }
    if (!is.null(condition)){
      condition <- condition[, 1]
      names(condition) <- rownames(tre.mt)
      info.condition <- paste(as.character(table(condition)), collapse = ",")
    }
    if (any(geno.snp == -1)) {
        non.na <- geno.snp > -1
        geno.snp <- geno.snp[non.na]
        tre.mt <- tre.mt[non.na, ]
        if (!is.null(condition)){
          condition <- condition[non.na]
        }
    }
    info.snp <- paste(info.snp, collapse =",")
    groups.snp.f <- factor(as.numeric(geno.snp))
    if(res) {
        tre.mt2md <- tre.mt
    } else{ # w/o covariates, undo sqrt transf and keep original MD
        tre.mt2md <- tre.mt^2 
    }   
    mdt <- md.trans(tre.mt2md, groups.snp.f)
 
    if (is.null(condition)){
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
        pv.snp <- pcqf(q = numer, lambda = lambda, df.i = dfnum, acc = item.acc)
        while (length(pv.snp) > 1) {
          item.acc <- item.acc * 10
          if(item.acc < 1e-3){
            stop("Accuracy requested below 1e-3")
          }
          pv.snp <- pcqf(q = numer, lambda = lambda, df.i = dfnum, acc = item.acc)
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
    } else {
        info.interaction <- paste(as.character(table(groups.snp.f:condition)), collapse = ",")
        ls <- names(table(condition))
        mdt_l1 <- md.trans(tre.mt2md[condition == ls[1], ],
                           groups.snp.f[condition == ls[1]])
        mdt_l2 <- md.trans(tre.mt2md[condition == ls[2], ],
                           groups.snp.f[condition == ls[2]])
        tre.mt <- scale(tre.mt, center = TRUE, scale = FALSE)
        fit <- stats::lm(tre.mt ~ groups.snp.f + condition + groups.snp.f:condition,
                  contrasts = list(groups.snp.f = "contr.sum", condition = "contr.sum"))
        R <- fit$residuals
        n <- nrow(R)
        UU <- car::Anova(fit, type = "II") 
        SS <- lapply(UU$SSP, function(x){sum(diag(x))})
        SSe <- sum(diag(UU$SSPE))
        f.tilde <- unlist(lapply(SS, function(x){x/SSe}))
        Df <- UU$df
        df.e <- fit$df.residual
        e <- eigen(stats::cov(R)*(n-1)/df.e, symmetric = T, only.values = T)$values
        e <- abs(e[abs(e) > 1e-12])
        item.acc <- 1e-14
        pvs <- mapply(pv.ss, ss = unlist(SS), df.i = Df, MoreArgs = list(lambda = e))
        Fs <- f.tilde*df.e/Df
        res.df <- data.frame(F_snp = Fs[1], F_cond = Fs[2], F_snpXcond = Fs[3],
                             pv_snp = pvs[1], pv_cond = pvs[2], pv_snpXcond = pvs[3], 
                             nb.groups = nlevels(groups.snp.f), md_snp = mdt$md, 
                             tr.first_snp = mdt$tr.first, tr.second_snp = mdt$tr.second,
                             md_cond_l1 = mdt_l1[[1]], md_cond_l2 = mdt_l2[[1]], 
                             tr.agree = sum(unlist(mdt_l1[2:3])%in%unlist(mdt_l2[2:3])),                                                                   
                             info_snp = info.snp, info_cond = info.condition, 
                             info_snpXcond = info.interaction, stringsAsFactors = FALSE) 
        if (svQTL) {
          bd <- vegan::betadisper(stats::dist(tre.mt), groups.snp.f:cond, type = "centroid")
          bd.perm <- permutest.betadisper(bd, control = permute::how(nperm = 2)) 
          res.df$F.svQTL <- bd.perm$F
        }
    }

    if (any(colnames(geno.df) == "LD")) {
        res.df$LD <- geno.df$LD
    }
    return(res.df)
}

pv.ss <- function(ss, lambda, df.i, acc = 1e-14, max.acc = 1e-8 , tol = 1e-3)
{
  start.acc <- acc
  lambda <- lambda[lambda/sum(lambda) > tol]
  pv <- pcqf(q = ss, lambda = lambda, df.i = df.i, acc = acc)
  while (length(pv) > 1) {
    acc <- acc * 10
    if(acc > max.acc){
      stop("Accuracy requested in CompQuadForm::farebrother above 1e-8.")
    }
    pv  <- pcqf(q = ss, lambda = lambda, df.i = df.i, acc = acc)
  }
  if (pv < start.acc) {
    pv <- start.acc
  }
  return(pv)
}

pcqf <- function (q, lambda, df.i, acc = 1e-14, min.time = 10)
{
    t0 <- Sys.time()
    pv <- CompQuadForm::farebrother(q, lambda = lambda, h = rep(df.i, length(lambda)), eps = acc) # Add suppressWarnings ?
    t1 <- Sys.time()
    time <- as.numeric(difftime(t1, t0, units = "secs"))
    if(time > min.time) {warning(sprintf("CompQuadForm::farebrother runs too slow: %.2f secs.", time))}
    
    if (pv$Qq < 0 || pv$Qq > 1) { # This is equivalent to ifaults 5 and 9. Can be solved reducing acc
        warning(sprintf("CompQuadForm::farebrother Qq is %s.", pv$Qq))
        return(pv)
    }
    if (! pv$ifault %in% c(0,4)) {
      stop(sprintf("CompQuadForm::farebrother ifault is %s.", pv$ifault))
    }
    if (pv$ifault %in% c(0, 4)) {
        return(pv$Qq)
    }
}
