
# Big input files
genotype.indexed.f <- "~/PhD/projects/Blueprint/run/input/snps.tsv.bgz"
gene.bed.f <- "~/PhD/projects/Blueprint/run/input/genes.bed"
gene.bed <- read.table(gene.bed.f, as.is=TRUE, sep="\t")
colnames(gene.bed) <- c("chr","start","end","geneId")
load(file="~/PhD/projects/sqtlseeker/prepare.trans.exp/betadisperVSme.RData")
colnames(tre.df) <- gsub(pattern="(S.{5}).*", replacement="\\1", perl=TRUE,colnames(tre.df)) 
gene.bed <- subset(gene.bed, geneId %in% tre.df$geneId) # Remove from genes.bed all genes that are not in tre.df


## Define parameters 
library(sQTLseekeR2)
genic.window <- 5000
gene.loc <- gene.bed
set.seed(123)
W <- 1

several <- function (edges){
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
  dd <- data.frame(nedges,redges)
  ordered <- rownames(dd[with(dd,order(-nedges,-redges)),])
  nedges <- nedges[ordered]
  return(list(edges,nedges))
}
f.calc <- function(tre.dist,snp){
  G <- sQTLseekeR2:::gower(tre.dist)
  X <- stats::model.matrix(~., data = data.frame(genotype = factor(snp)), contrasts.arg = list("genotype" = "contr.sum"))     
  H <- tcrossprod(tcrossprod(X, solve(crossprod(X))), X)
  numer <- crossprod(c(H), c(G))
  trG <- sum(diag(G))
  denom <- trG - numer                                                          
  f.snp <- as.numeric(numer/denom)  
}

Rstore <-c()
Nstore <-c()
for (gene.i in unique(tre.df$geneId)){
  cat(W,"\n")
  #gene.i <- unique(tre.df$geneId)[12] # For debugging
  tre.gene <- subset(tre.df, geneId == gene.i)
  message(tre.gene$geneId[1])
  
  gr.gene <- with(gene.loc[which(gene.loc$geneId == tre.gene$geneId[1]), ], GenomicRanges::GRanges(chr, IRanges::IRanges(start, end)))
  gr.gene <- GenomicRanges::resize(gr.gene, GenomicRanges::width(gr.gene) + 2 * genic.window, fix = "center")

  if(length(gr.gene) > 0){
    tre.gene <- tre.gene[, !is.na(tre.gene[1, ])]
    genotype.headers <- as.character(utils::read.table(genotype.indexed.f, as.is = TRUE, nrows = 1))
    com.samples <- intersect(colnames(tre.gene), genotype.headers)
    tre.dist <- sQTLseekeR2:::hellingerDist(tre.gene[, com.samples])
    res.df <- data.frame()
    gr.gene.spl <- gr.gene
    res.range <- data.frame()
    genotype.gene <- sQTLseekeR2:::read.bedix(genotype.indexed.f, gr.gene.spl)
    if(!is.null(genotype.gene)){
      snps.to.keep <- sQTLseekeR2:::check.genotype(genotype.gene[,com.samples], tre.gene[, com.samples])
      if(any(snps.to.keep == "PASS")){
          genotype.gene <- genotype.gene[snps.to.keep == "PASS", ]
          res <- LDfilter(genotype.gene, com.samples, tre.dist, th = 1)
          
          rownames(genotype.gene)<-NULL
          print(all.equal(res,genotype.gene))
          if (any(!is.na(res$LD))){
            ld <- unlist(res[!is.na(res$LD),"LD"])
            print(ld)
            cat("n",nrow(genotype.gene),"\n")
            cat("r",length(unlist(strsplit(ld,", ")))/nrow(genotype.gene),"\n")
          }
          
      }
    }
  }
  W <- W + 1 
}    

 