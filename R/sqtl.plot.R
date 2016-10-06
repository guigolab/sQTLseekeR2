##' For a specific gene and SNP, this function produce the graph with the transcript relative expression in the different genotype groups. It also shows the structure of the transcripts, i.e. the CDS/UTR configuration.
##'
##' The format of the transcript structure 'tr.str' should be : 'transId' the transcript ID, 'strand' the DNA strand, 'cdsStart' the start positions of the CDS regions separated by ",", 'cdsEnd' same with the end positions, 'urtStarts' and 'utrEnds' same with UTR regions. Look at the 'SplicingEventClassification' tutorial in the 'docs' folder of the github for some help on how to format such a file.
##'
##' @title Graphs for one sQTL
##' @param geneId the gene ID
##' @param snpId the SNP ID
##' @param tre.df the data.frame with the transcript relative expression (produced by 'prepare.trans.exp').
##' @param genotype.f the path to the genotype file (should be compressed and indexed by 'index.genotype').
##' @param gene.loc a data.frame with the genes location. Columns 'chr', 'start',
##' 'end' and 'geneId' are required.
##' @param genic.window the window(bp) around the gene in which the SNPs are tested. Default is 5000 (i.e. 5kb).
##' @param max.trans the maximum number of transcript shown. The least expressed transcript will be joined for visibility. Default is 6.
##' @param tr.str the data.frame with the structure of the transcripts. If NULL, not used.
##' @return a list of ggplot graphs.
##' @author Jean Monlong
##' @export
sqtl.plot <- function(geneId, snpId, tre.df, genotype.f, gene.loc, genic.window=5e3, max.trans=6, tr.str=NULL){
    if(!file.exists(paste0(genotype.f, ".tbi"))){
        genotype.f = paste0(genotype.f, ".bgz")
    }
    if(!file.exists(genotype.f)){
        stop("Can't find the genotype file: ", genotype.f)
    }

    ## Uglily appeases R checks
    trId = tre = type = strand = . = NULL

    splitTS <- function(df){
        cds = data.frame(type="CDS",
            start=as.integer(unlist(strsplit(df$cdsStarts, ","))),
            end = as.integer(unlist(strsplit(df$cdsEnds, ","))),
            stringsAsFactors=FALSE)
        utr = data.frame(type="UTR",
            start = as.integer(unlist(strsplit(df$utrStarts,","))),
            end = as.integer(unlist(strsplit(df$utrEnds,","))),
            stringsAsFactors=FALSE)
        data.frame(df[,c("trId","geneId","strand")], rbind(cds, utr))
    }

    gene.tre = tre.df[which(tre.df$geneId == geneId),]
    gtre.df = tidyr::gather(gene.tre, "sample", "tre", 3:ncol(gene.tre))

    if(length(unique(gtre.df$trId))>max.trans){
        gtre.rk = aggregate(tre~trId, data=gtre.df, mean)
        top.tr = gtre.rk$trId[which(rank(-gtre.rk$tre)<=max.trans)]
        gtre.df$trId[which(!(gtre.df$trId %in% top.tr))] = "others"
    }

    gene.gr = GenomicRanges::makeGRangesFromDataFrame(gene.loc[which(gene.loc$geneId == geneId),])
    gene.gr = GenomicRanges::resize(gene.gr, GenomicRanges::width(gene.gr)+2*genic.window, fix="center")
    genotype.gene = read.bedix(genotype.f, gene.gr)

    if(all(genotype.gene$snpId != snpId)){
        stop("Can't find SNP in this range and this file. Check 'genic.window' parameter.")
    }

    snp.geno = genotype.gene[which(genotype.gene$snpId == snpId),]
    snp.df = tidyr::gather(snp.geno, "sample", "geno", 5:ncol(snp.geno))
    gp.df = merge(snp.df, gtre.df)

    pdf.l = list()
    pdf.l$relexpression = ggplot2::ggplot(gp.df, ggplot2::aes(x=trId, y=tre, fill=trId)) + ggplot2::geom_boxplot() + ggplot2::facet_grid(.~geno) + ggplot2::theme_bw() + ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1)) + ggplot2::guides(fill=FALSE) + ggplot2::xlab("transcript") + ggplot2::ylab("relative expression") + ggplot2::ylim(0,1)

    if(!is.null(tr.str)){
        gene.str = tr.str[which(tr.str$geneId==geneId),]
        gene.str = dplyr::do(dplyr::group_by(gene.str, geneId, trId, strand), splitTS(.))

        pdf.l$structure = ggplot2::ggplot(gene.str, ggplot2::aes(x=factor(start), xend=factor(end), y=trId, yend=trId, size=type, colour=trId)) + ggplot2::geom_segment() + ggplot2::theme_bw() + ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, hjust=1)) + ggplot2::guides(colour=FALSE) + ggplot2::scale_size_manual(values=3:2) + ggplot2::xlab("position") + ggplot2::ggtitle(paste(gene.str$strand[1], "strand"))
    }

    return(pdf.l)
}

