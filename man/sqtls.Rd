% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sqtls.R
\name{sqtls}
\alias{sqtls}
\title{sQTL retrieval}
\usage{
sqtls(
  res.df,
  FDR = 0.05,
  method = "BH",
  md.min = 0.05,
  out.pdf = NULL,
  svQTL.removal = TRUE,
  FDR.svQTL = 0.05
)
}
\arguments{
\item{res.df}{a data.frame, output of 'sqtl.seeker' with the P-values
for each gene/SNP test.}

\item{FDR}{the False Discovery Rate to call an association significant. Default is 0.05.}

\item{method}{FDR approach: either \code{"BH"} (Benjamini-Hochberg) or \code{"qvalue"} (Storey).
Default is "BH".}

\item{md.min}{the minimum MD (Maximum Difference) in relative expression. Maximum
difference in relative expression (MD) gives an idea of the effect size of the association.
Default is 0.05.}

\item{out.pdf}{the name of the pdf file to create. If NULL (default), no pdf output
will be created. If non-NULL, the distribution of the P-values and a semi-volcano plot
(P-value vs MD) will be shown.}

\item{svQTL.removal}{if TRUE (and column 'pv.svQTL' is present in 'res.df') significant
sQTL which are also significant svQTLs are not reported.}

\item{FDR.svQTL}{the False Discovery Rate to call a svQTL, that may be removed from the final set of sQTLs.}
}
\value{
a subset of the input data.frame with only significant sQTLs and FDR estimates.
}
\description{
Retrieves sQTLs after multiple testing correction and svQTL removal (if requested).
The distribution of the P-values and a semi-volcano plot can be also displayed.
}
\details{
Multiple testing correction can be performed using either Benjamini-Hochberg (\code{\link{p.adjust}} method) 
or Storey (\code{\link{qvalue}} method) approach for FDR. 

If \code{svQTL.removal = TRUE} and svQTLs were tested in 'sqtl.seeker', gene/SNPs with
significant svQTL association (after multiple testing correction and similar FDR threshold)
are removed from the final set of significant sQTLs.
}
\author{
Jean Monlong, Diego Garrido-Martín
}
