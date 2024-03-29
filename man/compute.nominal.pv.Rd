% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sqtl.seeker.p.R
\name{compute.nominal.pv}
\alias{compute.nominal.pv}
\title{Compute nominal P-value}
\usage{
compute.nominal.pv(
  geno.df,
  tre.mt,
  permute = FALSE,
  seed = 1,
  item.acc = 1e-14,
  eigen.tol = 1e-12
)
}
\arguments{
\item{geno.df}{a data.frame of one row with the genotype information for each sample.}

\item{tre.mt}{a matrix of splicing ratios (samples x transcripts).}

\item{permute}{should the rows of the splicing ratio matrix be permuted. Default is \code{FALSE}.}

\item{seed}{if \code{permute} is TRUE, value provided to \code{\link{set.seed}}. Default is 1.}

\item{item.acc}{accuracy for P-value computation. Passed to \code{pcqf} function. Default is 1e-14.
Minimum accuracy allowed is 1e-14.}

\item{eigen.tol}{eigenvalues below this threshold are considered 0. Default is 1e-12.}
}
\value{
A data.frame containing a P-value for the association.
}
\description{
Computes nominal P-value for the association between the genotype of a locus (e.g. SNP) and the splicing ratios of a gene.
}
\author{
Diego Garrido-Martín
}
\keyword{internal}
