% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/varComp.R
\name{varComp}
\alias{varComp}
\title{varComp}
\usage{
varComp(Y, K, numThreads, outPath, outName)
}
\arguments{
\item{Y}{is the phenotype matrix, phenotypes x individual}

\item{K}{is obtained from the Kinship function}

\item{numThreads}{Number of parallel processes}

\item{outPath}{is a parameter that specifies the path of the result file}

\item{outName}{is a parameter that specifies the name of the result file}
}
\value{
Variance components Vg: genetic factor Ve: environment factor
}
\description{
varComp function will estimate variance components
}
\examples{
   Y <- as.matrix(read.table(GeneExpressionData))
   K <- Kinship(X, outPath, outName)

   VC <- varComp(Y, K, numThreads, outPath, outName)
}
