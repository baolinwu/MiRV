% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/D2GK.R
\name{GRVtest}
\alias{GRVtest}
\title{Revised implementation of GRV test for independence between paired data}
\usage{
GRVtest(GX, GY, Gower = FALSE)
}
\arguments{
\item{GX}{distance matrix or the computed Gower matrix for a set of outcomes}

\item{GY}{distance matrix or the computed Gower matrix for another set of outcomes}

\item{kernel}{indicator of whether the input is Gower matrix computed using D2Gs().}
}
\value{
test statistic and its p-value
}
\description{
Revised implementation of GRV test for independence between paired data
}
\references{
Minas, C., Curry, E., Montana, G., 2013. A distance-based test of association between paired heterogeneous genomic data. Bioinformatics 29, 2555–2563.

Guo,B. and Wu,B. (2017) On the fast small-sample kernel independence test for microbiome community-level association analysis. Biometrics, doi:10.1111/biom.12823.
}
