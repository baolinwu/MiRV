% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Zhan_misc.R
\name{KRVz}
\alias{KRVz}
\title{KRV test for overall association between host gene expression and microbiome composition from Zhan et al. (2017)}
\usage{
KRVz(K, L)
}
\arguments{
\item{K}{kernel matrix for a set of outcomes (no need to be centered)}

\item{L}{kernel matrix for another set of outcomes}
}
\value{
KRV p-value
}
\description{
Broadly applicable to test association between two sets of outcomes. We slightly modify the codes to remove the dependence
on PearsonDS library: simply use the builtin pgamma() function to compute p-value.
}
\references{
Zhan, X., Plantinga, A., Zhao, N., Wu, M.C. (2017). A fast small-sample kernel independence test for microbiome community-level association analysis. \emph{Biometrics}, DOI: 10.1111/biom.12684.

Original R codes: \url{http://onlinelibrary.wiley.com/wol1/doi/10.1111/biom.12684/suppinfo}

Guo,B. and Wu,B. (2017) On the fast small-sample kernel independence test for microbiome community-level association analysis. Biometrics, doi:10.1111/biom.12823.
}
