#' Create a centered Kernel matrix
#'
#' Given a symmetric kernel matrix K, we calculate \eqn{HKH}, where \eqn{H=I-11'/n} is a centering/projection matrix.
#' Mathematically it is equivalent to subtract the row and column means, which is computationally more efficient.
#' @param K kernel matrix
#' @return row and column centered kernel matrix
#' @export
#' @references
#' Guo,B. and Wu,B. (2017) On the fast small-sample kernel independence test for microbiome community-level association analysis. Biometrics, doi:10.1111/biom.12823.
K2C <- function(K){
  a1 = rowMeans(K)
  t(K-a1)-a1 + mean(a1)
}

#' Convert pairwise distance matrix to Gower matrix
#'
#' Adapted from the GRV R codes (Minas et al., 2013). Modified for better efficiency avoiding unnecessary matrix products.
#' Can show \eqn{G_{ij}= -0.5(D_{ij}^2-\bar{D^2}_{i,}-\bar{D^2}_{.,j}+\bar{D^2})}, where \eqn{\bar{D^2}_{i,},\bar{D^2}_{,j},\bar{D^2}} are
#' row, column and global averages. For symmetric distance matrix, \eqn{\bar{D^2}_{i,}=\bar{D^2}_{,i}}.
#' Note that G could be non-positive definite, and no further transformations are applied to G to make it non-negative definite.
#' 
#' @param D pairwise distance matrix
#' @return G centered Gower matrix
#' 
#' @export
#' @references
#' Minas, C., Curry, E., Montana, G., 2013. A distance-based test of association between paired heterogeneous genomic data. Bioinformatics 29, 2555–2563. \url{http://wwwf.imperial.ac.uk/~gmontana/grv.html}
#' 
#' Guo,B. and Wu,B. (2017) On the fast small-sample kernel independence test for microbiome community-level association analysis. Biometrics, doi:10.1111/biom.12823.
D2Gs <- function(D){
  D2 = as.matrix(D)^2
  Da = colMeans(D2); Dm = mean(Da)
  G = -0.5*( t(D2 - Da) - Da + Dm )
  return(G)
}

#' Convert pairwise distance matrix to kernels
#'
#' Adapted from the MiRKAT package (Zhao et al., 2015). Modified for better efficiency avoiding unnecessary matrix products, and output
#' kernel and its eigen decomposition for followup analysis.
#' Can show \eqn{K_{ij}= -0.5(D_{ij}^2-\bar{D^2}_{i,}-\bar{D^2}_{.,j}+\bar{D^2})}, where \eqn{\bar{D^2}_{i,},\bar{D^2}_{,j},\bar{D^2}} are
#' row, column and global averages. For symmetric distance matrix, \eqn{\bar{D^2}_{i,}=\bar{D^2}_{,i}}.
#' 
#' @param D pairwise distance matrix
#' @return
#' \describe{
#'   \item{K}{ computed (non-negative) definite kernel matrix}
#'   \item{Kh}{ computed from eigen decomposition of K, such that \eqn{K=K_hK_h^T} }
#'   \item{Km}{ the effective rank of K }
#' }
#' 
#' @export
#' @references
#' Chen, J., Chen, W., Zhao, N., Wu, M.C., Schaid, D.J., 2016. Small Sample Kernel Association Tests for Human Genetic and Microbiome Association Studies. Genet. Epidemiol. 40, 5--19.
#' 
#' Zhao,N. et al. (2015) Testing in Microbiome-Profiling Studies with MiRKAT, the Microbiome Regression-Based Kernel Association Test. \emph{AJHG}, 96(5): 797–807.
#' 
#' Guo,B. and Wu,B. (2017) On the fast small-sample kernel independence test for microbiome community-level association analysis. Biometrics, doi:10.1111/biom.12823.
D2Ks <- function(D){
  ## n = dim(D)[1]
  D2 = D^2
  Da = colMeans(D2); Dm = mean(Da)
  K = -0.5*( t(D2 - Da) - Da + Dm )
  ## centerM = diag(n ) - 1/n
  ## K = -0.5*centerM %*% (D*D) %*% centerM
  ## eK= eigen(K, symmetric = T)
  ## K = eK$vector %*% diag(abs(eK$values)) %*% t(eK$vector)
  ek = eigen(K, sym=TRUE)
  val = zapsmall(abs(ek$val))
  id = which(val>0); Km = length(id)
  Gt = t(ek$vec[,id])*sqrt(val[id]) 
  K = t(Gt)%*%Gt
  return( list(K=K, Kh=t(Gt), Km=Km ) )
}


## Compute p-value for GRV/KRV tests based on Gamma dist approximation
surv.GRV <- function(GX,GY, y){
  n = dim(GX)[1]
  ## constant = norm.matrix(GX)*norm.matrix(GY)
  constant = sqrt(sum(GX^2)*sum(GY^2))
  perm.info = perm.info.tr.AB(GX,GY,n)
  m = perm.info$mean.T
  s = sqrt(perm.info$variance.T)
  g = perm.info$skewness.T
  if(g >= 0){
    res = pgamma((y*constant - m)/s+2/g, shape=4/(g^2), scale = g/2, lower=FALSE)
  }else{
    res = pgamma(abs(2/g) -(y*constant - m)/s, shape=4/(g^2), scale = abs(g)/2,lower.tail=TRUE)
  }
  return(res)
}

#' Compute the kernel independence test statistic between two kernels efficiently
#'
#' For symmetric kernel matrices, we can efficiently compute the kernel test statistic as a form of correlation,
#' \eqn{(\sum_{i,j} X_{ij}Y_{ij})/\sqrt{(\sum_{ij}X_{ij}^2)(\sum_{ij}Y_{ij}^2)} }.
#' 
#' @param KX centered symmetric kernel matrix for a set of outcomes 
#' @param KY centered symmetric kernel matrix for another set of outcomes
#' @return kernel test statistic
#'
#' @export
#' @references
#' Minas, C., Curry, E., Montana, G., 2013. A distance-based test of association between paired heterogeneous genomic data. Bioinformatics 29, 2555–2563.
#' 
#' Zhan, X., Plantinga, A., Zhao, N., Wu, M.C. (2017). A fast small-sample kernel independence test for microbiome community-level association analysis. Biometrics, doi:10.1111/biom.12684.
#'
#' Guo,B. and Wu,B. (2017) On the fast small-sample kernel independence test for microbiome community-level association analysis. Biometrics, doi:10.1111/biom.12823.
KRV <- function(KX,KY){
  sum(KX*KY)/sqrt(sum(KX^2)*sum(KY^2))
}

#' Revised implementation of GRV test for independence between paired data
#'
#' @param GX distance matrix or the computed Gower matrix for a set of outcomes
#' @param GY distance matrix or the computed Gower matrix for another set of outcomes
#' @param kernel indicator of whether the input is Gower matrix computed using D2Gs().
#'
#' @return test statistic and its p-value
#' @export
#' @references
#' Minas, C., Curry, E., Montana, G., 2013. A distance-based test of association between paired heterogeneous genomic data. Bioinformatics 29, 2555–2563.
#' 
#' Guo,B. and Wu,B. (2017) On the fast small-sample kernel independence test for microbiome community-level association analysis. Biometrics, doi:10.1111/biom.12823.
GRVtest <- function(GX,GY, Gower=FALSE){
  if(!Gower){
    GX = D2Gs(as.matrix(GX))
    GY = D2Gs(as.matrix(GY))
  }
  Q = KRV(GX,GY)
  pval = surv.GRV(GX,GY, Q)
  return(c(statistic=Q, p.value=pval))
}

#' Revised implementation of KRV test for independence between paired data
#'
#' @param KX distance matrix or the computed kernel matrix for a set of outcomes
#' @param KY distance matrix or the computed kernel matrix for another set of outcomes
#' @param kernel indicator of whether the input is kernel matrix computed using D2Ks() and K2C().
#'
#' @return test statistic and its p-value
#' @export
#' @references
#' Minas, C., Curry, E., Montana, G., 2013. A distance-based test of association between paired heterogeneous genomic data. Bioinformatics 29, 2555–2563.
#' 
#' Zhan, X., Plantinga, A., Zhao, N., Wu, M.C. (2017). A fast small-sample kernel independence test for microbiome community-level association analysis. Biometrics, doi:10.1111/biom.12684.
#'
#' Guo,B. and Wu,B. (2017) On the fast small-sample kernel independence test for microbiome community-level association analysis. Biometrics, doi:10.1111/biom.12823.
KRVtest <- function(KX,KY, kernel=FALSE){
  if(!kernel){
    KX = K2C(D2Ks(as.matrix(KX))$K)
    KY = K2C(D2Ks(as.matrix(KY))$K)
  }
  Q = KRV(KX,KY)
  pval = surv.GRV(KX,KY, Q)
  return(c(statistic=Q, p.value=pval))
}



