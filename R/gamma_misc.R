## Calculate the first three permutation moments
## @references
## Kazi-Aoual, F., Hitier, S., Sabatier, R., Lebreton, J.-D. (1995). Refined Approximations to Permutation Tests for Multivariate Inference. Comput. Stat. Data Anal. 20, 643–656.
perm.mom.AB <- function(A,B){
  n = dim(A)[1]
  perm.1st.mom <- function(A1,B1)  return(A1*B1/(n-1))
  perm.2nd.mom <- function(A1,A2,A3,B1,B2,B3){
    p11 <- ((n-1)*A3*B3+(B1^2-B3)*(A1^2-A3)+2*(A2-A3)*(B2-B3)+4*A3*B3)/(n*(n-1))
    p22 <- (4*(n-3)*(2*A3-A2)*(2*B3-B2)+2*(2*A3-A1^2)*(2*B3-B1^2)*(n-3)+(2*A2+A1^2-6*A3)*(2*B2+B1^2-6*B3))/(n*(n-1)*(n-2)*(n-3))
    return(p11+p22)
  }
  perm.3rd.mom <- function(A1,A2,A3,A4,A5,A6,A7,A8,B1,B2,B3,B4,B5,B6,B7,B8){
    n2 <- n^2
    n3 <- n^3
    n4 <- n^4
    return((n2*(n+1)*(n2+15*n-4)*A5*B5 + 4*(n4-8*n3+19*n2-4*n-16)*A6*B6 + 24*(n2-n-4)*(A6*B8+B6*A8)+ 6*(n4-8*n3+21*n2-6*n-24)*A8*B8 + 12*(n4-n3-8*n2+36*n-48)*A7*B7 + 12*(n3-2*n2+9*n-12)*(A1*A3*B7 + A7*B1*B3) + 3*(n4 - 4*n3 - 2*n2+9*n-12)*A1*B1*A3*B3 + 24*((n3 - 3*n2 - 2*n+8)*(A7*B6 + A6*B7) + (n3 - 2*n2 - 3*n+12)*(A7*B8 + A8*B7)) + 12*(n2 - n + 4)*(A1*A3*B6 + B1*B3*A6) + 6*(2*n3 - 7*n2 - 3*n + 12)*(A1*A3*B8 + A8*B1*B3) - 2*n*(n-1)*(n2-n+4)*((2*A6+3*A8)*B5+(2*B6+3*B8)*A5) - 3*n*(n-1)*(n-1)*(n+4)*((A1*A3+4*A7)*B5+(B1*B3+4*B7)*A5) + 2*n*(n-1)*(n-2)*((A1^3 + 6*A1*A2 + 8*A4)*B5+(B1^3 + 6*B1*B2 + 8*B4)*A5) + (A1^3)*((n3-9*n2+23*n-14)*(B1^3)+6*(n-4)*B1*B2+8*B4) + 6*A1*A2*((n-4)*(B1^3)+(n3-9*n2+24*n-14)*B1*B2 + 4*(n-3)*B4) + 8*A4*((B1^3)+3*(n-3)*B1*B2+(n3-9*n2+26*n-22)*B4) - 16*((A1^3)*B6+A6*(B1^3)) - 6*(A1*A2*B6 + A6*B1*B2)*(2*n2-10*n+16) - 8*(A4*B6+A6*B4)*(3*n2-15*n+16)-((A1^3)*B8+A8*(B1^3))*(6*n2-30*n+24)-6*(A1*A2*B8+A8*B1*B2)*(4*n2-20*n+24) - 8*(A4*B8+A8*B4)*(3*n2-15*n+24) - (n-2)*(24*((A1^3)*B7+A7*(B1^3))+6*(A1*A2*B7+A7*B1*B2)*(2*n2-10*n+24)+8*(A4*B7+A7*B4)*(3*n2-15*n+24)+(3*n2-15*n+6)*((A1^3)*B1*B3+A1*A3*(B1^3))+6*(A1*A2*B1*B3+A1*A3*B1*B2)*(n2-5*n+6) + 48*(A4*B1*B3+A1*A3*B4)))/(n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5)))
  }
  Am2 <- A%*%A
  Am3 <- A%*%Am2
  dA <- diag(A)
  Ap2 <- A^2
  Ap3 <- Ap2*A
  A1 <- sum(dA)
  A2 = sum(A^2) ## A2 <- sum(diag(Am2))
  A3 = sum(dA^2) ##  A3 <- sum(diag(Ap2))
  A4 <- sum(diag(Am3))
  A5 = sum(dA^3) ## A5 <- sum(diag(Ap3))
  A6 = sum(Ap3) ## A6 <- sum(Ap3)
  A7 = sum(dA*diag(Am2)) ## A7 <- dA%*%diag(Am2)
  A8 = sum(t(dA*A)*dA) ## A8 <- dA%*%A%*%dA
  
  Bm2 <- B%*%B
  Bm3 <- B%*%Bm2
  dB <- diag(B)
  Bp2 <- B^2
  Bp3 <- Bp2*B
  B1 <- sum(dB)
  B2 <- sum(diag(Bm2))
  B3 <- sum(diag(Bp2))
  B4 <- sum(diag(Bm3))
  B5 <- sum(diag(Bp3))
  B6 <- sum(B^3)
  B7 <- dB%*%diag(Bm2)
  B8 <- dB%*%B%*%dB
  mom1 <- perm.1st.mom(A1, B1)
  mom2 <- perm.2nd.mom(A1,A2,A3,B1,B2,B3)
  mom3 <- perm.3rd.mom(A1,A2,A3,A4,A5,A6,A7,A8,B1,B2,B3,B4,B5,B6,B7,B8)
  mean.T <- mom1
  variance.T <- (mom2 - mom1^2)
  skewness.T <- (mom3 - 3 * variance.T * mom1 - mom1^3)/(variance.T^(3/2))
  return(list(mean.T=mean.T,variance.T=variance.T,skewness.T=skewness.T, moms=c(mom1,mom2,mom3)))
}

perm.C4 <- function(KX,KY, N=1e5){
  n = dim(KX)[1]
  Qb = rep(0, N)
  for(i in 1:N){
    ib = sample(1:n, size=n, rep=FALSE)
    Qb[i] = sum(KX*KY[ib,ib])
  }
  C4 = mean((Qb-mean(Qb))^4)
  return(list(C4=C4,Qb=Qb))
}



##' Calculate the first four permutation moments: mean, variance, skewness and kurtosis for RV statistic
##'
##' @param KX the first kernel matrix
##' @param KY the second kernel matrix
##'
##' @return a vector of mean, standard error, skewness, and kurtosis
##' @export
##' @references
##' Kazi-Aoual, F., Hitier, S., Sabatier, R., Lebreton, J.-D. (1995). Refined Approximations to Permutation Tests for Multivariate Inference. Comput. Stat. Data Anal. 20, 643–656.
##'
##' Guo,B. and Wu,B. (2017) On the fast small-sample kernel independence test for microbiome community-level association analysis. Biometrics, doi:10.1111/biom.12823.
RVparam <- function(KX,KY){
  const = sqrt(sum(KX^2)*sum(KY^2))
  aa = perm.mom.AB(KX,KY)
  s = aa$skewness
  mu = aa$mean/const
  sig = sqrt(aa$var)/const
  a0 = perm.C4(KX,KY, N=1e5)
  r = a0$C4/aa$var^2
  return(c(mu=mu,sig=sig,s=s,r=r))
}

#' Approximately model RV statistic Q with a shifted gamma distribution Gamma(k,\eqn{\theta})
#'
#' The idea is to model \eqn{Q-\lambda} with Gamma(k,\eqn{\theta}). We choose (k,\eqn{\theta},\eqn{\lambda}) to match
#' the mean, variance and skewness or kurtosis of Q (all computed under permutation distribution).
#' The Gamma PDF \eqn{\frac{1}{\Gamma(k)\theta^k}x^{k-1}e^{-x/\theta}}, mean \eqn{k\theta}, variance \eqn{k\theta^2},
#' skewness \eqn{2/\sqrt{k}}, kurtosis \eqn{6/k+3}.
#'
#' @param mu mean of permuted Q
#' @param sigma standard error of permuted Q
#' @param s skewness: \eqn{E(Q-\mu)^3/\sigma^3}
#' @param r kurtosis: \eqn{E(Q-\mu)^4/\sigma^4}
#' @return
#' \describe{
#'   \item{par3}{estimated (k,\eqn{\theta},\eqn{\lambda}) when matching skewness}
#'   \item{par4}{estimated (k,\eqn{\theta},\eqn{\lambda}) when matching kurtosis}
#' }
#' 
#' @export
#' @references
#' Minas, C., Curry, E., Montana, G., 2013. A distance-based test of association between paired heterogeneous genomic data. Bioinformatics 29, 2555–2563.
#' 
#' Zhan, X., Plantinga, A., Zhao, N., Wu, M.C. (2017). A fast small-sample kernel independence test for microbiome community-level association analysis. Biometrics, doi:10.1111/biom.12684.
#'
#' Guo,B. and Wu,B. (2017) On the fast small-sample kernel independence test for microbiome community-level association analysis. Biometrics, doi:10.1111/biom.12823.
RVgamma <- function(mu,sigma,s,r){
  ## match kurtosis
  k1 = 6/abs(r-3)
  t1 = sigma/sqrt(k1)
  lam1 = mu - 6*t1/(r-3)
  ## match skewness
  k2 = 4/s^2
  t2 = sigma*abs(s)/2
  lam2 = mu-2*sigma/s
  return(list(par3=c(k=k2,theta=t2,lambda=lam2), par4=c(k=k1,theta=t1,lambda=lam1)) )
}

#' Unified interface to conduct general RV test using various Gamma approximations
#'
#' Given centered kernel/Gower matrices, test their independence.
#' @param KX the first centered kernel/Gower matrix
#' @param KY the second centered kernel/Gower matrix
#'
#' @return a vector of three p-values: Gamma approx by matching skewness and kurtosis respectively, plus the permutation p-value
#' @export
#' 
#' @references
#' Guo,B. and Wu,B. (2017) On the fast small-sample kernel independence test for microbiome community-level association analysis. Biometrics, doi:10.1111/biom.12823.
RVtest <- function(KX,KY){
  const = sqrt(sum(KX^2)*sum(KY^2))
  Q = sum(KX*KY)/const
  aa = perm.mom.AB(KX,KY)
  s = aa$skewness
  mu = aa$mean/const
  sig = sqrt(aa$var)/const
  a0 = perm.C4(KX,KY, N=1e4)
  ppval = mean((a0$Qb/const)>Q)
  r = a0$C4/aa$var^2
  gk = RVgamma(mu,sig,s,r)
  if(s>0){
    pval3 = pgamma(Q-gk$par3[3],shape=gk$par3[1],scale=gk$par3[2],lower=FALSE)
  } else{
    pval3 = pgamma(-(Q-gk$par3[3]),shape=gk$par3[1],scale=gk$par3[2],lower=TRUE)
  }
  pval4 = pval3
  if(r>3){
    pval4 = pgamma(Q-gk$par4[3],shape=gk$par4[1],scale=gk$par4[2],lower=FALSE)
  } else if(r<3){
    pval4 = pgamma(-(Q-gk$par4[3]),shape=gk$par4[1],scale=gk$par4[2],lower=TRUE)
  }
  return(list(p.value=c(pval3=pval3,pval4=pval4,ppval=ppval), r=r) )
  ## return(c(pval3=pval3,pval4=pval4,ppval=ppval) )
}
