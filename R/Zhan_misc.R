### Trace function:
tr=function(x){return(sum(diag(x))) }

### Required library
## library(PearsonDS)

### KRV test function
### Input: Two n by n kernel matrices, where n>6
### Ouput: KRV test p-value


#' KRV test for overall association between host gene expression and microbiome composition from Zhan et al. (2017)
#'
#' Broadly applicable to test association between two sets of outcomes. We slightly modify the codes to remove the dependence
#' on PearsonDS library: simply use the builtin pgamma() function to compute p-value. 
#' 
#' @param K kernel matrix for a set of outcomes (no need to be centered)
#' @param L kernel matrix for another set of outcomes
#'
#' @return KRV p-value
#' @export
#' @references
#' Zhan, X., Plantinga, A., Zhao, N., Wu, M.C. (2017). A fast small-sample kernel independence test for microbiome community-level association analysis. \emph{Biometrics}, DOI: 10.1111/biom.12684.
#'
#' Original R codes: \url{http://onlinelibrary.wiley.com/wol1/doi/10.1111/biom.12684/suppinfo}
#'
#' Guo,B. and Wu,B. (2017) On the fast small-sample kernel independence test for microbiome community-level association analysis. Biometrics, doi:10.1111/biom.12823.
KRVz <- function(K,L){
  n=nrow(K)
  I.n=diag(1,n)
  I.1=rep(1,n)
  H=I.n-I.1%*%t(I.1)/n
  K=H%*%K%*%H
  L=H%*%L%*%H
  A=K/tr(K%*%K)  ## standard-version of K
  W=L/tr(L%*%L)
  
  Fstar=tr(A%*%W)
  mean.krv=tr(A)*tr(W)/(n-1)	## mean of KRV 
  
  T=tr(A);T2=tr(A%*%A);S2=sum(diag(A)^2)
  Ts=tr(W);T2s=tr(W%*%W);S2s=sum(diag(W)^2)
  temp1=2*((n-1)*T2-T^2)*((n-1)*T2s-Ts^2)/(n-1)^2/(n+1)/(n-2)
  temp21=n*(n+1)*S2- (n-1)*(T^2+2*T2)
  temp22=n*(n+1)*S2s- (n-1)*(Ts^2+2*T2s)
  temp23=(n+1)*n*(n-1)*(n-2)*(n-3)
  temp2=temp21*temp22/temp23
  variance.krv=temp1+temp2		## variance of KRV

  T3=tr(A%*%A%*%A);S3=sum(diag(A)^3);U=sum(A^3);R=t(diag(A))%*%diag(A%*%A);B=t(diag(A))%*%A%*%diag(A)
  T3s=tr(W%*%W%*%W);S3s=sum(diag(W)^3);Us=sum(W^3);Rs=t(diag(W))%*%diag(W%*%W);Bs=t(diag(W))%*%W%*%diag(W)
  t1=n^2*(n+1)*(n^2+15*n-4)*S3*S3s
  t2=4*(n^4-8*n^3+19*n^2-4*n-16)*U*Us
  t3=24*(n^2-n-4)*(U*Bs+B*Us)
  t4=6*(n^4-8*n^3+21*n^2-6*n-24)*B*Bs
  t5=12*(n^4-n^3-8*n^2+36*n-48)*R*Rs
  t6=12*(n^3-2*n^2+9*n-12)*(T*S2*Rs+R*Ts*S2s)
  t7=3*(n^4-4*n^3-2*n^2+9*n-12)*T*Ts*S2*S2s
  t81=(n^3-3*n^2-2*n+8)*(R*Us+U*Rs);t82=(n^3-2*n^2-3*n+12)*(R*Bs+B*Rs)
  t8=24*(t81+t82)
  t9=12*(n^2-n+4)*(T*S2*Us+U*Ts*S2s)
  t10=6*(2*n^3-7*n^2-3*n+12)*(T*S2*Bs+B*Ts*S2s)
  t11=-2*n*(n-1)*(n^2-n+4)*((2*U+3*B)*S3s+(2*Us+3*Bs)*S3)
  t12=-3*n*(n-1)^2*(n+4)*((T*S2+4*R)*S3s+(Ts*S2s+4*Rs)*S3)
  t13=2*n*(n-1)*(n-2)*((T^3+6*T*T2+8*T3)*S3s+(Ts^3+6*Ts*T2s+8*T3s)*S3)
  t14=T^3*((n^3-9*n^2+23*n-14)*Ts^3+6*(n-4)*Ts*T2s+8*T3s)
  t15=6*T*T2*((n-4)*Ts^3+(n^3-9*n^2+24*n-14)*Ts*T2s+4*(n-3)*T3s)
  t16=8*T3*(Ts^3+3*(n-3)*Ts*T2s+(n^3-9*n^2+26*n-22)*T3s)
  t17=-16*(T^3*Us+U*Ts^3)-6*(T*T2*Us+U*Ts*T2s)*(2*n^2-10*n+16)
  t18=-8*(T3*Us+U*T3s)*(3*n^2-15*n+16)-(T^3*Bs+B*Ts^3)*(6*n^2-30*n+24)
  t19=-6*(T*T2*Bs+B*Ts*T2s)*(4*n^2-20*n+24)-8*(T3*Bs+B*T3s)*(3*n^2-15*n+24)
  t201=24*(T^3*Rs+R*Ts^3)+6*(T*T2*Rs+R*Ts*T2s)*(2*n^2-10*n+24)
  t202=8*(T3*Rs+R*T3s)*(3*n^2-15*n+24)+(3*n^2-15*n+6)*(T^3*Ts*S2s+T*S2*Ts^3)
  t203=6*(T*T2*Ts*S2s+Ts*T2s*T*S2)*(n^2-5*n+6)+48*(T3*Ts*S2s+T3s*T*S2)
  t20=-(n-2)*(t201+t202+t203)
  temp31=t1+t2+t3+t4+t5+t6+t7+t8+t9+t10+t11+t12+t13+t14+t15+t16+t17+t18+t19+t20
  temp32=n*(n-1)*(n-2)*(n-3)*(n-4)*(n-5)
  mom3=temp31/temp32
  skewness.krv= (mom3-3*mean.krv*variance.krv-mean.krv^3)/variance.krv^1.5 ## skewness of KRV
  
  m1=mean.krv
  m2=variance.krv
  m3=skewness.krv
  shape=4/m3^2
  scale=sqrt(m2)*m3/2
  location=m1-2*sqrt(m2)/m3
  ## PIIIpars=list(shape,location,scale)
  ## pv=1-ppearsonIII(Fstar, params=PIIIpars)
  ## return(pv)
  
  ## Modified 4/17/2017 by B.Wu
  ssgn = sign(scale)
  gscale = abs(scale)
  ## pv = pgamma(Fstar-location,shape=shape,scale=scale,lower=FALSE)
  pv = pgamma(ssgn * (Fstar - location), shape = shape, scale = gscale, lower.tail=(scale<0))
  return(pv)
}

