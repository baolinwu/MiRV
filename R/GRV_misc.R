## GRV Test
## downloaded from http://wwwf.imperial.ac.uk/~gmontana/grv.html. date 4/11/2017
## Date: 21/8/2012
## Author: Christopher Minas

"trace.matrix"<- function(x){
  return(sum(diag(x)))
}
"norm.matrix" <- function(mat){
  return(sqrt(trace.matrix(mat%*%t(mat))))
}

"create.G" <- function(n,D){
  D <- as.matrix(D)
  A <- -(1/2)*D^2
  ones <- matrix(data=1,ncol=n,nrow=n)
  centmat <- (diag(n)-ones/n)
  G <- centmat%*%A%*%centmat
  return(G)
}


"perm.info.tr.AB" <- function(A,B,n){
  "perm.1st.moment.tr.AB" <- function(A1,B1,n){
    return(A1*B1/(n-1))
  }
  "perm.2nd.moment.tr.AB" <- function(A1,A2,A3,B1,B2,B3,n){
    p11 <- ((n-1)*A3*B3+(B1^2-B3)*(A1^2-A3)+2*(A2-A3)*(B2-B3)+4*A3*B3)/(n*(n-1))
    p22 <- (4*(n-3)*(2*A3-A2)*(2*B3-B2)+2*(2*A3-A1^2)*(2*B3-B1^2)*(n-3)+(2*A2+A1^2-6*A3)*(2*B2+B1^2-6*B3))/(n*(n-1)*(n-2)*(n-3))
    return(p11+p22)
  }
  "perm.3rd.moment.tr.AB" <- function(A1,A2,A3,A4,A5,A6,A7,A8,B1,B2,B3,B4,B5,B6,B7,B8,n){
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
  A2 <- trace.matrix(Am2)
  A3 <- trace.matrix(Ap2)
  A4 <- trace.matrix(Am3)
  A5 <- trace.matrix(Ap3)
  A6 <- sum(Ap3)
  A7 <- dA%*%diag(Am2)
  A8 <- dA%*%A%*%dA
  Bm2 <- B%*%B
  Bm3 <- B%*%Bm2
  dB <- diag(B)
  Bp2 <- B^2
  Bp3 <- Bp2*B
  B1 <- sum(dB)
  B2 <- trace.matrix(Bm2)
  B3 <- trace.matrix(Bp2)
  B4 <- trace.matrix(Bm3)
  B5 <- trace.matrix(Bp3)
  B6 <- sum(B^3)
  B7 <- dB%*%diag(Bm2)
  B8 <- dB%*%B%*%dB
  mom1 <- perm.1st.moment.tr.AB(A1, B1, n)
  mom2 <- perm.2nd.moment.tr.AB(A1,A2,A3,B1,B2,B3,n)
  mom3 <- perm.3rd.moment.tr.AB(A1,A2,A3,A4,A5,A6,A7,A8,B1,B2,B3,B4,B5,B6,B7,B8,n)
  mean.T <- mom1
  variance.T <- (mom2 - mom1^2)
  skewness.T <- (mom3 - 3 * variance.T * mom1 - mom1^3)/(variance.T^(3/2))
  return(list(mean.T=mean.T,variance.T=variance.T,skewness.T=skewness.T))
}

"pdf.GRV" <- function(DX,DY,ys){
  "pearson.type.3.dgamma" <- function(g,y){
    if(g >= 0){
      ans <- dgamma(y+2/g, shape=4/(g^2), scale = g/2, log = FALSE)
    }else{
      ans <- dgamma(abs(2/g) -y, shape=4/(g^2), scale = abs(g)/2, log = FALSE)
    }
    return(ans)
  }
  DX <- as.matrix(DX)
  DY <- as.matrix(DY)
  n <- nrow(DX)
  GX <- create.G(n,DX)
  GY <- create.G(n,DY)
  constant <- norm.matrix(GX)*norm.matrix(GY)
  perm.info <- perm.info.tr.AB(GX,GY,n)
  ## constant <- norm.matrix(GX)*norm.matrix(GY)
  m <- perm.info$mean.T
  s <- sqrt(perm.info$variance.T)
  g <- perm.info$skewness.T
  return((constant/s)*sapply(1:length(ys), function(i) pearson.type.3.dgamma(g,(ys[i]*constant - m)/s)))
}


"cdf.GRV" <- function(GX,GY,n,y){
  constant <- norm.matrix(GX)*norm.matrix(GY)
  perm.info <- perm.info.tr.AB(GX,GY,n)
  m <- perm.info$mean.T
  s <- sqrt(perm.info$variance.T)
  g <- perm.info$skewness.T
  if(g >= 0){
    ans <- pgamma((y*constant - m)/s+2/g, shape=4/(g^2), scale = g/2)
  }else{
    ans <- pgamma(abs(2/g) -(y*constant - m)/s, shape=4/(g^2), scale = abs(g)/2,lower.tail=FALSE)
  }
  return(ans)
}



"GRV" <- function(GX,GY){
  return(trace.matrix(GX%*%GY)/(norm.matrix(GX)*norm.matrix(GY)))
}

#' GRV test for independence beteween paired heterogeneous genomic data from Minas et al. (2013)
#'
#' @param DX distance matrix for a set of outcomes
#' @param DY distance matrix another set of outcomes
#' @return test statistic and its associated p-value
#'
#' @export
#' @references
#' Minas, C., Curry, E., Montana, G., 2013. A distance-based test of association between paired heterogeneous genomic data. Bioinformatics 29, 2555–2563.
"GRV.test" <- function(DX,DY){
  DX <- as.matrix(DX)
  DY <- as.matrix(DY)
  n <- nrow(DX)
  GX <- create.G(n,DX)
  GY <- create.G(n,DY)
  grv0 <- GRV(GX,GY)
  p.grv0 <- 1 - cdf.GRV(GX,GY,n,grv0)
  return(c(grv.statistic=grv0,grv.p.value=p.grv0))
}
