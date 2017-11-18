# MiRV
The developed R package MiRV is for testing the association of multiple outcomes with the microbiome profiles.

You can install the package with devtools::install_github("baolinwu/MiRV").


# Reference
Guo,B. and Wu,B. (2017) On the fast small-sample kernel independence test for microbiome community-level association analysis. *Biometrics*, doi:10.1111/biom.12823.


# Usage
-----

``` r
library(vegan)
library(MiRV)
## Simulate microbiome profiles of 1000 OTUs with Dirich-Multinom
## simulate Dirichlet from Gamma dist
Mi = matrix(NA,200,1e3)
alpha = rchisq(1e3,30)
Gam = matrix(rgamma(200*1e3, shape=alpha), 200,1e3, byrow=TRUE)
ak = Gam/rowSums(Gam)
for(i in 1:200){
  Mi[i,] = rmultinom(1,1000,ak[i,])
}
## Simulate host gene expression profiles of 1000 genes
Y = matrix(rnorm(200*1e3), 200,1e3)
Z = scale(rowSums(Mi[,1:20]),TRUE,TRUE)
for(i in 1:10) Y[,i] = Y[,i] + rnorm(1)*Z*1.5
## compute dist and kernel
Dx = as.matrix(dist(scale(Y,TRUE,TRUE)))
s2 = median(Dx[upper.tri(Dx)])
Kx = exp(-Dx/s2)
Kx = D2Ks(Dx)$K
K = K2C(Kx)
Dm = as.matrix(vegdist(Mi))
Lm = D2Ks(Dm)$K
L = K2C(Lm)
KRVz(K,L)
KRVtest(K,L,TRUE)
GRVtest(K,D2Gs(Dm),TRUE)
RVtest(K,L)
```
