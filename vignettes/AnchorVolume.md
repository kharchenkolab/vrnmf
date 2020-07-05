---
title: "AnforFree algorithm"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{AnchorVolume}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



Generally, structured matrix facrotization problem is formulated as decomposition of original matrix $X$ in a product of unknown non-negative matrix $C$ and matrix $D$ of lower rank: 

$$
X = C \cdot D
$$
In practice, if $X$ has significantly more rows than columns or if observations are noisy and represent statistical sampling than reformulation of the problem in "covariance domain" simplifies the problem:

$$
Y = C \cdot E \cdot C^{T},
$$
where $Y = X \cdot X^{T}, E = D \cdot D^{T}.$

The goal of _AnchorFree_ approach is to tri-factorize matrix $Y$ in a product of matrices $C$ and $E$. Under a relatively mild assumption on spread of column vectors of matrix $C$, matrix $D$ has minimum volume across all possible factorizations. _AnchorFree_ seeks to find a factorization pair of matrices $(C,D)$ such that matrix $D$ would have minimum determinant (as a proxy of volume) using alternating linear programming. 

Application of the method is examplified below for a simulated dataset. For that, we will use simulated datasets with matrices `X`, its noisy version `X.noise` and decomposition matrices `C` and `D` such that `X = C*D`. To simulate dataset with parameters of interest one can use function `sim_factors()`.



```r
library(vrnmf)
#data(simnmf)
```