# AnforFree algorithm
 
## Introduction

  Generally, structured matrix factorization problem is formulated as a decomposition of original matrix <img src="https://render.githubusercontent.com/render/math?math=X"> in a product of unknown non-negative matrix <img src="https://render.githubusercontent.com/render/math?math=C"> and matrix <img src="https://render.githubusercontent.com/render/math?math=D"> of lower rank: 
  
<img src="https://render.githubusercontent.com/render/math?math=\begin{aligned}X = C \cdot D \end{aligned}">

  In practice, if <img src="https://render.githubusercontent.com/render/math?math=X"> has significantly more rows than columns or if observations are noisy and represent statistical sampling than reformulation of the problem in "covariance domain" simplifies inference:
  
<img src="https://render.githubusercontent.com/render/math?math=Y = C \cdot E \cdot C^{T},">

  where <img src="https://render.githubusercontent.com/render/math?math=Y = X \cdot X^{T}, E = D \cdot D^{T}.">

  The goal of _AnchorFree_ approach is to tri-factorize matrix <img src="https://render.githubusercontent.com/render/math?math=Y"> in a product of matrices <img src="https://render.githubusercontent.com/render/math?math=C"> and <img src="https://render.githubusercontent.com/render/math?math=E">. Under a relatively mild assumption on spread of column vectors of matrix <img src="https://render.githubusercontent.com/render/math?math=C">, matrix <img src="https://render.githubusercontent.com/render/math?math=E"> has minimum volume across all possible factorizations. _AnchorFree_ seeks to find a factorization pair of matrices <img src="https://render.githubusercontent.com/render/math?math=(C, E)"> such that column vectors of <img src="https://render.githubusercontent.com/render/math?math=C"> represent unit simplex and matrix <img src="https://render.githubusercontent.com/render/math?math=E"> has minimum determinant (as a proxy of volume) using alternating linear programming. 

## Application

Application of the method is examplified below for a simulated dataset. For that, we will use dataset with matrices `X`, its noisy version `X.noise` and decomposition matrices `C` and `D` such that `X = C*D` (to simulate a dataset with parameters of interest one can use function `sim_factors`). First, we upload the library and the dataset:
  
  ```r
  library(vrnmf)
  
  data(simnmf)
  
  names(simnmf)
  ```
  
  ```
  ## [1] "X.noise" "X"       "C"       "D"
  ```

Summary below shows that the original matrix `X` of the size <img src="https://render.githubusercontent.com/render/math?math=192 \cdot 3000"> is a product of two matrices of lower rank 10:
  
  
  ```r
  str(simnmf$X)
  ```
  
  ```
  ##  num [1:192, 1:3000] 0.000493 0.000695 0.000347 0.000248 0.000699 ...
  ```
  
  ```r
  str(simnmf$C)
  ```
  
  ```
  ##  num [1:192, 1:10] 0.00387 0.01099 0.00302 0 0.008 ...
  ```
  
  ```r
  str(simnmf$D)
  ```
  
  ```
  ##  num [1:10, 1:3000] 0.03793 0.001 0.00138 0.00229 0.001 ...
  ```

As a first step, input matrix `X` is preprocessed using function `vol_preprocess()` to estimate and re-scale co-occurence matrix `Y = X*t(X)` and its singular value decomposition:
  
  
  ```r
  vol <- vol_preprocess(t(simnmf$X), pfactor = 1e+5)
  ```

A parameter `pfactor` re-scales co-occurence matrix, but does not affect outcome. Currently, due to technical issues in implementation, too small `pfactor` can collapse volume to zero and we recommend to increase the factor to avoid numerical issues.

Next, `AnchorFree` function decomposes co-occurence matrix into matrices `C` and `E` of rank `n.comp` by minimizing volume of `E`:
  
  
  ```r
  vol.anchor <- AnchorFree(vol, n.comp = 10)
  
  str(vol.anchor$C)
  ```
  
  ```
  ##  num [1:192, 1:10] 0.003008 0.002325 0.000633 0.010146 0.008196 ...
  ```
  
  ```r
  str(vol.anchor$E)
  ```
  
  ```
  ##  num [1:10, 1:10] 304 129 138 132 135 ...
  ```

Comparison of original matrix `simnmf$C` and inferred matrix `volres$C` rescaled back by `vol$col.factors` shows that they have highly similar pairs of column vectors:
  
  
  ```r
  C <- vol.anchor$C*vol$col.factors
  apply(cor(simnmf$C, C), 1, max)
  ```
  
  ```
  ##  [1] 0.9777554 0.9750203 0.9792922 0.9779889 0.9802875 0.9818999 0.9833044 0.9738819 0.9752943
  ## [10] 0.9783594
  ```

Having inferred `C`, the second matrix `D` can be obtained throught linear regression with constraints or using the function `infer_intensities()`:
  
  
  ```r
  D <- infer_intensities(C, simnmf$X)
  ```

