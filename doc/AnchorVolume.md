
# Anchorfree Algorithm

## Preliminary: Loading Simulated Data

We have pre-generated a dataset of simulated NMF matrices that you can load as a list directly using the package `simulatedNMF`, which is available through a `drat` repository on GitHub. Note that the size of the 'simulatedNMF' package is approximately 8 MB. This package may be installed as follows:


```r
install.packages('simulatedNMF', repos='https://kharchenkolab.github.io/drat/', type='source')
```

(Please see the [drat documentation](https://dirk.eddelbuettel.com/code/drat.html) for more comprehensive explanations and vignettes regarding `drat` repositories.)


The following command load this dataset for the following tutorial:

```r
simnmf <- simulatedNMF::simnmf
```


## Introduction  
  
  Generally, the structured matrix factorization problem is formulated as a decomposition of original matrix <img src="https://render.githubusercontent.com/render/math?math=X"> in a product of unknown non-negative matrix <img src="https://render.githubusercontent.com/render/math?math=C"> and matrix <img src="https://render.githubusercontent.com/render/math?math=D"> of lower rank: 
  
    
<img src="https://render.githubusercontent.com/render/math?math=\begin{aligned}X = C \cdot D \end{aligned}">

  In practice, if <img src="https://render.githubusercontent.com/render/math?math=X"> has significantly more rows than columns or if observations are noisy and represent statistical sampling, then reformulation of the problem in the "covariance domain" simplifies inference:
  
<img src="https://render.githubusercontent.com/render/math?math=\begin{aligned}Y = C \cdot E \cdot C^{T},\end{aligned}">

  where <img src="https://render.githubusercontent.com/render/math?math=Y = X \cdot X^{T}, E = D \cdot D^{T}.">
  
  The goal of the _AnchorFree_ approach is to tri-factorize matrix <img src="https://render.githubusercontent.com/render/math?math=Y"> in a product of matrices <img src="https://render.githubusercontent.com/render/math?math=C"> and <img src="https://render.githubusercontent.com/render/math?math=E">. Under a relatively mild assumption on the spread of column vectors of matrix <img src="https://render.githubusercontent.com/render/math?math=C">, matrix <img src="https://render.githubusercontent.com/render/math?math=E"> has the minimum volume across all possible factorizations. _AnchorFree_ seeks to find a factorization pair of matrices <img src="https://render.githubusercontent.com/render/math?math=(C, E)"> such that the column vectors of <img src="https://render.githubusercontent.com/render/math?math=C"> represent the unit simplex and the matrix <img src="https://render.githubusercontent.com/render/math?math=E"> has the minimum determinant (as a proxy of volume) using alternating linear programming. 

## Application

Application of the method is exemplified below for a simulated dataset. For that, we will use dataset with matrices `X`, its noisy version `X.noise` and decomposition matrices `C` and `D` such that `X = C*D` (to simulate a dataset with parameters of interest one can use function `sim_factors`). First, we load the library and dataset:
  
  
  ```r
  library(vrnmf)
  
  simnmf <- simulatedNMF::simnmf
  
  names(simnmf)
  #> [1] "X.noise" "X"       "C"       "D"
  ```

Summary below shows that the original matrix `X` of the size <img src="https://render.githubusercontent.com/render/math?math=192 \cdot 3000"> is a product of two matrices of lower rank 10:
  
  
  ```r
  str(simnmf$X)
  #>  num [1:192, 1:3000] 0.000493 0.000695 0.000347 0.000248 0.000699 ...
  
  str(simnmf$C)
  #>  num [1:192, 1:10] 0.00387 0.01099 0.00302 0 0.008 ...
  
  str(simnmf$D)
  #>  num [1:10, 1:3000] 0.03793 0.001 0.00138 0.00229 0.001 ...
  ```

As a first step, the input matrix `X` is preprocessed using the function `vol_preprocess()` to estimate and rescale the co-occurrence matrix `Y = X*t(X)` and its singular value decomposition:
  
  
  ```r
  vol <- vol_preprocess(t(simnmf$X), pfactor = 1e+5)
  ```

A parameter `pfactor` re-scales the co-occurrence matrix, but it does not affect the outcome. Currently, due to technical issues in the current implementation, too small `pfactor` can collapse volume to zero and we recommend increasing the factor to avoid numerical issues.

Next, the `AnchorFree` function decomposes the  co-occurrence matrix into matrices `C` and `E` of rank `n.comp` by minimizing the volume of `E`:
  
  
  ```r
  vol.anchor <- AnchorFree(vol, n.comp = 10)
  
  str(vol.anchor$C)
  #>  num [1:192, 1:10] 0.003008 0.002325 0.000633 0.010146 0.008196 ...
  
  str(vol.anchor$E)
  #>  num [1:10, 1:10] 304 129 138 132 135 ...
  ```

Comparison of the original matrix `simnmf$C` and the inferred matrix `vol.anchor$C` shows that they have highly similar pairs of column vectors:
  
  
  ```r
  C <- vol.anchor$C*vol$col.factors
  apply(cor(simnmf$C, C), 1, max)
  #>  [1] 0.9999967 1.0000000 1.0000000 0.9993503 1.0000000 1.0000000 1.0000000
  #>  [8] 0.9999061 1.0000000 0.9999527
  ```

Having inferred `C`, the second matrix `D` can be obtained either through linear regression with constraints or using the function `infer_intensities()`:
  
  
  ```r
  D <- infer_intensities(C, simnmf$X)
  ```

