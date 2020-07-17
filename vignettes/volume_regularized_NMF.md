# Volume_regularized_NMF

## Introduction  
  
  Generally, structured matrix factorization problem is formulated as a decomposition of original matrix <img src="https://render.githubusercontent.com/render/math?math=X"> in a product of unknown non-negative matrix <img src="https://render.githubusercontent.com/render/math?math=C"> and matrix <img src="https://render.githubusercontent.com/render/math?math=D"> of lower rank: 
  
<img src="https://render.githubusercontent.com/render/math?math=X = C \cdot D">

  In practice, if <img src="https://render.githubusercontent.com/render/math?math=X"> has significantly more rows than columns or if observations are noisy and represent statistical sampling than reformulation of the problem in "covariance domain" simplifies inference:
  
<img src="https://render.githubusercontent.com/render/math?math=Y = C \cdot E \cdot C^{T},">

where <img src="https://render.githubusercontent.com/render/math?math=Y = X \cdot X^{T}, E = D \cdot D^{T}.">
  
  Under a relatively mild assumption on spread of column (row) vectors of matrix <img src="https://render.githubusercontent.com/render/math?math=C"> that belong to column (row) unit simplex, matrix <img src="https://render.githubusercontent.com/render/math?math=E"> has minimum volume across all possible factorizations. _Vrnmf_, as compared to _AnchorFree_, considers noisy version of the problem. _Vrnmf_ seeks to find a factorization pair of matrices <img src="https://render.githubusercontent.com/render/math?math=(C,E)"> that balances goodness of fit of <img src="https://render.githubusercontent.com/render/math?math=\|Y-C \cdot E \cdot C^{T}\|_{F}^{2}"> and volume fo matrix <img src="https://render.githubusercontent.com/render/math?math=E">. The method uses alternating optimization of the following objective function:
  
<img src="https://render.githubusercontent.com/render/math?math=\|Y-C \cdot E \cdot C^{T}\|_{F}^{2} %2B \lambda \cdot Volume(E).">

## Application in "covariance" domain

Application of the method is examplified below for a simulated dataset. For that, we will use dataset with matrices `X`, its noisy version `X.noise` and decomposition matrices `C` and `D` such that `X = C*D` (to simulate a dataset with parameters of interest one can use function `sim_factors`). First, we upload the library and the dataset:
  
  
  ```r
  library(vrnmf)
  
  data(simnmf)
  
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

As a first step, input matrix `X` is preprocessed using function `vol_preprocess()` to estimate and re-scale co-occurence matrix `Y = X*t(X)` and its singular value decomposition:
  
  
  ```r
  vol <- vol_preprocess(t(simnmf$X))
  ```

Next, `vrnmf` function decomposes co-occurence matrix into matrices `C` and `E` of rank `n.comp = 10` and using weight of volume `wvol = 2e-2` (and additional options of intra-block Nesterov extrapolation `extrapolate` and inter-block Nesterov-style `accelerate`):
  
  
  ```r
  volres <- volnmf_main(vol, n.comp = 10, wvol = 2e-2, 
                      n.iter = 3e+3, vol.iter = 1e+1, c.iter = 1e+1, 
                      extrapolate = TRUE, accelerate = TRUE, verbose = FALSE)
  #> run standard nmf.. done
  #> run volume-regularized nmf.. done
  ```

The output `volres` contains predicted matrix `C` and other auxiliary information:

```r
names(volres)
#>  [1] "C"      "R"      "Q"      "C.init" "R.init" "Q.init" "C.rand" "R.rand" "Q.rand" "rec"

str(volres$C)
#>  num [1:192, 1:10] 0.012967 0 0.000464 0.002328 0 ...
```

Comparison of original matrix `simnmf$C` and inferred matrix `volres$C` rescaled back by `vol$col.factors` shows that they have highly similar pairs of column vectors:
  
  
  ```r
  C <- volres$C*vol$col.factors
  apply(cor(simnmf$C, C), 1, max)
  #>  [1] 0.9638703 0.9760804 0.9851034 0.9940353 0.9865532 0.9939367 0.9970672 0.9696916 0.9928067
  #> [10] 0.9885614
  ```

Having inferred `C`, the second matrix `D` can be obtained throught linear regression with constraints or using the function `infer_intensities()`:
  
  
  ```r
  D <- infer_intensities(C, simnmf$X)
  ```

## Structural and algorithmic parameters

The method has a number of variations in the type of structural constraints on `C` or `R`:
* `C.constraint` parameter constraints row vectors (`row`) or column vectors (`col`) of `C` to unit simplex (sum of vector elements equal to 1). The former option is able to fit a general matrix factorization with non-negative matrix `C`, but the latter imposes a specific structural constraint usually meet in hyperspectral unmixing.
* `R.constraint` parameter constraints matrix `R` to be non-negative (`pos` otherwise `no`). 

The objective function is nonconvex and thus alternating optimization does not guarantee to find globally optimal solution. Additional parameters that affect convergence are below:

* Parameter `wvol` regulates weight of volume.  Small values makes the problem similar to standard NMF, while large values start collapsing eigenvalues of matrix `E` to zero. Carefully choosen `wvol` is critical for good performance.
* Initialization matrices. Initial matrices `C`, `R` and `Q` can be obtained using NMF optimization without volume term  (default option, parameters `do.nmf` and `iter.nmf`) or initialized manually using `C.init`, `R.init` and `Q.init`. 
* Number of iterations. Beyond overall number `n.iter` of iterations that alternate optimization of `C`, `R` and `Q`, parameters `vol.iter` and `c.iter` regulate number of iterations to update `R` and `C` at each alternating step respectively. These parameters control convergence rate. Additionally, `R.majorate` control precision of approximation of `R` at each iteration at the expense of time.

