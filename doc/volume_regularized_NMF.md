
# Volume-regularized NMF

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
  
  Generally, structured matrix factorization problem is formulated as a decomposition of original matrix <img src="https://render.githubusercontent.com/render/math?math=X"> in a product of unknown non-negative matrix <img src="https://render.githubusercontent.com/render/math?math=C"> and matrix <img src="https://render.githubusercontent.com/render/math?math=D"> of lower rank: 
  
<img src="https://render.githubusercontent.com/render/math?math=X = C \cdot D">

  In practice, if <img src="https://render.githubusercontent.com/render/math?math=X"> has significantly more rows than columns or if observations are noisy and represent statistical sampling, then reformulation of the problem in the "covariance domain" simplifies inference:
  
<img src="https://render.githubusercontent.com/render/math?math=Y = C \cdot E \cdot C^{T},">

where <img src="https://render.githubusercontent.com/render/math?math=Y = X \cdot X^{T}, E = D \cdot D^{T}.">
  
  Under a relatively mild assumption on spread of column (row) vectors of matrix <img src="https://render.githubusercontent.com/render/math?math=C"> that belong to column (row) unit simplex, matrix <img src="https://render.githubusercontent.com/render/math?math=E"> has a minimum volume across all possible factorizations. _Vrnmf_, as compared to _AnchorFree_, considers a noisy version of the problem. _Vrnmf_ seeks to find a factorization pair of matrices <img src="https://render.githubusercontent.com/render/math?math=(C,E)"> that balances the goodness of fit of <img src="https://render.githubusercontent.com/render/math?math=\|Y-C \cdot E \cdot C^{T}\|_{F}^{2}"> and the volume of the matrix <img src="https://render.githubusercontent.com/render/math?math=E">. The method uses alternating optimization of the following objective function:
  
<img src="https://render.githubusercontent.com/render/math?math=\|Y-C \cdot E \cdot C^{T}\|_{F}^{2} %2B \lambda \cdot Volume(E).">

## Application in "covariance" domain

An application of the method is exemplified below for a simulated dataset. For that, we will use a dataset with matrices `X`, its noisy version `X.noise` and decomposition matrices `C` and `D` such that `X = C*D` (to simulate a dataset with parameters of interest one can use function `sim_factors`). First, we load the library and the dataset:
  
  
  ```r
  library(vrnmf)
  
  simnmf <- simulatedNMF::simnmf
  
  names(simnmf)
  #> [1] "X.noise" "X"       "C"       "D"
  ```

The summary below shows that the original matrix `X` of the size <img src="https://render.githubusercontent.com/render/math?math=192 \cdot 3000"> is a product of two matrices of lower rank 10:
  
  
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
  vol <- vol_preprocess(t(simnmf$X))
  ```

Next, the `vrnmf` function decomposes co-occurrence matrix into matrices `C` and `E` of rank `n.comp = 10` and using the weight of volume `wvol = 2e-2` (and additional options of intra-block Nesterov extrapolation `extrapolate` and inter-block Nesterov-style `accelerate`):
  
  
  ```r
  volres <- volnmf_main(vol, n.comp = 10, wvol = 2e-2, 
                      n.iter = 3e+3, vol.iter = 1e+1, c.iter = 1e+1, 
                      extrapolate = TRUE, accelerate = TRUE, verbose = FALSE)
  #> run standard nmf..
  #> done
  #> 
  #> run volume-regularized nmf..
  #> done
  #> 
  ```

The output `volres` contains the predicted matrix `C` and other auxiliary information:

```r
names(volres)
#>  [1] "C"      "R"      "Q"      "C.init" "R.init" "Q.init" "C.rand" "R.rand"
#>  [9] "Q.rand" "rec"

str(volres$C)
#>  num [1:192, 1:10] 0 0.00752 0.00329 0.00629 0.00776 ...
```

Comparison of original matrix `simnmf$C` and inferred matrix `volres$C` rescaled back by `vol$col.factors` shows that they have highly similar pairs of column vectors:
  
  
  ```r
  C <- volres$C*vol$col.factors
  apply(cor(simnmf$C, C), 1, max)
  #>  [1] 0.9163232 0.9032780 0.9729961 0.9935809 0.9873883 0.9933038 0.9944044
  #>  [8] 0.2875271 0.9791225 0.9859103
  ```

Having inferred `C`, the second matrix `D` can be obtained either through linear regression with constraints or using the function `infer_intensities()`:
  
  
  ```r
  D <- infer_intensities(C, simnmf$X)
  ```

## Structural and algorithmic parameters

The method has a number of variations in the type of structural constraints on `C` or `R`:
* The `C.constraint` parameter constrains row vectors (`row`) or column vectors (`col`) of `C` to unit simplex (sum of vector elements equal to 1). The former option is able to fit a general matrix factorization with non-negative matrix `C`, but the latter imposes a specific structural constraint usually met in hyperspectral unmixing.
* `R.constraint` parameter constrains matrix `R` to be non-negative (`pos` otherwise `no`). 

The objective function is nonconvex and thus alternating optimization does not guarantee to find a globally optimal solution. Additional parameters that affect convergence are:

* Parameter `wvol`, which regulates weight of volume.  Small values make the problem similar to standard NMF, while large values start collapsing eigenvalues of matrix `E` to zero. A carefully chosen `wvol` parameter is critical for good performance.
* Initialization matrices: Initial matrices `C`, `R` and `Q` can be obtained using NMF optimization without volume term  (default option, parameters `do.nmf` and `iter.nmf`) or initialized manually using `C.init`, `R.init` and `Q.init`. 
* Number of iterations. Beyond the overall number `n.iter` of iterations that alternate optimization of `C`, `R` and `Q`, the parameters `vol.iter` and `c.iter` regulate the number of iterations to update `R` and `C` at each alternating step respectively. These parameters control the convergence rate. Additionally, `R.majorate` controls the precision of the approximation of `R` at each iteration at the expense of time.

