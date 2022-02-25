[![<kharchenkolab>](https://circleci.com/gh/kharchenkolab/vrnmf.svg?style=svg)](https://app.circleci.com/pipelines/github/kharchenkolab/vrnmf)
[![CRAN status](https://www.r-pkg.org/badges/version/vrnmf)](https://cran.r-project.org/package=vrnmf)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/vrnmf)](https://cran.r-project.org/package=vrnmf)
  
# vrnmf: Volume-regularized NMF

The R package `vrnmf` implements a set of methods to perform non-negative matrix decomposition with minimum volume constraints. A general problem is to decompose a non-negative matrix <img src="https://render.githubusercontent.com/render/math?math=X_{nm}"> in a product of non-negative matrix <img src="https://render.githubusercontent.com/render/math?math=C_{nr}"> and matrix <img src="https://render.githubusercontent.com/render/math?math=D_{rm}"> of lower rank r: <img src="https://render.githubusercontent.com/render/math?math=X = C\cdot D">. In case of additional non-negativity constraints on the matrix <img src="https://render.githubusercontent.com/render/math?math=D">, the problem is known as NMF. 

This problem, and NMF as a particular case, is not identifiable in the general case, meaning that there are potentially many different solutions that deliver the same decomposition quality [[1]](#1). This both makes interpretation of factorized matrices challenging and limits applications of NMF to instrumental dimensionality reduction. However, recent theoretical advances have shown that the issue can be overcome under a relatively mild assumption based on "spread". That is, the column vectors of C are known as "sufficiently spread"[[2-3]](#2) if the matrix C is non-negative and the matrix C has sufficiently spread column vectors then the volume minimization of a matrix D delivers a correct and unique, up to a scale and permutation, solution (C, D). 

The _AnchorFree_ approach enables efficient estimation of matrix C by reformulating the problem in the covariance domain using the application of the volume minimization criterion [[4]](#4). A short walkthrough can be found [here](https://github.com/kharchenkolab/vrnmf/blob/main/doc/volume_regularized_NMF.md)

A more general formulation of the problem that accounts for noise in matrix X, such that only approximately <img src="https://render.githubusercontent.com/render/math?math=X \approx CD">, is called volume-regularized NMF (_vrnmf_). To balance goodness of matrix approximation and matrix D volume, _vrnmf_ minimizes the following objective function [[5-6]](#2):

<img src="https://render.githubusercontent.com/render/math?math=F = \| X-CD \|_{F}^{2} %2B \lambda \cdot Vol(D)"> 

We provide implementation of `vrnmf` approach and devise its reformulation in covariance domain.




## Walkthrough

### Volume-regularized NMF:

* [Markdown version](https://github.com/kharchenkolab/vrnmf/blob/main/doc/AnchorVolume.md)

### Anchorfree Algorithm:

* [Markdown version](https://github.com/kharchenkolab/vrnmf/blob/main/doc/volume_regularized_NMF.md)


## Installation 


To install the stable version from CRAN, use:

```R
install.packages('vrnmf')
```


To install the latest version of vrnmf, use:

```R
install.packages('devtools')
devtools::install_github('kharchenkolab/vrnmf')
```


## References

#### R package

The R package can be cited as:

```
Ruslan Soldatov, Peter Kharchenko, Viktor Petukhov, and Evan Biederstedt (2021). vrnmf:
Volume-regularized structured matrix factorization. R package version
1.0.2. https://github.com/kharchenkolab/vrnmf
```

#### Publication

If you find this software useful for your research, please cite the corresponding [paper](https://www.science.org/doi/10.1126/science.aba7408):

```
Vladimir B. Seplyarskiy Ruslan A. Soldatov, et al. 
Population sequencing data reveal a compendium of mutational processes in the human germ line.
Science, 12 Aug 2021. doi: 10.1126/science.aba7408
```

#### README references

<a id="1">[1]</a> 
K. Huang, N. D. Sidiropoulos and A. Swami.
Non-Negative Matrix Factorization Revisited: Uniqueness and Algorithm for Symmetric Decomposition.
IEEE Transactions on Signal Processing, vol. 62, no. 1, pp. 211-224, Jan.1, 2014, doi: 10.1109/TSP.2013.2285514.

<a id="2">[2]</a> 
C.-H. Lin, W.-K. Ma, W.-C. Li, C.-Y. Chi, and A. Ambikapathi.
Identifiability of the simplex volume minimization criterion for blind hyperspectral unmixing: The no-pure-pixel case.
IEEE Trans. Geosci.
Remote Sens.,  vol. 53, no. 10, pp. 5530–5546, Oct 2015.

<a id="3">[3]</a> 
X. Fu, W.-K. Ma, K. Huang, and N. D. Sidiropoulos.
Blind separation of quasi-stationary sources: Exploiting convex geometry in covariance domain.
IEEE Trans. Signal Process., vol. 63, no. 9, pp. 2306–2320, May 2015.

<a id="4">[4]</a> 
X. Fu, K. Huang, N. D. Sidiropoulos, Q. Shi, and M. Hong.
Anchorfree correlated topic modeling.
IEEE Trans. Patt. Anal. Machine Intell., vol. to appear, DOI: 10.1109/TPAMI.2018.2827377, 2018.

<a id="5">[5]</a> 
X. Fu, K. Huang, B. Yang, W. Ma and N. D. Sidiropoulos.
Robust Volume Minimization-Based Matrix Factorization for Remote Sensing and Document Clustering.
IEEE Transactions on Signal Processing, vol. 64, no. 23, pp. 6254-6268, 1 Dec.1, 2016, doi: 10.1109/TSP.2016.2602800.

<a id="6">[6]</a> 
A. M. S. Ang and N. Gillis.
Algorithms and Comparisons of Nonnegative Matrix Factorizations With Volume Regularization for Hyperspectral Unmixing.
IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, vol. 12, no. 12, pp. 4843-4853, Dec. 2019, doi: 10.1109/JSTARS.2019.2925098.
