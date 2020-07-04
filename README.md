# vrnmf
Volume-regularized NMF.

Package implements a set of methods to perform non-negative matrix decomposition with minimum volume contraints. A general problem is to decompose a non-negative matrix <img src="https://render.githubusercontent.com/render/math?math=X_{nm}"> in a product of non-negative matrix <img src="https://render.githubusercontent.com/render/math?math=C_{nr}"> and matrix <img src="https://render.githubusercontent.com/render/math?math=D_{rm}"> of lower rank r: <img src="https://render.githubusercontent.com/render/math?math=X = C\cdot D">. In case of additional non-negativity constraint on matrix <img src="https://render.githubusercontent.com/render/math?math=D"> the problem is known as NMF. 

This problem, and NMF as a particualr case, is not identifiable in general case, meaning that there are potentially many different solutions that deliver the same decomposition quality [[1]](#1). It makes interpretation of factorized matrices challenging and limits applications of NMF to instrumental dimensionality reduction. However, recent theoretical advances showed that the issue can be overcome under a relatively mild assumption on "spread" of column vectors of C known as sufficient spreading [[2-3]](#2). If matrix C is non-negative and have sufficiently spread column vectors than volume minimization of a matrix D delivers a correct and unique, up to a scale and permutation, solution (C, D). 




_AnhorFree_ approach enables efficient estimation of matrix C by reformulating the problem in covariance domain following by application of volume minimization criterion (CITES!). A short walkthrough can be found hereLINK.

A more general formulation of the problem that accounts for noise in matrix X, such that only approximately <img src="https://render.githubusercontent.com/render/math?math=X \approx CD">, is called volume-regularized NMF (_vrnmf_). To balance goodness of matrix approximation and matrix D volume, _vrnmf_ minimizes the following objective function (cites):

<img src="https://render.githubusercontent.com/render/math?math=F = \| X-CD \|_{F}^{2} %2B \lambda \cdot Vol(D)"> 

We provide implementation of _vrnmf_ approach and devise its reformulation in covariance domain. A detailed walkthrough is available hereLINK.

"...the **go to** statement should be abolished..." [[1]](#1).

"...the **go to** statement should be abolished..." [[2]](#Diatal).

TODO:
- links.
- citations.
- ADD VISUAL GEOMETRIC INTERPRETATION.

## Installation instructions

```{r setup}
install.packages("devtools")
devtools::install_github("kharchenkolab/vrnmf")
library(vrnmf)
```

## References
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

