# vrnmf
volume-regularized NMF.

Package implements a set of methods to perform non-negative matrix decomposition with minimum volume contraints. 

A general problem is to decompose a non-negative matrix <img src="https://render.githubusercontent.com/render/math?math=X_{nm}"> in a product of non-negative matrix <img src="https://render.githubusercontent.com/render/math?math=C_{nr}"> and matrix <img src="https://render.githubusercontent.com/render/math?math=D_{rm}"> of lower rank r: <img src="https://render.githubusercontent.com/render/math?math=X = C\cdot D">. In case of additional non-negativity constraint on matrix <img src="https://render.githubusercontent.com/render/math?math=D"> the problem is known as NMF.

This problem, and NMF as a particualr case, is not identifiable in general case, meaning that there are potentially many different solutions that deliver the same decomposition quality. It makes interpretation of factorized matrices challenging and limits applications of NMF to instrumental dimensionality reduction. However, recent theoretical advances showed that the issue can be overcome under a relatively mild assumption on "spread" of column vectors of C known as sufficient spreading (SECOND MATRIX? CITES!). If matrix C is non-negative and have sufficiently spread column vectors than volume minimization of a matrix D delivers a correct and unique, up to a scale and permutation, solution (C, D). 

5) An approach AnhorFree enables efficient estimate of matrix C by reformulating the problem in covariance domain (cite). A short walkthrough is here.

6) A more general formulation of the problem that accounts for noise in matrix X, such that only approximately X ~ CD, is called volume-regularized NMF (vrnmf). Vrnmf balances goodness of matrix approximation and size volume (cites):
||X-CD|| + lambda*volume(D).

7) Here we provide implementation of vrnmf and reformulation of vrnmf in covariance domain. A detailed walkthrough is available here.

TODO:
- formulas.
- links.
- citations.
- ADD ABOUT GEOMETRIC INTERPRETATION - SOLVE EXACT NMF IS THE SAME AS FIND CONE OF POINTS.

## Installation instructions

```{r setup}
install.packages("devtools")
devtools::install_github("kharchenkolab/vrnmf")
library(vrnmf)
```
