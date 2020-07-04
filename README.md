# vrnmf
volume-regularized NMF.

Package implements a set of methods to perform non-negative matrix decomposition with minimum volume contraints. 

A general problem is to decompose a non-negative matrix $X_{nm}$ in a product of non-negative matrix $C_{nr}$ and matrix $D_{rm}$ of lower rank $r$: $X = C\timesD$. In case of non-negativity constraint on matrix $D$ the problem is known as NMF.

(how the problem is called with only one non-negative matrix?)
3) However, the problem, and NMF as a particualr case, is not identifiable in general case (cite). There are potentially infite number of solutions. It makes interpretation of obtained results challenging and limit NMF applicability to dimensionality reduction.

4) Recent theoretical advances showed that the issue can be overcome under a relatively mild assumption on "spread" of column vectors of C known as sufficient spreading (CHECK IF THERE IS A N ASSUMPTION ON THE SECOND MATRIX, cites). If matrix C is non-negative and have sufficiently spread column vectors than volume minimization of a matrix D delivers a correct and unique solution (C, D). In that way,

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