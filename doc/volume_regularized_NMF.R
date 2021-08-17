## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(vrnmf)

simnmf <- simulatedNMF::simnmf

names(simnmf)

## -----------------------------------------------------------------------------
str(simnmf$X)

str(simnmf$C)

str(simnmf$D)

## -----------------------------------------------------------------------------
vol <- vol_preprocess(t(simnmf$X))

## -----------------------------------------------------------------------------
volres <- volnmf_main(vol, n.comp = 10, wvol = 2e-2, 
                    n.iter = 3e+3, vol.iter = 1e+1, c.iter = 1e+1, 
                    extrapolate = TRUE, accelerate = TRUE, verbose = FALSE)

## -----------------------------------------------------------------------------
names(volres)

str(volres$C)

## -----------------------------------------------------------------------------
C <- volres$C*vol$col.factors
apply(cor(simnmf$C, C), 1, max)


## -----------------------------------------------------------------------------
D <- infer_intensities(C, simnmf$X)

