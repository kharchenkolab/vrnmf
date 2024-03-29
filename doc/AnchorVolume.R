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
vol <- vol_preprocess(t(simnmf$X), pfactor = 1e+5)

## -----------------------------------------------------------------------------
vol.anchor <- AnchorFree(vol, n.comp = 10)

str(vol.anchor$C)

str(vol.anchor$E)

## -----------------------------------------------------------------------------
C <- vol.anchor$C*vol$col.factors
apply(cor(simnmf$C, C), 1, max)

## -----------------------------------------------------------------------------
D <- infer_intensities(C, simnmf$X)

