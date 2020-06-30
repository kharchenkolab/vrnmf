
# preprocess data
vol.preprocess <- function(X,col.norm="sd",row.norm=NULL){
  row.factors <- rep(1,nrow(X))
  col.factors <- rep(1,ncol(X))

  if (col.norm=="sd"){
    col.factors <- apply(X,2,sd)
  }
  X.process <- t(t(X)/col.factors)/row.factors

  P0 <- t(X.process)%*%X.process
  pfactor <- max(P0)
  P <- P0/pfactor # how to normalize P to avoid inf problems?
  dimr <- svd(P)

  return( list(X.process = X.process, row.factor=row.factors, col.factors=col.factors,
               P0=P0, P = P, pfactor = pfactor, U = dimr$u, eigens = dimr$d) )
}
