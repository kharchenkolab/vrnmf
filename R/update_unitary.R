# find Q: A~B*Q (Procrustes problem)
volnmf.procrustes <- function(A,B){
  sv <- svd( t(A)%*%B )
  sv$v%*%t(sv$u)
}
