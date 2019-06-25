Katz_cent <- function(network){
  n = nrow(network)
  I = diag(n)
  K = I - 0.01*t(network)
  inv_K = solve(K)
  In = matrix(1,nrow(network),1)
  Katz_cent = inv_K%*%In
  names(Katz_cent) = rownames(network)
  return(Katz_cent)
}
