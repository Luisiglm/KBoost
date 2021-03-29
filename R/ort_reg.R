#' Function to calculate a regression when regressors are orthogonal and have norm 1.
#'
#' @param X an NxT numeric matrix with the orthogonal normalized regressors (e.g. KPCs)
#' @param Y an NxG numeric matrix with the expression values of G genes and N obersvations.
#' @param v a number between 0 and 1 it corresponds to the shrinkage parameter.
#' @export
#' @return coefficients and log marginal likelihoods.
ort_reg = function(X,Y,v){
  # Y = X%*%b .
  # here t(X)%*%X = diag(T) (T is the number of columns of X). % denotes matrix multiplication in R.
  # b = solve(t(X)%*%X)%*%t(X)%*%Y = diag(T)%*%t(X)%*%Y.
  b = diag(dim(X)[2])%*%(t(X)%*%Y)
  f = v*X%*%b
  # Calculate the log marginal likelihood.
  llik = log_lik(f,Y)
  # Maybe we won't return f as it might be using too much memory?## CHANGE HERE ## We STOPPED RETURNING F.
  results = list(b = b, llik = llik)
  return(results)
}
