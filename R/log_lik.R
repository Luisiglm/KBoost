#' Function to calculate the Log-Marginal Likelihood of a Function of the Kernel PC regression model.
#'
#' @param Y an NxG numeric matrix with the expression values of G genes and N observations.
#' @param F an NxG numeric matrix with the predicted expression values of G genes and N observations.
#' @export
#'
log_lik= function(Y,F){

  # Calculate the sum of squared errors.

  llik = colSums((F-Y)^2)/dim(Y)[1]

  # Elevate to the power of (n/2). n is the number of observations.
  llik = log(llik)*(-(dim(Y)[1])/2)

  return(llik)
}
