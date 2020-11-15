#' Function to calculate the RBF Kernel of a matrix X with width g.
#'
#' @param x an Nx1 numeric matrix with N observations.
#' @param g a positive scalar with the width parameter.
#' @export
#'
RBF_K = function(x,g){

  # formula K(i,j) = exp(-Euc[i,j]/g), Euc is the Squared Euclidean distance.
  # Calculate the Squared Euclidian Distance between Every Observation.
  # Euc[i,j] = (x[i]-x[j])^2 = x[i]^2 - 2x[i]*x[j]-x[j]^2.
  # In matrix operations it yields the following:

  x2 = x^2

  X = matrix(x,length(x),length(x))

  X2 = matrix(x2, length(x),length(x))

  # Take the exponential and divide it by g.

  K = exp(-(X2-2*X*t(X)+t(X2))/g)

  return(K)
}
