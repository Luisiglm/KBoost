#' Function to initialize the log prior of a model given our assumptions.
#'
#' @param prior_weights an GxK numeric matrix with the prior probability of a TF in our set to regulate a gene.
#' @export
#'
init_prior = function(prior_weights){
  # if Prior has 1 values it is better they are removed as it will yield sets of models with -Inf probabilities.
  # we will return p_o and Model as a standard all False binary matrix.
  G = dim(prior_weights)[1]
  K = dim(prior_weights)[2]
  prior = matrix(0,G,K)
  p = rowSums(log(1-prior_weights))
  for (j in 1:K){
    # Initialize a prior with all TFs missing.
    prior[,j] = p

  }
  model = matrix(FALSE, G,K)
  return(list(prior = prior, model= model))
}
