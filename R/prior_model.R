#' Function to calculate the log prior of a model given our assumptions.
#'
#' @param prior_weights an GxK numeric matrix with the prior probability of a TF in our set to regulate a gene.
#' @param model an GxK binary matrix with values of True and False if a TF was added to a model.
#' @param prior a Gx1 numeric matrix with the prior probability of the current model.
#' @export
#' @return prior probability of a model
prior_model = function(prior_weights, model, prior){
  # At the start of each boosting iteration we will calculate the new model priors.
  G = dim(prior_weights)[1]
  K = dim(prior_weights)[2]
  # We will return prior modified.
  # we will add the log prior weight if TF is not in current and delete the log(1-Prior[i,j]).
  # Otherwise it remains the same.
  for (j in 1:K){
    # For every gene we have a binary model with the TFs added in previous iterations.
    for (i in 1:G){
      # If TF was NOT there before.
      if (!model[i,j]){

        prior[i,j] = prior[i,j] + log(prior_weights[i,j]) - log(1-prior_weights[i,j])
      }

    }
  }
  # Great! Now, Return prior..
  return(prior)
}
