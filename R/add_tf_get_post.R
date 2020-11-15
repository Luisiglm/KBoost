#' Function to add a new tf to an existing model and calculate posteriors.
#'
#' @param priors a GxK numeric matrix with the log-priors of TF regulation models.
#' @param lliks a GxK numeric matrix with the log-marginal likelihoods of TF regulation models.
#' @param model a GxK binary matrix with the models at the previous iteration.
#' @param f an NxG matrix with the predicted gene expression values.
#' @param reg a list produced from the ort_reg function.
#' @param kpca a list with the Kernel Principal Components.
#' @param TFs indexes of the genes in the system integers that are TFs.
#' @param v the shrinkage parameter.
#' @export
#'
add_tf_get_post = function(priors,lliks,model,f,reg, kpca, TFs,v){
  # For each gene will select the new TF with the highest posterior.
  # the number of genes are the rows.
  G = dim(priors)[1]
  # the number of TFs are the columns.
  K = dim(priors)[2]
  posteriors = lliks + priors
  # Store the priors of the best models.
  prior_best = matrix(0,G,1)
  # For each gene in G do.
  for (i in 1:G){
    # Find the Tf with the highest posterior.
    best = which.max(posteriors[i,])
    # update model.
    model[i,best] = TRUE
    # update prior
    prior_best[i] = priors[i,best]
    # update f.
    # adjust best for potential indexes issues related to tf_check.
    if (TFs[best]<i){
      idx = i - 1
    } else {
      idx = i
    }
    #### USE COEFFICIENTS INSTEAD### We also added v to the input.
    f[,i] = f[,i]+v*kpca[[best]]%*%reg[[best]]$b[,idx]
  }
  return(list(model = model,posteriors = posteriors, prior_best = prior_best, f = f))
}

