#' Function to perform boosting iterations in KBoost.
#'
#' @param X an NxG matrix with the expression values.
#' @param f an NxG matrix with the predicted gene expression values.
#' @param kpca a lits with the Kernel Principal Components.
#' @param v a double between 0 and 1 that corresponds to the shrinkage parameter.
#' @param TFs a matrix with integers with the columns of X that are TFs.
#' @export
#' @return list with log likelihoods and residuals
greedy_uni_boosting = function(f, X, TFs, kpca,v){
  # Calculate pseudo-residuals.
  pse = X-f
  # Pre-allocate memory for outputs.
  res = list()
  llik = matrix(-Inf,dim(X)[2],length(TFs))
  # Cool Beans! Now, we'll do a for loop along the TFs.
  for (j in seq_len(length(TFs))){
    # do an orthogonal regression with the KPCAs and the pseudo-residuals.
    res[[j]] = ort_reg(kpca[[j]],pse[,-TFs[j]],v)
    # get log_likelihoods.
    llik[-TFs[j],j] = res[[j]]$llik

  }
  return(list(res = res, llik = llik))
}
