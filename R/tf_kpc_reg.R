#' Function to perform a univariate RBF KPC regression per TF in X.
#'
#' @param X an NxG numeric matrix with the expression values of G genes and N obersvations.
#' @param TFs a Kx1 numeric matrix with integers that are the indexes of columns of X with genes that encode TFs.
#' @param g a positive scalar with the width parameter for the RBF kernel.
#' @param v a double between 0 and 1 with the shrinkage parameter.
#' @export
#'
tf_kpc_reg = function(X,TFs,g,v){
  # Create a list to store the values of the KPC.
  kpca = list()
  res = list()
  llik = matrix(-Inf,dim(X)[2],length(TFs))
  for (j in seq_along(TFs)){
    # Get the RBF Kernel.
    k = RBF_K(X[,TFs[j]],g)
    # Get the KPCs.
    kpca[[j]] = KPC(k, 0.01)
    # Do the orthogonal regression. Exclude the Gene for the same TF. Autoregulation is not included here.
    res[[j]] = ort_reg(kpca[[j]],X[,-TFs[j]],v)
    llik[-TFs[j],j] = res[[j]]$llik
  }
  return(list(results = res, kpca = kpca, llik = llik))
}
