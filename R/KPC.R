#' Function to calculate the principal components of a kernel.
#'
#' @param K an NxN numeric matrix with the Kernel matrix.
#' @param thr a positive scalar which is a threshold to discard eigen-vectors based on eigen-values.
#' @export
#' @return the kernel principal components
#' @examples
#' x = rnorm(100,0,1)
#' k = RBF_K(x,1)
#' k_ = kernel_normal(k)
#' kpca = KPC(k,1e-8)
KPC = function(K,thr){
  # Perform the eigen decomposition.
  Ei = eigen(K)
  # Ei is a structure with eigen values and eigen vectors.
  Pass = (Re(Ei$values)>thr)
  # Project the eigen-vectors on the Kernel matrix.
  if (sum(Pass) > 0){
      KPC = K%*%Re(Ei$vectors[,Pass])
       # Set their euclidean norm to one.
      for (i in seq_len(dim(KPC)[2])){
      KPC[,i] = KPC[,i] / sum(KPC[,i]^2)^(1/2)
      }
      return(KPC)
  } else {
    stop("No Eigen-values lower than thr. You can try using a lower thr")
  }


}
