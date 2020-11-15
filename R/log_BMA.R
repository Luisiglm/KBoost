#' Function to perform Bayesian Model Averaging with log un-normalized posteriors.
#'
#' @param posterior_list a list with GxK numeric matrix with the log-posteriors of TF regulation models.
#' @param model_list a list with GxK logical matrix that indicates if a TF was part of a regulation model.
#' @param G the number of genes.
#' @param K the number of TFs.
#' @param ite the number of iterations.
#' @export
#
log_BMA = function(posterior_list, model_list, G,K,ite){
  # Pre-allocate memory for the GRN. this is the output.
  bma = matrix(0,G,K)
  # we are going to access the elements of a each list. Then concatenate them in a matrix.
  for (i in 1:length(posterior_list)){
    if (i ==1){
      post = posterior_list[[i]]
    } else {
      post = cbind(post,posterior_list[[i]])
    }

  }
  # cool beans.
  # We will do bma per gene.
  for (i in 1:G){
    # we will remove the -infinity values cause they can cause numeric problems. They're value is zero in natural space they don't affect the sum.
    idx_inf = is.infinite(post[i,])
    p = post[i,!(idx_inf)]
    # access the log posteriors and repeat them into a matrix.
    p = matrix(p, length(p), length(p))
    # take the differences.
    p = p - t(p)
    # a difference in log space is a division in natural space.
    p = exp(p)
    # now let's average!
    p = colSums(p)
    p = p^-1
    # Make a new variable p_ with zeros.
    p_ = matrix(0,dim(post)[2],1)
    p_[!idx_inf] = p
    # Great! now we have the averaged model posteriors. We need to add them according the Tfs that belong in each model.
    j_o = 1
    j_end = K
    rang_ = j_o:j_end
    # Loop over the iterations to perform the model averaging.
    for (t in 1:ite){
      # Access the model at iteration t.
      model_t = model_list[[t]]
      # Keep only the model for gene i.
      model_t = model_t[i,]
      # At each iteration we have evaluated each TF. So we will do another loop. Maybe not too efficient?
      for (j in 1:K){
        # copy model_t
        model_j = model_t
        # add the TF j.
        model_j[j] = TRUE
        # add the normalized model to the model_j
        bma[i,model_j] = bma[i,model_j]+p_[rang_[j]]
      }
      # Now we need to move the window j_o and j_end.
      j_o = j_end + 1
      j_end = j_end + K
      rang_ = j_o:j_end
    }
  }
  return(bma)


}
