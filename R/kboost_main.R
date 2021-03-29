#' Function to perform a univariate RBF KPC regression per TF in X.
#'
#' @param X an NxG numeric matrix with the expression values of G genes and N obersvations.
#' @param TFs a Kx1 numeric matrix with integers that are the indexes of columns of X with genes that encode TFs.
#' @param g a positive scalar with the width parameter for the RBF kernel. (default g = 40)
#' @param v a double between 0 and 1 with the shrinkage parameter. (default v =0.1)
#' @param prior_weights a GxK matrix with the prior weights for the prior probability that a TF regulates a gene.(default 0.5 for all values)
#' @param ite an integer with the number of iterations (default ite = 3)
#' @export
#' @return list with Kboost results.
#'
kboost_main = function(X,TFs,g,v,prior_weights,ite){
  # do first part of the algorithm.
  # initialize f. In pre-processing we have set the mean of X to zero.
  pc = proc.time()[3]
  N = dim(X)[1]
  G = dim(X)[2]
  K = length(TFs)
  f = matrix(0,N,G)
  # cool beans! Let's do the first part.
  # the function will run a for loop, calculate RBF_kernel per TF, the Kernel principal components, and the regression results.
  part_1 = tf_kpc_reg(X,TFs,g,v)
  # part_1 is a list which contains a) kpca a list with the K, kernel principal components and res, with the regression results.
  # We can put the regression results in a matrix for clarity.
  kpca = part_1$kpca
  lliks = part_1$llik
  ## Maybe clear part_1 after??
  # Initialize the prior and model.
  p_m = init_prior(prior_weights)
  prior = p_m$prior
  model = p_m$model
  ## clear p_m after??
  # We need to make tf_check to avoid indexing issues.
  # update priors for the potential models.
  prior_new = prior_model(prior_weights, model, prior)
  # add a tf per gene and get posteriors.### CHANGE NEW VERSION WE ADD V.
  added_tfs = add_tf_get_post(prior_new,lliks,model,f, part_1$results, kpca,TFs,v)
  # cool beans! We need to store the results. we keep the models and posteriors in two lists.
  model_list =  list()
  posterior_list = list()
  model_list[[1]] = model
  model = added_tfs$model
  # update f
  f = added_tfs$f
  # We will update the model after, as this is necessary for the next iterations.
  posterior_list[[1]] = added_tfs$posteriors
  # update the prior.
  prior = matrix(added_tfs$prior_best,G,K)
  # If we have more than 1 iterations, as is usually the case, we will run a for loop a perform boosting.
  if (ite>1){
    for (i in 2:ite){
      # We will use the function greedy_uni_boosting.
      res = greedy_uni_boosting(f, X,TFs, kpca, v)
      # update priors for the potential models.
      prior_new = prior_model(prior_weights, model, prior)
      # add a tf per gene and get posteriors. The same tf can be added multiple times. ## CHANGE NEW VERSION WE ADD V.
      added_tfs = add_tf_get_post(prior_new,res$llik,model,f, res$res, kpca,TFs,v)
      # update fs.
      f = added_tfs$f
      # store old model.
      model_list[[i]] = model
      # update model.
      model = added_tfs$model
      # update posteriors list.
      posterior_list[[i]] =  added_tfs$posteriors
      # update prior.
      prior = matrix(added_tfs$prior_best,G,K)
    }
  }
  # Cool beans! Now we need to combine the iterations. we have the log_BMA function that does it.
  BMA = log_BMA(posterior_list, model_list, G,K,ite)
  # that's what we are returning son! Thanks for using KBoost, we hope you enjoyed the ride!
  # we will do the heuristic post-processing.
  BMA_proc = net_refine(BMA)
  return(list(GRN= BMA_proc, GRN_UP = BMA,model = model, g = g, v = v, prior = prior,TFs = TFs, prior_weights = prior_weights, run_time = proc.time()[3] - pc ))


}
