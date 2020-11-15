#'Function to automatically build and run KBoost on data with human gene expression values and Symbol gene names.
#' @param X an NxG numeric matrix with the expression values of G genes and N obersvations. The gen_names can come as column names.
#' @param gen_names a set of SYMBOL gene names that correspond to the names of the columns of X.
#' @param g a positive scalar with the width parameter for the RBF kernel. (default g = 40)
#' @param v a double between 0 and 1 with the shrinkage parameter. (default v =0.1)
#' @param ite an integer with the number of iterations (default ite = 3)
#' @param pos_weight a scalar between 0 and 1. This is the prior probability of observing a TF regulate a gene given that this interaction was observed before.
#' @param neg_weight a scalar between 0 and 1. This is the prior probability of observing a TF regulate a gene given that this interaction was NOT observed before.
#' @export
KBoost_human_symbol = function(X,gen_names,g, v,ite,pos_weight,neg_weight){
  # check the input making sure that:
  # 1) the input has gene names.
  if (missing(gen_names) && is.null(colnames(X))){
    stop("No gene names were found. Gene names can come either as a vector as an argument or as column names of X.")
  } else if (missing(gen_names)){
    gen_names = colnames(X)
  } else if (length(gen_names)!=dim(X)[2]){
        if (length(gen_names)>dim(X)[2]){
          stop("There are more gene names than columns of X. Please double check the inputs.")
        } else {
          stop("There are less gene names than columns of X. Please double check the inputs.")
        }
  }
  # 2) there are Tfs in the input. # The function tests this automatically.
  TFs = get_tfs_human(gen_names)
  # If it all passed. It's time to build a priors and run KBoost if no errors were flagged.
  # If the user did not submit any negative or positive weights. We can add the default value.
  # We will check that pos_Weight and neg_weight should be between 0 and 1.
  if (missing(pos_weight)){
    pos_weight = 0.6
  } else if (pos_weight>=1 || pos_weight<=0){
    stop("<pos_weight> needs to be lower than 1 and higher than 0.")
  }
  if (missing(neg_weight)){
    neg_weight = 0.5
  } else if (neg_weight>=1 || neg_weight<=0){
    stop("<neg_weight> needs to be lower than 1 and higher than 0.")
  }
  # Make sure that neg_weight is lower than pos_weight.
  if (neg_weight>pos_weight){
    stop("<neg_weight> needs to be lower than <pos_weight>.")
  }
  prior_weights = get_prior_Gerstein(gen_names,TFs,pos_weight,neg_weight)
  grn = kboost(X,TFs,g,v,prior_weights,ite)
  return(grn)
}
