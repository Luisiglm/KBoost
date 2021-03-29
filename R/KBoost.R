#' A function to run KBoost.
#'
#' @param X an NxG numeric matrix with the expression values of G genes and N obersvations.
#' @param TFs a Kx1 numeric matrix with integers that are the indexes of columns of X with genes that encode TFs.
#' @param g a positive scalar with the width parameter for the RBF kernel. (default g = 40)
#' @param v a double between 0 and 1 with the shrinkage parameter. (default v =0.1)
#' @param prior_weights it can be a scalar or a GxK matrix with the prior weights for the prior probability that a TF regulates a gene.(default 0.5 for all values)
#' @param ite an integer with the number of iterations (default ite = 3)
#' @export
#' @return results for kboots.
#' @examples
#' data(D4_multi_1)
#' Net = kboost(D4_multi_1)
#'
kboost = function(X,TFs,g,v,prior_weights,ite){
  # First we will check the input and make sure everything is in the right form.
  inpts = check_input(X, TFs,g,v,prior_weights,ite)
  # Now that we have the inputs. We can proceed.
  print('KBoost has checked your input and will proceed. Once we are finished, a random completion message will appear.')
  grn = kboost_main(inpts$X, inpts$TFs, inpts$g, inpts$v, inpts$prior_weights, inpts$ite)
  # Let's return grn.
  # format the results a bit if gene names were given.
  if (!is.null(colnames(X))){
    grn = add_names(grn,colnames(X))
  }

  print_completion_message()
  return(grn)
  # Print completion message!
  # Cool beans! Thanks for using KBoost! Have a great day!
}
