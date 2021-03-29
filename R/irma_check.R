#' Function to produce the AUPR and AUROC Results on the DREAM4 Multifactorial Challenge.
#' @param g a number larger than 0 that is the width parameter for the RBF Kernel
#' @param v a number between 0 and 1 that is the shrinkage parameter
#' @param ite an integer with number of iterations.
#' @export
#' @return list with aurocs and auprs for IRMA datasets
#' @examples
#' res = irma_check()
#'
irma_check = function(v,g,ite){
  # Pre-allocate memory for the results.
  aurocs = matrix(0,2,1)
  auprs = matrix(0,2,1)
  # Run the data.
  grn_on = kboost(KBoost::irma_on,v = v,g = g, ite= ite)
  grn_off = kboost(KBoost::irma_off,v = v,g = g, ite= ite)
  # Read the gold standard.
  g_mat = KBoost::IRMA_Gold
  # Now obtain the AUROCS and the AUPRS.
  a = AUPR_AUROC_matrix(grn_on$GRN,g_mat,FALSE,seq_len(5))
  aurocs[1] = a$AUROC
  auprs[1] = a$AUPR
  a = AUPR_AUROC_matrix(grn_off$GRN,g_mat,FALSE,seq_len(5))
  aurocs[2] = a$AUROC
  auprs[2] = a$AUPR
  return(list(aurocs = aurocs, auprs= auprs))

}
