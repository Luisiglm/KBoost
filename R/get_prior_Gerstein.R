#' Function to build a prior automatically using a previously built Network on ChIP-Seq.
#'@param gen_names the gene names in Symbol nomenclature.
#'@param TFs the indexes of gene names which are TFs.
#'@param pos_weight the prior weight for edges that were previously found in the Gerstein network.
#'@param neg_weight the prior weight for edges that were not found in the Gerstein network.
#'@export
#'
get_prior_Gerstein = function(gen_names, TFs, pos_weight, neg_weight){
  # We need to load the prior first.
  P = KBoostv::Gerstein_Prior_ENET_2
  # Great! Now, let's build the prior_weights matrix.
  prior_weights = matrix(neg_weight, length(gen_names), length(TFs))
  # Cool beans!!! Let's see if we find any previously found edges! These edges will have a higher prior probability.
  # Find unique TFs in Gerstein.
  tfs_gerstein = unique(P[,1])
  # Get the overlap between Gerstein and ours.
  gerst_ours = which(is.element(gen_names[TFs],tfs_gerstein))
  # If there's no overlap we will just return the prior_weights.
  if (length(gerst_ours)==0){
    print("This is not an error message. There was no overlap between the TFs in the prior and yours. The object <prior_weights> is a matrix with all weights equal to <neg_weight>. It can still be used for KBoost :).")
  } else {
    message = paste("There are ", toString(length(gerst_ours)), sep = " ")
    message = paste(message,"TFs that appear also in the prior network.", sep = " ")
    print(message)
    for (j in length(gerst_ours)){
      # identify the rows in P with TF j.
      tf_j = gen_names[TFs[gerst_ours[j]]]
      # Now match it with P!
      p_j = P[P[,1]==tf_j,2]
      # Match the gene names in the second column.
      match_g = which(is.element(gen_names,p_j))
      # Assing the weights in the prior.
      prior_weights[match_g,gerst_ours[j]] = pos_weight

    }
  }
  # return the prior weights.
  rownames(prior_weights) = gen_names
  colnames(prior_weights) = gen_names[TFs]
  return(prior_weights)
}
