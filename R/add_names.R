#' Function to add names to network for the user.
#'@param grn a GRN object from KBoost.
#'@param gen_names a vector with the gene names.
#'@export
#'
add_names = function(grn,gen_names){
  # Add the gene names to the processed network.
  rownames(grn$GRN) = gen_names
  colnames(grn$GRN) = gen_names[grn$TFs]
  # Add the gene names to the un-processed network.
  rownames(grn$GRN_UP) = gen_names
  colnames(grn$GRN_UP) = gen_names[grn$TFs]
  # Add the gene names to the Prior.
  rownames(grn$prior) = gen_names
  colnames(grn$prior) = gen_names[grn$TFs]

  return(grn)
}
