#' Function to add names to network for the user.
#'@param grn a GRN object from KBoost.
#'@param gen_names a vector with the gene names.
#'@export
#'@return grn a GRN object with elements with user-defined gene names.
#'@examples
#'     data(D4_multi_1)
#'     Net = kboost(D4_multi_1)
#'     g_names = matrix("G",100,1)
#'     for (i in seq_along(g_names)){
#'         g_names[i] = paste(g_names[i],toString(i), sep = "")
#'     }
#'     Net = add_names(Net,g_names)
#'
add_names <- function(grn,gen_names){
    # Add the gene names to the processed network.
    rownames(grn$GRN) <- gen_names
    colnames(grn$GRN) <- gen_names[grn$TFs]
    # Add the gene names to the un-processed network.
    rownames(grn$GRN_UP) <- gen_names
    colnames(grn$GRN_UP) <- gen_names[grn$TFs]
    # Add the gene names to the Prior.
    rownames(grn$prior) <- gen_names
    colnames(grn$prior) <- gen_names[grn$TFs]
    return(grn)
}
