#' Function to automatically assign Human TFs given a list of Symbols.
#'@param gen_names a vector or matrix with the Symbol Gene Names of the system.
#'@export
#'@return indexes of gen_names who are TFs.
#'@examples
#' gen_names = c("TP53","MDM2","FOXM1","ESR1","CTCF")
#' tfs = get_tfs_human(gen_names)
#'
get_tfs_human <- function(gen_names){
    # Load TFs data.
    Human_TFs = KBoost::Human_TFs
    # If there wasn't any match then display error message.
    TFs <- which(is.element(gen_names,Human_TFs))
    # If length TFs is zero report error.
    if (length(TFs)<=1){
        stop("We didn't find any TFs in the list of gene names submitted. Please double check that you are using gene Symbols")
    } else{
        return(TFs)
    }
}
