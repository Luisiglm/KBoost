#'Function for KBoost on data from a human sample annotated with Symbol names.
#' @param X an NxG matrix with the expression values of G genes and N samples.
#' @param gen_names SYMBOL gene names corresponding to the columns of X.
#' @param g a positive no., width parameter for the RBF kernel. (default g = 40)
#' @param v a double between 0 and 1, the shrinkage parameter. (default v =0.1)
#' @param ite an integer with the number of iterations (default ite = 3)
#' @param pos_weight no. between 0 and 1. Prior that a TF regulate a gene.
#' @param neg_weight no. between 0 and 1, for TF gene pairs not seen before.
#' @export
#' @return list with results of KBoost on a dataset with Symbol gene names.
#' @examples
#' X = rnorm(50,0,1)
#' X = matrix(X,10,5)
#' gen_names = c("TP53","YY1","CTCF","MDM2","ESR1")
#' grn = KBoost_human_symbol(X,gen_names,pos_weight = 0.6, neg_weight =0.4)
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
    grn = add_names(grn,gen_names)
    return(grn)
}
