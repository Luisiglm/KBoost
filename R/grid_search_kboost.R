#' Function to perform a grid search and find the hyperparameters.
#' @param dataset  1 for IRMA or 2 for DREAM4 multifactorial.
#' @param vs The range of values of v. All values need to be between 0 and 1.
#' @param gs The range of values of g. All values need to be larger than 0.
#' @param ite An integer that is the number of iterations, fixed in this case.
#' @export
#' @return list with auprs and aurocs of different values of vs and gs and ite.
#' @examples
#' res = grid_search_kboost(1,c(0.1,0.5,1),c(1,10,60,100),3)
#'
grid_search_kboost = function(dataset, vs, gs , ite){
    # Pre-allocate memory to store AUROCs and AUPRs.
    aurocs = list()
    auprs = list()
    if (dataset == 1){
        # There's two datasets in IRMA
        d_size = 2
    } else if (dataset ==2){
        d_size = 5
    }
    # We will generate colnames and rownames for clarity.
    r_n = matrix("g=",length(gs),1)
    c_n = matrix("v=", length(vs),1)
    for (i in seq_len(length(gs))){
        r_n[i] = paste(r_n[i], toString(gs[i]), sep = "")
    }
    for (j in seq_len(length(vs))){
        c_n[j] = paste(c_n[j], toString(vs[j]), sep = "")
    }
    for (i in seq_len(d_size)){
        # allocate memory.
        aurocs[[i]] = matrix(0,length(gs),length(vs))
        colnames(aurocs[[i]]) = c_n
        rownames(aurocs[[i]]) = r_n
        auprs[[i]] = matrix(0,length(gs),length(vs))
        colnames(auprs[[i]]) = c_n
        rownames(auprs[[i]]) = r_n
    }
    # Loop around the range of vs.
    for (i in seq_len(length(gs))){
        #Loop around the range of gs.
        for (j in seq_len(length(vs))){
        # Check what dataset was indicated.
            if (dataset ==2){
                res = d4_mfac(v = vs[j],g = gs[i], ite = ite)
            } else if (dataset == 1){
                res = irma_check(v = vs[j], g = gs[i], ite = ite)
            }
            # Store the results.
            for (p in seq_len(length(aurocs))){
                aurocs[[p]][i,j] = res$aurocs[p]
                auprs[[p]][i,j] = res$auprs[p]
            }
        }
    }
    return(list(auprs = auprs, aurocs = aurocs))
}
