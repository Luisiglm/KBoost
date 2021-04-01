#' Function to do a heuristic post-processing that improves accuracy. Each column is multiplied by its variance.
#'@param Net a GRN with TFs in the columns.
#'@export
#'@return the network with Slavek and Arodz heuristic
#'@examples
#' Net =rbeta(10000,1,2)
#' Net = matrix(Net,100,100)
#' net_ref = net_refine(Net)

net_refine = function(Net){
    for (i in seq_len(dim(Net)[2])){
        Net[,i] = Net[,i]*stats::var(Net[,i])
    }
    Net = Net/max(Net)
    return(Net)
}
