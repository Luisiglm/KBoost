#' Function to do a heuristic post-processing that improves accuracy. Each column is multiplied by its variance.
#'@param Net a GRN with TFs in the columns.
#'@export

net_refine = function(Net){

    for (i in 1:dim(Net)[2]){

      Net[,i] = Net[,i]*stats::var(Net[,i])

    }

    Net = Net/max(Net)

    return(Net)

}
