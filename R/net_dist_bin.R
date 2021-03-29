#' Function to calculate the distance between nodes.
#'
#' @param GRN An inferred networks with the predictive probabilities that a transcription factor regulates a gene.
#' @param TFs A vector with indexes of the rows of GRN which correspond to TFs.
#' @param thr A scalar between 0 and 1 that is used select the edges with large posterior probabilities.
#' @return a matrix with the distances between edges.
#' @examples
#' data(D4_multi_1)
#' Net = kboost(D4_multi_1)
#' dist = net_dist_bin(Net$GRN,Net$TFs,0.1)
#'
#'@export

net_dist_bin  <- function(GRN,TFs,thr){

  GRN = 1*(GRN>=thr)
  # Initialize distance matrix
  d_mat = matrix(Inf,dim(GRN)[1],dim(GRN)[2])
  diag(d_mat) = 0
  # Number of TFs
  k = length(TFs)
  # for each TF
  for (j in seq_len(k)){
    # First level of adjacency.
    d_mat[GRN[,j]>0,j] = GRN[GRN[,j]>0,j]
    # Find Targets larger than zero.
    checked = TFs[j]
    targts = which(is.element(TFs,which(GRN[,j]>0)))
    if (j>1){
      prev_dist = which(is.element(targts,TFs[seq_len(j-1)]))
      if (length(prev_dist)>0){
        for (t in seq_len(length(prev_dist))){
          d_prev = d_mat[,prev_dist[t]]
          d_prev[TFs[prev_dist[t]]] = d_mat[TFs[prev_dist[t]],j]
          d_mat[d_mat[,j]>d_prev,j] = d_prev[d_mat[,j]>d_prev]+1
        }
      }
    }
    # Add distances sequencially
    if (length(targts)>0){
      lev = 2
      Go = TRUE
        while (Go){
          for (i in seq_len(length(targts))){
            g = GRN[,targts[i]]*lev
            d_mat[d_mat[,j]>g&g>0,j] = g[(d_mat[,j]>g&g>0)]
            # Access second degree targets
            next_targts_1 =  which(is.element(TFs,which(GRN[,targts[i]]>0)))
            if (i == 1) {
              next_targts = next_targts_1
            }else {
              next_targts = unique(next_targts,next_targts_1)
            }
          }
          # update list of TFs which have been checked before.
            checked = unique(c(targts,checked))
            targts = next_targts[which(!is.element(next_targts,checked))]
          # Increase Level
            lev = lev +1
            if (length(targts)==0){
              Go = FALSE
            }
          }
        }
  }
  return(d_mat)
}
