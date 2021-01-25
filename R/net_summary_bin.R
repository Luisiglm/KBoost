#' Function to summarize the GRN filtered with a threshold,
#' the function produces a table version of the GRN, outdegree and indegree centrality, Katz centrality and closeness centrality.
#' @param GRN An inferred networks with the predictive probabilities that a transcription facor regulates a gene.
#' @param TFs A vector with indexes of the rows of GRN which correspond to TFs.
#' @param thr a scalar between 0 and 1, edges with posterior probabilities lower than thr will be discarded.
#' @param a a scalar for the Katz and PageRank centrality measures. Default the inverse of the largest eigenvalue of GRN.
#' @param b a scalar for the Katz and PageRank centrality measures. Default is 1.
#'
#' @examples
#' library(KBoost)
#' data(D4_multi_1)
#' Net = kboost(D4_multi_1)
#' Net_Summary = net_summary_bin(Net$GRN)
#'
#' @export

net_summary_bin = function(GRN,TFs,thr,a,b){

  if (class(GRN)[1]!="matrix"){

    stop("GRN needs to be an GxK matrix with G number of genes and K TFs")

  }
  # Check X values are numerical

  else if (class(GRN[1,1])=="character" ){

    stop("The values of X need to be numeric")

  }

  if (missing(TFs)){
    TFs = 1:dim(GRN)[2]
  }

  if (missing(thr)){
    thr = 0.2
  }

  if (length(TFs)!=dim(GRN)[2]){
    stop("TFs needs to be a K vector with indexes that correspond to the TFs in G.")
  }

  # check if the GRN has gene names

  if (length(rownames(GRN))==0){
    # Create generic names for genes

    gene_names = matrix("G",dim(GRN)[1],1)

    for (i in 1:dim(GRN)[1]){

      gene_names[i] = paste("G",toString(i), sep = "")
    }

  } else{
    gene_names = rownames(GRN)
  }


  # Filter network according to posterior.

  GRN = 1*(GRN>=thr)

  G = dim(GRN)[1]

  K = dim(GRN)[2]

  # Generate Table result.

  edges_GRN = matrix(GRN,dim(GRN)[1]*dim(GRN)[2],1)

  names_GRN = matrix("i",dim(GRN)[1]*dim(GRN)[2] ,2)

  block = 1:dim(GRN)[1]

  for (i in 1:dim(GRN)[2]){

    names_GRN[block,2] = gene_names
   names_GRN[block,1] = gene_names[TFs[i]]

    block = block + dim(GRN)[1]
  }

  GRN_table = data.frame(names_GRN,edges_GRN)

  colnames(GRN_table) = c("TF","Target","edge")

  o= order(GRN_table[,3], decreasing = TRUE)

  GRN_table = GRN_table[o,]

  # Calculate the Outdegree Centrality

  Outdegree = colSums(GRN)

  names(Outdegree) = gene_names[TFs]

  o = order(Outdegree,decreasing = TRUE)

  Outdegree = Outdegree[o]

  # Calculate the Indegree Centrality

  Indegree = rowSums(GRN)

  names(Indegree) = gene_names

  o = order(Indegree,decreasing = TRUE)

  Indegree = Indegree[o]



  # Calculate the Distance Matrix
  dist_mat = net_dist_bin(GRN,TFs,thr)

  rownames(dist_mat) = gene_names

  colnames(dist_mat) = gene_names[TFs]
  # Calculate the Closeness Centrality


  dist_mat_2=(1/dist_mat)

  diag(dist_mat_2) = 0

  Close_centr =  colSums(dist_mat_2)

  names(Close_centr) = gene_names[TFs]

  o = order(Close_centr,decreasing = TRUE)

  Close_centr = Close_centr[o]



return(list(GRN_table = GRN_table, Outdegree = Outdegree, Indegree = Indegree, Close_centr = Close_centr, dist_mat=  dist_mat ))

}
