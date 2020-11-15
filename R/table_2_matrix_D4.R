#' Function to produce the gold standard of the DREAM4 Multifactorial Challenge in matrix format.
#' @param g_table the network in table format. The first column is the Tf, the second column the gene, and the third indicates if there is an interaction.
#' @param G the number of genes.
#' @export
#'
tab_2_matrix_D4 = function(g_table,G){
  # Pre-allocate memory for the matrix.
  g_mat  = matrix(0,G,G)
  # Run a loop and set the network.
  for (i in 1:dim(g_table)[1]){
    g_mat[g_table[i,2],g_table[i,1]] = g_table[i,3]
  }
  # Return g_mat.
  return(g_mat)
}
