#' Function to write output in DREAM4 Challenge Format.
#' @param GRN a GxK gene regulatory network.
#' @param TFs a K set of indixes of G that are TFs.
#' @param filename a string with the filename.
#' @export
#' @return a file with the network written as a file.
#' @examples
#' data(D4_multi_1)
#' Net = Kboost(D4_multi_1)
#' write_GRN_D4(Net$GRN, 1:100, "D4_multi_1_network.txt")
write_GRN_D4 = function(GRN,TFs, filename){
  G = dim(GRN)[1]
  K = dim(GRN)[2]
  # if no colnames and rownames are set.
  # Table form results.
  table_grn = matrix("i",G*K-K,3)
  table_no = matrix(0,G*K-K,1)
  if (is.null(rownames(GRN))){
    # Create generic names.
    gen_names = matrix("G",G,1)
    for (i in 1:G){
      gen_names[i]= paste(gen_names[i],toString(i),sep = "")
    }
  } else {
    gen_names = rownames(GRN)
  }
  # start counter
  count_  = 1
  for (j in 1:K){
      for (i in 1:G){
        # write in table.
        if (TFs[j]!=i){
          table_grn[count_, ] = c(gen_names[TFs[j]], gen_names[i] , toString(GRN[i,j]))
          table_no[count_] = GRN[i,j]
          count_ = count_  +1
        }

      }
  }

  o = order(table_no, decreasing = TRUE)

  utils::write.table(table_grn[o,], filename, col.names = FALSE, row.names = FALSE, quote = FALSE)

}
