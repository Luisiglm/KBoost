#' Function to produce the AUPR and AUROC Results on the DREAM4 Multifactorial Challenge.
#' @param g a number larger than 0 that is the width parameter for the RBF Kernel
#' @param v a number between 0 and 1 that is the shrinkage parameter
#' @param ite an integer with number of iterations.
#' @param write_res a binary to indicate if the tables should be written.
#' @export
#'
d4_mfac = function(v,g,ite, write_res){
  # Pre-allocate memory for the results.
  aurocs = matrix(0,5,1)
  auprs = matrix(0,5,1)
  # Run the data.
  grn_1 = kboost(KBoostv::D4_multi_1,v = v,g = g, ite= ite)
  grn_2 = kboost(KBoostv::D4_multi_2,v = v,g = g, ite= ite)
  grn_3 = kboost(KBoostv::D4_multi_3,v = v,g = g, ite= ite)
  grn_4 = kboost(KBoostv::D4_multi_4,v = v,g = g, ite= ite)
  grn_5 = kboost(KBoostv::D4_multi_5,v = v,g = g, ite= ite)
  # Now format the gold standards.
  g_mat1 = tab_2_matrix_D4(KBoostv::G_D4_multi_1,100)
  g_mat2 = tab_2_matrix_D4(KBoostv::G_D4_multi_2,100)
  g_mat3 = tab_2_matrix_D4(KBoostv::G_D4_multi_3,100)
  g_mat4 = tab_2_matrix_D4(KBoostv::G_D4_multi_4,100)
  g_mat5 = tab_2_matrix_D4(KBoostv::G_D4_multi_5,100)
  # Now obtain the AUROCS and the AUPRS.
  a = AUPR_AUROC_matrix(grn_1$GRN,g_mat1,TRUE,1:100)
  aurocs[1] = a$AUROC
  auprs[1] = a$AUPR
  a = AUPR_AUROC_matrix(grn_2$GRN,g_mat2,TRUE,1:100)
  aurocs[2] = a$AUROC
  auprs[2] = a$AUPR
  a = AUPR_AUROC_matrix(grn_3$GRN,g_mat3,TRUE,1:100)
  aurocs[3] = a$AUROC
  auprs[3] = a$AUPR
  a = AUPR_AUROC_matrix(grn_4$GRN,g_mat4,TRUE,1:100)
  aurocs[4] = a$AUROC
  auprs[4] = a$AUPR
  a = AUPR_AUROC_matrix(grn_5$GRN,g_mat5,TRUE,1:100)
  aurocs[5] = a$AUROC
  auprs[5] = a$AUPR
  if (missing(write_res)){
    write_res = FALSE
  }
  if (write_res){
    file_ = "KBoost_D4_"
    write_GRN_D4(grn_1$GRN,TFs = 1:100, filename= paste(file_, "_1.txt", sep = ""))
    write_GRN_D4(grn_2$GRN,TFs = 1:100, filename =paste(file_, "_2.txt", sep = ""))
    write_GRN_D4(grn_3$GRN,TFs = 1:100, filename = paste(file_, "_3.txt", sep = ""))
    write_GRN_D4(grn_4$GRN,TFs = 1:100, filename =paste(file_, "_4.txt", sep = ""))
    write_GRN_D4(grn_5$GRN,TFs = 1:100, filename =  paste(file_, "_5.txt", sep = ""))
  }
  return(list(aurocs = aurocs, auprs= auprs))
}