#' A function to format a check that the user input to KBoost is correct.
#'
#' @param X an NxG numeric matrix with the expression values of G genes and N obersvations.
#' @param TFs a Kx1 numeric matrix with integers that are the indexes of columns of X with genes that encode TFs.
#' @param g a positive scalar with the width parameter for the RBF kernel. (default g = 40)
#' @param v a double between 0 and 1 with the shrinkage parameter. (default v =0.1)
#' @param prior_weights a GxK matrix with the prior weights for the prior probability that a TF regulates a gene.(default 0.5 for all values)
#' @param ite an integer with the number of iterations (default ite = 3)
#' @export
#'
check_input = function(X,TFs,g,v,prior_weights,ite){
  ## Check Input and
  # Fill in the default values for certain parameters
   ############ CHECK MATRIX
  # Check X is a matrix
  if (class(X)[1]!="matrix"){

    stop("X needs to be an NxG matrix with N equal observations and G number of genes")

  } else if (class(X[1,1])!="numeric"){ # check that X is numeric.

    stop("The values of X need to be numeric")

  }
  ## We will scale X.
  X = scale(X)
  S = colSums(X)
  # We will check if there are any Nas or Infs.
  if (sum(is.infinite(S))>0){
    stop("After scaling some values of X are infinite.This can be because some of the variances in X are zero or some values are infinite. Check X. ")
  } else if (sum(is.nan(S))>0){
    stop("After scaling some values of X are not a number (NAN).This can be because some of the variances in X are zero or some values are infinite or NAN. Check X. ")
  }
  # Cool beans!
  G = dim(X)[2]

  N = dim(X)[1]
  ########## CHECK TFs
  # Let's check that TFs. It needs to be a vector of integers.
  # Check the length of TFs is not larger than the number of columns of X.
  if (missing(TFs)){

    TFs = 1:G
    K = length(TFs)

  } else {
      if (length(TFs)> dim(X)[2]){
        stop("TFs needs to be shorter than the number of columns of X. It is a vector of indexes of columns of X which are TFs.")
      }
      # Let's check that all members of TF are unique integers which are part of X.
      # We will do a fast check first.
      if (min(TFs)<=0){
        stop("TFs s a vector of indexes of columns of X which are TFs. There are values lower than 0, in R matrix indexing starts at 1.")
      } else if (max(TFs)>dim(X)[2]){
         stop("TFs s a vector of indexes of columns of X which are TFs. The values in TFs are larger than the number of columns.")
      }
       # check TFs are numeric or integers.
      if (class(TFs[1])!="integer"){

          stop("TFs need to be a matrix with integers corresponding to the columns in X that are TFs")

      }
     K = length(TFs)

    if (K==1){

      stop("If only 1 TF is used the network will be a vector column of  ones")

    }

  }

####### CHECK prior_Weights!
# Check if prior was specified.

if (missing(prior_weights)){

  prior_weights = matrix(0.5,G,K)

} else {

  if (min(prior_weights)<0 || max(prior_weights)>1){

    stop("the prior network has to have values between 0 and 1")

  }
  if (length(prior_weights)==1){

    if (prior_weights==0||prior_weights==1){

      stop("A prior network exactly equal to 0 or 1 for every edge would yield a posterior network equal to the prior network")

    }

    prior_weights = matrix(prior_weights,G,K)


  } else  if (length(prior_weights)!=G*K||class(prior_weights)!="matrix"||dim(prior_weights)[1]!=G||dim(prior_weights)[2]!=K){

    stop("the prior network needs to be a GxK matrix, where G is the number of genes and K the number of TFs, or a scalar between 0 and 1")

  } else Uni_P = unique(as.vector(prior_weights))

  if (length(Uni_P)<3){

    if (length(Uni_P)==1){

      if (Uni_P == 0||Uni_P==1){


        stop("Here, A prior network exactly equal to 0 or 1 for every edge would yield a posterior network equal to the prior network")

      }

    }
    }
  }

###### CHECK g!
if (missing(g)){

  g = 40

} else {

  if (class(g)!="numeric"||length(g)>1||g<=0){

    stop("g needs to be a number greater than 0")

   }
}
###### CHECK v!
if (missing(v)){

  v = 0.05
} else if (length(v)>1|| v>1|| v<=0){

  stop("v needs to be a number between 0 and 1")

}
###### CHECK ite (Maximum number of iterations)!
if (missing(ite)){

  ite = 3

} else if (ite<=0||length(ite)>1){

  stop("Max needs to be an integer greater than 0")
}

return(list(X = X, TFs= TFs, g =g, v =v, prior_weights = prior_weights, ite = ite))

}


