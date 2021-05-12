#' Function to perform Kernel Principal Component Boosting
#' @param X A matrix with the explanatory variables.
#' @param Y a matrix with the variable to predict.
#' @param g a positive number with the width parameter for the RBF Kernel.
#' @param v a number between 0 and 1 that corresponds to the shrinkage parameter.
#' @param ite an integer with the number of iterations.
#' @param thr a threshold to discard Kernel principal components whose eigenvalue
#' @export
#' @return function an sum of squared errors.
#' @examples
#' data(D4_multi_1)
#' Y = scale(matrix(D4_multi_1[,91],100,1))
#' X = scale(D4_multi_1[,-91])
#' res = kernel_pc_boosting(X,Y, g= 40, v = 0.5, ite = 3, thr = 1e-10)
kernel_pc_boosting <- function(X,Y,g,v,ite,thr){
    # First we will calculate the RBF kernel and the Kernel Principal Components.
    kpca <- list()
    K <- dim(X)[2]
    N <- length(Y)
    # Make sure that the input is correct.
    if (N != dim(X)[1]){
        stop("The dimensions of X and Y are different")
    }
    if (ite<1){
        stop("ite needs to be an integer larger than 0")
    }
    if (g<0){
        stop("g needs to tbe larger than 0.")
    }
    if (v<0 || v>1){
        stop("v needs to be larger than 0 but smaller than 1.")
    }
    if(thr<0 || thr>1){
        stop("v needs to be larger than 0 but smaller than 1.")
    }
    # Cool Beans pal!
    # Loop over the variables in X.
    for (j in seq_len(K)){
        # RBF Kernel.
        k <- RBF_K(X[,j],g)
        # Normalize kernel
        k <- kernel_normal(k)
        # PC on Kernel.
        kpca[[j]] <- KPC(k,thr)
    }
    # Now we will perform the Boosting Algorithm and only add the variables that perform the best.
    # Initiate f.
    f <- matrix(mean(Y),N,1)
    # Iterate in a for loop.
    for (i in seq_len(ite)){
        pse <- Y-f
        #  Perform a regression with each variable individually.
        reg <- list()
        for (j in seq_len(K)){
            # Do the regression on the pseudo residuals.
            reg[[j]] <- ort_reg(kpca[[j]],pse,v)
            if (j == 1){
                best <- j
                best_llik <- reg[[j]]$llik
            } else if (best_llik>reg[[j]]$llik){
                best <- j
                best_llik <- reg[[j]]$llik
            }
        }
        #update f.
        f <- f + v*kpca[[best]]%*%reg[[best]]$b
    }
    # get sum of squares, the object reg contains llik = log(sse/N)*(-N/2)
    sse <- exp((reg[[best]]$llik)/(-N/2))
    return(list(f  <- f, sse <- sse))
}
