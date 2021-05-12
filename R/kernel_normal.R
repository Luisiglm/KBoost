#' A function to perform feature normalization in kernel space.
#'
#' @param K an NxN numeric matrix with the kernel function with N observations.
#' @export
#' @return feature centred kernel.
#' @examples
#' x = rnorm(100,0,1)
#' k = RBF_K(x,40)
#' k_ = kernel_normal(k)
kernel_normal <- function(K){
    # We first need to perform the row sum
    K_term2 <- matrix(rowSums(K),dim(K)[1],dim(K)[1])
    K_term1 <- t(K_term2)
    S <- 1/(dim(K)[1])*sum(K_term1)
    Kc <- K - 1/dim(K)[1]*(K_term1+K_term2)+1/(dim(K)[1]^2)*S
    return(Kc)
}
