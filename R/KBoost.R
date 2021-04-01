#' A function to run KBoost.
#'
#' @param X an NxG matrix with the expression values of G genes and N obvs..
#' @param TFs a Kx1 numeric matrix with integers of columns of X that are TFs.
#' @param g a positive no., width parameter for RBF kernel. (default g = 40)
#' @param v a no. between 0 and 1 with the shrinkage parameter. (default v =0.1)
#' @param prior_weights it can be a scalar or GxK. (default 0.5)
#' @param ite an integer for the maximum number of iterations (default 3)
#' @export
#' @return a list with the results for kboost, with fields:
#' GRN a matrix with the posterior edge probability after network refinement.
#' GRN_UP a matrix with the posterior edges before refinement.
#' model a matrix with logical values for the TFs selected for each gene.
#' g the width parameter for the RBF kernel.
#' v the shrinkage parameter.
#' prior the prior of each model.
#' TFs a matrix with integers of each gene that is a TF.
#' prior_weights the prior_weights with which KBoost was run.
#' run_time a sacalar with the running time.
#' @examples
#' data(D4_multi_1)
#' Net = kboost(D4_multi_1)
#'
kboost = function(X,TFs,g,v,prior_weights,ite){
    # First we will check the input and make sure everything is in the right form.
    inpts = check_input(X, TFs,g,v,prior_weights,ite)
    # Now that we have the inputs. We can proceed.
    message('KBoost has checked your input and will proceed. Once we are finished, a random completion message will appear.')
    grn = kboost_main(inpts$X, inpts$TFs, inpts$g, inpts$v, inpts$prior_weights, inpts$ite)
    # Let's return grn.
    # format the results a bit if gene names were given.
    if (!is.null(colnames(X))){
        grn = add_names(grn,colnames(X))
    }
    print_completion_message()
    return(grn)
    # Print completion message!
    # Cool beans! Thanks for using KBoost! Have a great day!
}
######################################################################
##### AUXILIARY FUNCTIONS

## A function to format a check that the user input to KBoost is correct.
# X an NxG numeric matrix with the expression values of G genes and N obersvations.
# TFs a Kx1 numeric matrix with the indexes of columns of X that are TFs.
# g a positive scalar with the width parameter for the RBF kernel. (default g = 40)
# v a double between 0 and 1 with the shrinkage parameter. (default v =0.1)
# prior_weights a GxK matrix with the prior weights (default 0.5 for all values)
# ite an integer with the number of iterations (default ite = 3)
# return list with input revised and default values added
#

check_input = function(X,TFs,g,v,prior_weights,ite){
    # Check Input and
    # Fill in the default values for certain parameters
    #CHECK MATRIX
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
    # CHECK TFs
    # Let's check that TFs. It needs to be a vector of integers.
    # Check the length of TFs is not larger than the number of columns of X.
    if (missing(TFs)){
        TFs = seq_len(G)
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
    }
    if (K==1){
        stop("If only 1 TF is used the network will be a vector column of  ones")
    }
    # CHECK prior_Weights!
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
        } else  if (class(prior_weights)[1]!="matrix"||dim(prior_weights)[1]!=G||dim(prior_weights)[2]!=K){
            stop("the prior network needs to be a GxK matrix, where G is the number of genes and K the number of TFs, or a scalar between 0 and 1")
        } else {
        Uni_P = unique(as.vector(prior_weights))
        }
        if (length(Uni_P)<3 && length(Uni_P)==1 && (Uni_P == 0||Uni_P==1)){
            stop("Here, A prior network exactly equal to 0 or 1 for every edge would yield a posterior network equal to the prior network")
        }
    }
    # CHECK g!
    if (missing(g)){
        g = 60
    } else {
        if (class(g)!="numeric" || length(g)>1 || g<=0){
            stop("g needs to be a number greater than 0")
        }
    }
    # CHECK v!
    if (missing(v)){
        if (N>10){
            v = 10/N
        } else{
            v = 0.5
        }
    } else if (length(v)>1 || v>1 || v<=0){
        stop("v needs to be a number between 0 and 1")
    }
        # CHECK ite (Maximum number of iterations)!
    if (missing(ite)){
        ite = 3
    } else if (ite<=0 || length(ite)>1){
        stop("Max needs to be an integer greater than 0")
    }
    return(list(X = X, TFs= TFs, g =g, v =v, prior_weights = prior_weights, ite = ite))
}
##############################################################################
#### Function to perform boosting iterations in KBoost.
#
#  X an NxG matrix with the expression values.
# f an NxG matrix with the predicted gene expression values.
# kpca a lits with the Kernel Principal Components.
# v a double between 0 and 1 that corresponds to the shrinkage parameter.
# TFs a matrix with integers with the columns of X that are TFs.
#
# return list with log likelihoods and residuals
### greedy_uni_boosting

greedy_uni_boosting = function(f, X, TFs, kpca,v){
    # Calculate pseudo-residuals.
    pse = X-f
    # Pre-allocate memory for outputs.
    res = list()
    llik = matrix(-Inf,dim(X)[2],length(TFs))
    # Cool Beans! Now, we'll do a for loop along the TFs.
    for (j in seq_len(length(TFs))){
        # do an orthogonal regression with the KPCAs and the pseudo-residuals.
        res[[j]] = ort_reg(kpca[[j]],pse[,-TFs[j]],v)
        # get log_likelihoods.
        llik[-TFs[j],j] = res[[j]]$llik
    }
    return(list(res = res, llik = llik))
}
###############################################################################


# Function to initialize the log prior of a model given our assumptions.
# prior_weights the prior matrix with the prior for each edge.
# returns: prior the prior probaility for a model with no TFs.
#          model a logical matrix with all values false (No TFs)
## init_prior

init_prior = function(prior_weights){
    # if Prior has 1 values it is better they are removed as it will yield sets of models with -Inf probabilities.
    # we will return p_o and Model as a standard all False binary matrix.
    G = dim(prior_weights)[1]
    K = dim(prior_weights)[2]
    prior = matrix(0,G,K)
    p = rowSums(log(1-prior_weights))
    for (j in seq_len(K)){
        # Initialize a prior with all TFs missing.
        prior[,j] = p
    }
    model = matrix(FALSE, G,K)
    return(list(prior = prior, model= model))
}
###############################################################################


# Function to run the KBoost function.
# X an NxG numeric matrix with the expression values of G genes and N samples.
# TFs a Kx1 numeric matrix with  indexes of columns of X that are TFs.
# g width parameter for the RBF kernel. (default g = 40)
# v a double between 0 and 1 with the shrinkage parameter. (default v =0.1)
# prior_weights a GxK matrix with the prior weights (default 0.5 for all values)
# ite an integer with the number of iterations (default ite = 3)

## kboost main

kboost_main = function(X,TFs,g,v,prior_weights,ite){
    # do first part of the algorithm.
    # initialize f. In pre-processing we have set the mean of X to zero.
    pc = proc.time()[3]
    N = dim(X)[1]
    G = dim(X)[2]
    K = length(TFs)
    f = matrix(0,N,G)
    # cool beans! Let's do the first part.
    # the function will run a for loop, calculate RBF_kernel per TF, the Kernel principal components, and the regression results.
    part_1 = tf_kpc_reg(X,TFs,g,v)
    # part_1 is a list which contains a) kpca a list with the K, kernel principal components and res, with the regression results.
    # We can put the regression results in a matrix for clarity.
    kpca = part_1$kpca
    lliks = part_1$llik
    ## Maybe clear part_1 after??
    # Initialize the prior and model.
    p_m = init_prior(prior_weights)
    prior = p_m$prior
    model = p_m$model
    ## clear p_m after??
    # We need to make tf_check to avoid indexing issues.
    # update priors for the potential models.
    prior_new = prior_model(prior_weights, model, prior)
    # add a tf per gene and get posteriors.### CHANGE NEW VERSION WE ADD V.
    added_tfs = add_tf_get_post(prior_new,lliks,model,f, part_1$results, kpca,TFs,v)
    # cool beans! We need to store the results. we keep the models and posteriors in two lists.
    model_list =  list()
    posterior_list = list()
    model_list[[1]] = model
    model = added_tfs$model
    # update f
    f = added_tfs$f
    # We will update the model after, as this is necessary for the next iterations.
    posterior_list[[1]] = added_tfs$posteriors
    # update the prior.
    prior = matrix(added_tfs$prior_best,G,K)
    # If we have more than 1 iterations, as is usually the case, we will run a for loop a perform boosting.
    if (ite>1){
        for (i in 2:ite){
            # We will use the function greedy_uni_boosting.
            res = greedy_uni_boosting(f, X,TFs, kpca, v)
            # update priors for the potential models.
            prior_new = prior_model(prior_weights, model, prior)
            # add a tf per gene and get posteriors. The same tf can be added multiple times. ## CHANGE NEW VERSION WE ADD V.
            added_tfs = add_tf_get_post(prior_new,res$llik,model,f, res$res, kpca,TFs,v)
            # update fs.
            f = added_tfs$f
            # store old model.
            model_list[[i]] = model
            # update model.
            model = added_tfs$model
            # update posteriors list.
            posterior_list[[i]] =  added_tfs$posteriors
            # update prior.
            prior = matrix(added_tfs$prior_best,G,K)
        }
    }
    # Cool beans! Now we need to combine the iterations. we have the log_BMA function that does it.
    BMA = log_BMA(posterior_list, model_list, G,K,ite)
    # that's what we are returning son! Thanks for using KBoost, we hope you enjoyed the ride!
    # we will do the heuristic post-processing.
    BMA_proc = net_refine(BMA)
    return(list(GRN= BMA_proc, GRN_UP = BMA,model = model, g = g, v = v, prior = prior,TFs = TFs, prior_weights = prior_weights, run_time = proc.time()[3] - pc ))
}
################################################################################

# Function to add a new tf to an existing model and calculate posteriors.
#
# priors a GxK numeric matrix with the log-priors of TF regulation models.
# lliks a GxK numeric matrix with the log-marginal likelihoods of TF models.
# model a GxK binary matrix with the models at the previous iteration.
# f an NxG matrix with the predicted gene expression values.
# reg a list produced from the ort_reg function.
# kpca a list with the Kernel Principal Components.
# TFs indexes of the genes in the system integers that are TFs.
# v the shrinkage parameter.
# return list with posterior and tf to add in the next iteration.
#

add_tf_get_post = function(priors,lliks,model,f,reg, kpca, TFs,v){
    # For each gene will select the new TF with the highest posterior.
    # the number of genes are the rows.
    G = dim(priors)[1]
    # the number of TFs are the columns.
    K = dim(priors)[2]
    posteriors = lliks + priors
    # Store the priors of the best models.
    prior_best = matrix(0,G,1)
    # For each gene in G do.
    for (i in seq_len(G)){
        # Find the Tf with the highest posterior.
        best = which.max(posteriors[i,])
        # update model.
        model[i,best] = TRUE
        # update prior
        prior_best[i] = priors[i,best]
        # update f.
        # adjust best for potential indexes issues related to tf_check.
        if (TFs[best]<i){
            idx = i - 1
        } else {
            idx = i
        }
        #### USE COEFFICIENTS INSTEAD### We also added v to the input.
        f[,i] = f[,i]+v*kpca[[best]]%*%reg[[best]]$b[,idx]
    }
    return(list(model = model,posteriors = posteriors, prior_best = prior_best, f = f))
}
################################################################################

# Function to perform Bayesian Model Averaging with log un-normalized posteriors.
#
# posterior_list a list with GxK matrix with log-posteriors of TF models.
# model_list a list with GxK logical matrix indicating a TF belongs to a model.
# G the number of genes.
# K the number of TFs.
# ite the number of iterations.
# return the Bayesian model averaging result.
#
# log BMA

log_BMA = function(posterior_list, model_list, G,K,ite){
    # Pre-allocate memory for the GRN. this is the output.
    bma = matrix(0,G,K)
    # we are going to access the elements of a each list. Then concatenate them in a matrix.
    for (i in seq_len(length(posterior_list))){
        if (i ==1){
            post = posterior_list[[i]]
        } else {
            tempo = cbind(post,posterior_list[[i]])
            post = tempo
            # we will clear tempo because it is no longer used.
            rm(tempo)
        }
    }
    # cool beans.
    # We will do bma per gene.
    for (i in seq_len(G)){
        # we will remove the -infinity values cause they can cause numeric problems. They're value is zero in natural space they don't affect the sum.
        idx_inf = is.infinite(post[i,])
        p = post[i,!(idx_inf)]
        # No we have all the log-posteriors. We can do some tricks to add them in this form.
        # It is much faster to use the exponential though. However there is a chance that it will be a number under the
        # precision of the computer. We will multiply them by a constant, c, that will garantee this will not happen.
        # we will use 1e-30 as an arbitrary threshold.
        if (min(p)<log(1e-30)){
            c = log(1e-30) - min(p)
        } else {
            c = 0
        }
        # We get the exponential of the log posteriors.
        p = exp(c + p)
        # And we can get the sum.
        p = p/sum(p)
        # Make a new variable p_ with zeros.
        p_ = matrix(0,dim(post)[2],1)
        p_[!idx_inf] = p
        # Great! now we have the averaged model posteriors. We need to add them according the Tfs that belong in each model.
        j_o = 1
        j_end = K
        rang_ = j_o:j_end
        # Loop over the iterations to perform the model averaging.
        for (t in seq_len(ite)){
            # Access the model at iteration t.
            model_t = model_list[[t]]
            # Keep only the model for gene i.
            model_t = model_t[i,]
            # At each iteration we have evaluated each TF. So we will do another loop. Maybe not too efficient?
            for (j in seq_len(K)){
            # copy model_t
            model_j = model_t
            # add the TF j.
            model_j[j] = TRUE
            # add the normalized model to the model_j
            bma[i,model_j] = bma[i,model_j]+p_[rang_[j]]
        }
        # Now we need to move the window j_o and j_end.
        j_o = j_end + 1
        j_end = j_end + K
        rang_ = j_o:j_end
        }
    }
    return(bma)
}

###############################################################################
# Function to calculate the log marginal likelihood.
# Y the observed values
# F the predicted values.
# return the log marginal likelihood
#
## marginal log likelihood
log_lik= function(Y,F){
    # Calculate the sum of squared errors.
    llik = colSums((F-Y)^2)/dim(Y)[1]
    # Elevate to the power of (n/2). n is the number of observations.
    llik = log(llik)*(-(dim(Y)[1])/2)
    return(llik)
}
################################################################################
#
## function for regression when X are orthogonal and v is a shrinkage parameter.
# X a matrix with explanatory variables that are orthogonal
# Y the variable to predict.
# v a shrinkage parameter.
# returns: a list with b: the regression coefficients and the marginal loglik.
#
ort_reg = function(X,Y,v){
    # Y = X%*%b .
    # here t(X)%*%X = diag(T) (T is the number of columns of X). % denotes matrix multiplication in R.
    # b = solve(t(X)%*%X)%*%t(X)%*%Y = diag(T)%*%t(X)%*%Y.
    b = diag(dim(X)[2])%*%(t(X)%*%Y)
    f = v*X%*%b
    # Calculate the log marginal likelihood.
    llik = log_lik(f,Y)
    # Maybe we won't return f as it might be using too much memory?## CHANGE HERE ## We STOPPED RETURNING F.
    results = list(b = b, llik = llik)
    return(results)
}
################################################################################
#
# A Function to calculate the prior given a set of prior weights and a model.
# prior_weights a matrix with the prior probability to include a Tf in a model.
# model a logical matrix whose elements are TRUE if a TF belongs to a  model.
# prior the results of the prior probability of the previous iteration
# returns the prior probability for the models.
#
prior_model = function(prior_weights, model, prior){
    # At the start of each boosting iteration we will calculate the new model priors.
    G = dim(prior_weights)[1]
    K = dim(prior_weights)[2]
    # We will return prior modified.
    # we will add the log prior weight if TF is not in current and delete the log(1-Prior[i,j]).
    # Otherwise it remains the same.
    for (j in seq_len(K)){
        # For every gene we have a binary model with the TFs added in previous iterations.
        for (i in seq_len(G)){
        # If TF was NOT there before.
            if (!model[i,j]){
                prior[i,j] = prior[i,j] + log(prior_weights[i,j]) - log(1-prior_weights[i,j])
            }
        }
    }
    # Great! Now, Return prior..
    return(prior)
}
################################################################################
# function that prints completion message
#

print_completion_message = function(){
    # We have a list with ten random completion messages.
    messages = list()
    messages[[1]]= "Cool beans friend! Thanks for using KBoost."
    messages[[2]]= "KBoost has finished KBoosting your data! Thanks for using KBoost amig@!."
    messages[[3]]= "You sure sound smart for using KBoost! Thanks!"
    messages[[4]]= "KBoost has finished running.(We're feeling serious today)"
    messages[[5]]= "Gracias por usar KBoost! (Feeling Spanish today :D)"
    messages[[6]]= "Ihr Netzwerk ist bereit. Danke.(Feeling German today :D)"
    messages[[7]]= "La sua GRN e prontissima! Grazie per utilizare KBoost.(Feeling Italian today :D)"
    messages[[8]]= "La GRN esta a punt! Moltes gracies per utilizar KBoost. (Feeling Catalan today :D)"
    messages[[9]]=  "Fardig! Takk Takk! (Feeling Swedish today :D)"
    messages[[10]]= "Klaar! dankjewel vriend! (Feeling Dutch today :D)"
    messages[[11]] = "Ta ya lista la tu GRN. Muches gracies por usar KBoost! (Feeling Asturian today :D)"
    messages[[12]] = "Go raibh math agat! (Feeling Irish today :D) "
    messages[[13]] = "Oh My God Hun! You're such a star for using KBoost like I can't cope. xoxo. (We're feeling pretty today :D)"
    # Randomly choose a message.
    idx =  sample(seq_len(length(messages)),1)
    message(messages[[idx]])
}
################################################################################
# Function to do the first iteration of KBoost it performs a univariate RBF KPC regression per TF in X.
# X the gene expression matrix.
# TF a matrix with integers of which columns of X are transcription factors.
# g the width of the RBF kernel.
# v the shrinkage parameter.
# returns: the kernel principal components, the log marginal likelihoods and the coefficients in list res.
tf_kpc_reg = function(X,TFs,g,v){
    # Create a list to store the values of the KPC.
    kpca = list()
    res = list()
    llik = matrix(-Inf,dim(X)[2],length(TFs))
    for (j in seq_along(TFs)){
        # Get the RBF Kernel.
        k = RBF_K(X[,TFs[j]],g)
        # Normalize the RBF Kernel.
        k = kernel_normal(k)
        # Get the KPCs.
        kpca[[j]] = KPC(k, 1e-4)
        # Do the orthogonal regression. Exclude the Gene for the same TF. Autoregulation is not included here.
        res[[j]] = ort_reg(kpca[[j]],X[,-TFs[j]],v)
        llik[-TFs[j],j] = res[[j]]$llik
    }
    return(list(results = res, kpca = kpca, llik = llik))
}
