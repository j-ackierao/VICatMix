#' runVICatMixVarSelAvg
#'
#' @param data - data frame/matrix with N rows of observations, and P columns of covariates
#' @param K - maximum number of clusters
#' @param alpha - Dirichlet prior parameter
#' @param a - hyperprior variable selection parameter
#' @param maxiter - maximum number of interations
#' @param tol - convergence parameter
#' @param inits - number of initialisations in the co-clustering matrix. Default 25
#' @param loss - loss function used in co-clustering matrix. Default VoIcomp. Options are "VoIavg", "VoIcomp" and "medv".
#' @param parallel - whether to use parallelisation to run initialisations
#' @param var_threshold - treshold for selection rate for determining selected variables across the averaged model. Options are 0 < n <= 1 for a threshold. Default is 0.95.
#' @param cores - Uses same default as mcclapply; user can specify number of cores for parallelisation if parallel = TRUE. Package automatically uses the user's parallel backend if one has already been registered.
#'
#' @returns A list containing the averaged labels, the averaged selected variables, as well as the results of each mixture model initialisation.
#' 
#' 
#' 
#' 
#' @importFrom mcclust.ext minVI
#' @importFrom mcclust medv
#' @importFrom mcclust comp.psm
#' @export
runVICatMixVarSelAvg <- function(data, K, alpha, a, maxiter = 2000, tol = 0.00000005, inits = 25, 
                           loss = "VoIcomp", var_threshold = 0.95, parallel = FALSE, cores=getOption('mc.cores', 2L)){
  
  resultforpsm <- list()
  resultforvars <- list()
  
  if (parallel == FALSE){
    for (i in 1:inits){
      mix <- runVICatMixVarSel(data, K, alpha, a , maxiter, tol)
      resultforpsm[[i]] <- mix$model$labels
      resultforvars[[i]] <- mix$model$c
    }
  }
  
  if (parallel == TRUE){
    if (!requireNamespace("doParallel", quietly = TRUE)) {
      stop(
        "Package 'doParallel' must be installed if parallel = TRUE",
        call. = FALSE
      )
    }
    `%dopar%` <- foreach::`%dopar%`
    `%dorng%` <- doRNG::`%dorng%`
    if (! foreach::getDoParRegistered()) {
      cl <- parallel::makePSOCKcluster(cores)
      doParallel::registerDoParallel(cl)
      message('Registered doParallel with ',
              cores, ' workers')
    } else {
      message('Using ', foreach::getDoParName(), ' with ',
              foreach::getDoParWorkers(), ' workers')
    }
    
    if (requireNamespace("doRNG", quietly = TRUE)) {
      foreach::foreach(i = 1:inits) %dorng% {
        mix <- runVICatMixVarSel(data, K, alpha, a, maxiter, tol)
        resultforpsm[[i]] <- mix$model$labels
        resultforvars[[i]] <- mix$model$c
      }
    } else {
      message("Package 'doRNG' is strongly recommended for parallelisation. Using doPar:")
      foreach::foreach(i = 1:inits) %dopar% { 
        mix <- runVICatMixVarSel(data, K, alpha, a, maxiter, tol)
        resultforpsm[[i]] <- mix$model$labels
        resultforvars[[i]] <- mix$model$c
      }
    }
  }
  
  p1 <- t(matrix(unlist(resultforpsm), dim(data)[1], inits))
  psm <- mcclust::comp.psm(p1)
  
  if (loss == "VoIavg"){
    labels_avg <- mcclust.ext::minVI(psm, method = 'avg',max.k = K)$cl
  } else if (loss == "VoIcomp"){
    labels_avg <- mcclust.ext::minVI(psm, method = 'comp',max.k = K)$cl
  } else if (loss == "medv"){
    labels_avg <- mcclust::medv(psm)
  }
  
  #Selected variables
    psmvariables <- t(matrix(unlist(resultforvars), dim(data)[2], inits))
    psmselect <- vector(mode = 'numeric', length = dim(data)[2])
    psmvariables[psmvariables > 0.5] <- 1
    psmvariables[psmvariables <= 0.5] <- 0
    for (j in 1:dim(data)[2]){
      psmselect[j] <- sum(psmvariables[,j]) / inits
    }
    psmselect[psmvariables > var_threshold] <- 1
    psmselect[psmvariables <= var_threshold] <- 0

  
  final_result <- list()
  final_result$labels_avg <- labels_avg
  final_result$varsel_avg <- psmselect
  final_result$init_results <- resultforpsm
  
  final_result
  
}