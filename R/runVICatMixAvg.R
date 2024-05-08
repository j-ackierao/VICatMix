#' runVICatMixAvg
#'
#' @param data - data frame/matrix with N rows of observations, and P columns of covariates
#' @param K - maximum number of clusters
#' @param alpha - Dirichlet prior parameter
#' @param maxiter - maximum number of interations
#' @param tol - convergence parameter
#' @param inits - number of initialisations in the co-clustering matrix. Default 25
#' @param loss - loss function used in co-clustering matrix. Default VoIcomp. Options are "VoIavg", "VoIcomp" and "medv".
#' @param parallel - whether to use parallelisation to run initialisations
#' @param cores - Uses same default as mcclapply; user can specify number of cores for parallelisation if parallel = TRUE. Package automatically uses the user's parallel backend if one has already been registered.
#'
#' @returns A list containing the averaged labels, as well as the results of each mixture model initialisation.
#' 
#' 
#' 
#' 
#' @importFrom mcclust.ext minVI
#' @importFrom mcclust medv
#' @importFrom mcclust comp.psm
#' @export
runVICatMixAvg <- function(data, K, alpha, maxiter = 2000, tol = 0.00000005, inits = 25, 
                           loss = "VoIcomp", parallel = FALSE, cores=getOption('mc.cores', 2L)){
  
  resultforpsm <- list()
  
  if (parallel == FALSE){
    for (i in 1:inits){
      mix <- runVICatMix(data, K, alpha, maxiter, tol)
      resultforpsm[[i]] <- mix$model$labels
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
          mix <- runVICatMix(data, K, alpha, maxiter, tol)
          resultforpsm[[i]] <- mix$model$labels
        }
      } else {
        message("Package 'doRNG' is strongly recommended for parallelisation. Using doPar:")
        foreach::foreach(i = 1:inits) %dopar% { 
          mix <- runVICatMix(data, K, alpha, maxiter, tol)
          resultforpsm[[i]] <- mix$model$labels
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
  
  final_result <- list()
  final_result$labels_avg <- labels_avg
  final_result$init_results <- resultforpsm
  
  final_result
    
}