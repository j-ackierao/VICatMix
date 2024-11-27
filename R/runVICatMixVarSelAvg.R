#' runVICatMixVarSelAvg
#'
#' An extension of `runVICatMixVarSel` to incorporate model
#' averaging/summarisation over multiple initialisations.
#'
#' @param data A data frame or data matrix with N rows of observations, and P
#'   columns of covariates.
#' @param K Maximum number of clusters desired. Must be an integer greater than 1.
#' @param alpha The Dirichlet prior parameter. Recommended to set this to a
#'   number < 1. Must be > 0.
#' @param a Hyperparameter for variable selection hyperprior. Default is 2.
#' @param maxiter The maximum number of iterations for the algorithm. Default is
#'   2000.
#' @param tol A convergence parameter. Default is 5x10^-8.
#' @param outcome Optional outcome variable. Default is NA; having an outcome
#'   triggers semi-supervised profile regression.
#' @param inits The number of initialisations included in the co-clustering
#'   matrix. Default is 25.
#' @param loss The loss function to be used with the co-clustering matrix.
#'   Default is VoIcomp. Options are "VoIavg", "VoIcomp" and "medv".
#' @param var_threshold Threshold for selection proportion for determining
#'   selected variables under the averaged model. Options are 0 < n <= 1 for a
#'   threshold. Default is 0.95.
#' @param parallel Logical value indicating whether to run initialisations in
#'   parallel. Default is FALSE.
#' @param cores User can specify number of cores for parallelisation if parallel
#'   = TRUE. Package automatically uses the user's parallel backend if one has
#'   already been registered.
#' @param verbose Default FALSE. Set to TRUE to output ELBO values for each iteration.
#'
#' @returns A list with the following components: (maxNCat refers to the maximum
#'   number of categories for any covariate in the data) 
#'   \item{labels_avg}{A numeric N-vector listing the cluster assignments for 
#'   the observations in the averaged model.} 
#'   \item{varsel_avg}{A numeric P-vector with a variable selection indicator 
#'   for the covariates in the averaged model.}
#'   \item{init_results}{A list where each entry is the cluster assignments for
#'   one of the initialisations included in the model averaging.}
#'   \item{init_varsel_results}{A list where each entry is the expected value
#'   for the variable selection parameters ('c') for one of the initialisations
#'   included in the model averaging.}
#'
#' @seealso \code{\link{runVICatMixVarSel}}
#'
#' @examples
#' # example code
#'
#' \donttest{set.seed(12)
#' generatedData <- generateSampleDataBin(500, 4, c(0.1, 0.2, 0.3, 0.4), 40, 10)
#' result <- runVICatMixVarSelAvg(generatedData$data, 10, 0.01, inits = 10)
#'
#' print(result$labels_avg)
#' print(result$varsel_avg)}
#'
#'
#'
#' @importFrom mcclust.ext minVI
#' @importFrom mcclust medv
#' @importFrom mcclust comp.psm
#' @export
runVICatMixVarSelAvg <- function(data, K, alpha, a = 2, maxiter = 2000, tol = 0.00000005, outcome = NA, inits = 25, 
                           loss = "VoIcomp", var_threshold = 0.95, parallel = FALSE, cores=getOption('mc.cores', 2L), verbose = FALSE){
  
  resultforpsm <- list()
  resultforvars <- list()
  
  if (parallel == FALSE){
    for (i in 1:inits){
      mix <- runVICatMixVarSel(data, K, alpha, a , maxiter, tol, outcome, verbose)
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
      mix_par <- foreach::foreach(i = 1:inits) %dorng% {
        mix <- runVICatMixVarSel(data, K, alpha, a, maxiter, tol, outcome, verbose)
      }
      resultforpsm <- lapply(mix_par, function(x) x$model$labels)
      resultforvars <- lapply(mix_par, function(x) x$model$c)
    } else {
      message("Package 'doRNG' is strongly recommended for parallelisation. Using doPar:")
      mix_par <- foreach::foreach(i = 1:inits) %dopar% { 
        mix <- runVICatMixVarSel(data, K, alpha, a, maxiter, tol, outcome, verbose)
      }
      resultforpsm <- lapply(mix_par, function(x) x$model$labels)
      resultforvars <- lapply(mix_par, function(x) x$model$c)
    }
  }
  
  p1 <- t(matrix(unlist(resultforpsm), dim(data)[1], inits))
  psm <- mcclust::comp.psm(p1)
  
  if (loss == "VoIavg"){
    labels_avg <- minVI(psm, method = 'avg',max.k = K)$cl
  } else if (loss == "VoIcomp"){
    labels_avg <- minVI(psm, method = 'comp',max.k = K)$cl
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
    psmselect[psmselect > var_threshold] <- 1
    psmselect[psmselect <= var_threshold] <- 0

  
  final_result <- list()
  final_result$labels_avg <- labels_avg
  final_result$varsel_avg <- psmselect
  final_result$init_results <- resultforpsm
  final_result$init_varsel_results <- resultforvars
  
  final_result
  
}
