#' runVICatMix
#'
#' Perform a run of the VICatMix model on a data-frame with no variable
#' selection imposed.
#'
#' @param data A data frame or data matrix with N rows of observations, and P
#'   columns of covariates.
#' @param K Maximum number of clusters desired.
#' @param alpha The Dirichlet prior parameter. Recommended to set this to a
#'   number < 1.
#' @param maxiter The maximum number of iterations for the algorithm. Default is
#'   2000.
#' @param tol A convergence parameter. Default is 5x10^-8.
#' @param verbose Default FALSE. Set to TRUE to output ELBO values for each iteration.
#'
#'
#' @returns A list with the following components: (maxNCat refers to the maximum
#'   number of categories for any covariate in the data) 
#'   \item{labels}{A numeric vector listing the cluster assignments for the 
#'   observations.} 
#'   \item{ELBO}{A numeric vector tracking the value of the ELBO in every 
#'   iteration.}
#'   \item{Cl}{A numeric vector tracking the number of clusters in every 
#'   iteration.}
#'   \item{model}{A list containing all variational model parameters and the 
#'   cluster labels:
#'  \describe{
#'      \item{alpha}{A K-length vector of Dirichlet parameters for alpha.}
#'      \item{eps}{A K x maxNCat x P array of Dirichlet parameters for epsilon.}
#'      \item{labels}{A numeric vector listing the cluster assignments for the 
#'      observations.}
#'      \item{rnk}{A N x K matrix of responsibilities for the latent variables 
#'      Z.}
#'    }
#'  }
#'   \item{factor_labels}{A data frame showing how variable categories
#'   correspond to numeric factor labels in the model.}
#'
#' @examples
#' # example code
#' result <- runVICatMix(zoo, 10, 0.01)
#'
#' print(result$labels)
#'
#'
#'
#' @importFrom klaR kmodes
#' @export
runVICatMix <- function(data, K, alpha, maxiter = 2000, tol = 0.00000005, verbose = FALSE){
  
  #Process dataset
  X <- as.data.frame(data)
  X[colnames(X)] <- lapply(X[colnames(X)], as.factor)
  
  #Data frame mapping factor labels to numbers
  factor_labels = data.frame()
  for (i in 1:length(colnames(X))){
    factorlist <- data.frame(factors = unique(X[,i]), value = as.numeric(unique(X[,i])))
    factordict <- cbind(data.frame(variable = colnames(X)[i]), factorlist)
    factor_labels <- rbind(factor_labels, factordict)
  }
  
  #Check for any completely empty factors as these may cause problems
  categories <- lapply(1:ncol(X), function(j)levels(X[, j])) 
  cat_lengths <- sapply(categories, length)
  if(any(cat_lengths == 0)) {
    stop("Column(s) ", paste0(colnames(X)[cat_lengths == 0], collapse = ","), " is empty")
  }
  
  #Create numeric matrix for data
  X <- data.matrix(X)
  
  N = dim(X)[1]
  D = dim(X)[2]
  
  nCat <- as.vector(apply(X, 2, max)) #number of categories in each variable
  maxNCat <- max(nCat)
  ELBO = Cl = rep(-Inf,maxiter*2) #record ELBO and no of clusters at each iteration
  
  #Define priors
  
  prior = list(alpha = alpha) #prior for clustering pi
  prior$eps = matrix(0, D, maxNCat) #prior for clustering phi
  for(d in 1:D){
    prior$eps [d,1:nCat[d]] <- 1/nCat[d]
  }
  
  #Initialising cluster labels
  clusterInit <- klaR::kmodes(X, modes = K)$cluster #k modes: analogous to k means
  
  EPSreshape = t(prior$eps) 
  dim(EPSreshape) = c(1,maxNCat,D) 
  model = list(alpha = rep(prior$alpha, K),
               eps = EPSreshape[rep(1,K),,],
               labels = clusterInit) 

  model$eps <- .firstepsCalc(K, maxNCat, D, N, prior$eps, X, clusterInit)
  #Update the epsilons based on the initial cluster assignment
  
  for (iter in 1:maxiter){
    model = .expectStep(X, model) #Expectation step
    .GlobalEnv$maxNCat <- maxNCat
    ELBO[iter * 2-1] = .ELBOCalc(X, model, prior) #ELBO
    Cl[iter] = length(unique(model$labels)) #Counts number of non-empty clusters
    
    model = .maxStep(X, model, prior) #Maximisation step
    ELBO[iter * 2] = .ELBOCalc(X, model, prior) #ELBO
    Cl[iter] = length(unique(model$labels)) #Counts number of non-empty clusters
    
    if(verbose){
      cat("Iteration number ", iter, ": ", ELBO[iter * 2-1], ", ", ELBO[iter * 2], "\n", sep = "")
    }
    if(.check_convergence(ELBO, iter, maxiter, tol)) break
    
  }
  
  output <- list(labels = model$labels, ELBO = ELBO[1:(iter*2)], Cl = Cl[1:iter], model = model, factor_labels = factor_labels)
  return(output)
  
}

#' @keywords internal

.expectStep <- function(X, model){
  #Model should contain all current parameters: parameters alpha, epsilon, labels
  
  alpha <- model$alpha
  eps <- model$eps
  
  N = dim(X)[1]
  D = dim(X)[2]
  K = length(model$alpha)
  nCat <- as.vector(apply(X, 2, max)) #number of categories in each variable
  maxNCat <- max(nCat)
  
  Elogpi <- digamma(alpha) - digamma(sum(alpha)) #k dim vector, expectation of log pi k
  
  Elogphi <- .ElogphiCalc(eps, K, D, N, maxNCat, X)
  
  logrhonk <- .logrhonkCalc(Elogpi, Elogphi, K, D, N) #calculate rho_nk
  lse <- matrixStats::rowLogSumExps(logrhonk)
  rnk <- .rnkCalc(logrhonk, lse, N, K)
  
  labels <- apply(rnk, 1, which.max) #k with the highest responsibility is the cluster that zn is assigned to
  
  model$rnk <- rnk #update responsibilities in model
  model$labels <- labels #update labels in model
  
  return(model)
  
}

#' @keywords internal
#' 

.maxStep <- function(X, model, prior){
  #Model should contain all current parameters: parameters alpha, epsilon, rnk, labels, need prior
  prioralpha <- prior$alpha
  prioreps <- prior$eps
  rnk <- model$rnk
  
  N = dim(X)[1]
  D = dim(X)[2]
  K = length(model$alpha)
  nCat <- as.vector(apply(X, 2, max)) #number of categories in each variable
  maxNCat <- max(nCat)
  
  #Parameters for pi update - Dirichlet
  alpha <- prioralpha + colSums(rnk)
  
  #Parameters for phi update - Dirichlet
  eps <- .epsCalc(K, maxNCat, D, N, prioreps, X, rnk)
  
  model$alpha <- alpha #update alpha* in model
  model$eps <- eps #update epsilon* in model
  return(model)
}

#' @keywords internal
#' 

.ELBOCalc <- function(X, model, prior){
  N = dim(X)[1]
  D = dim(X)[2]
  K = length(model$alpha)
  
  EPSreshape = t(prior$eps) 
  dim(EPSreshape) = c(1,maxNCat,D)
  prior2 = list(alpha = rep(prior$alpha, K),
                eps = EPSreshape[rep(1,K),,])
  prioralpha <- prior2$alpha
  prioreps <- prior2$eps
  
  alpha <- model$alpha
  eps <- model$eps
  rnk <- model$rnk
  
  nCat <- as.vector(apply(X, 2, max)) #number of categories in each variable
  maxNCat <- max(nCat)
  
  Elogpi <- digamma(alpha) - digamma(sum(alpha)) #Taken from E step
  Elogphi <- .ElogphiCalc(eps, K, D, N, maxNCat, X)
  ElogphiL <- .ElogphiLCalc(eps, K, D, maxNCat, nCat)
  
  #(log) normalising constants of Dirichlet
  Cprioralpha <- lgamma(sum(prioralpha)) - sum(lgamma(prioralpha))
  Cpostalpha <- lgamma(sum(alpha)) - sum(lgamma(alpha))
  Cprioreps <- .CpriorepsCalc(prioreps, K, D, nCat)
  Cposteps <- .CpostepsCalc(eps, K, D, nCat)
  
  #Calculations
  
  sumDElogphi <- .sumDElogphiCalc(Elogphi, K, D, N) #nxK matrix of sum of E(log phi) over i = 1, ..., D, used in calculation of first expression
  
  #Matrix of epsilon parameters -1, where all 0's remain 0 (as these are unused parameters)
  priorepsminusone <- .priorepsminusoneCalc(prioreps, K, D, maxNCat)
  epsminusone <- .epsminusoneCalc(eps, K, D, maxNCat)
  
  matExp1 <- rnk * sumDElogphi
  
  Exp1 <- sum(matExp1) #E(logp(X|Z,phi))
  
  matExp2 <- rnk * matrix(rep(Elogpi,N), ncol = K, byrow = TRUE)
  
  Exp2 <- sum(matExp2) #E(logp(Z|pi))
  
  Exp3 <- sum((prioralpha - 1)*Elogpi) + Cprioralpha #E(logp(pi)), sum will sum over all k. CHECKED
  
  Exp4 <- sum((priorepsminusone)*ElogphiL) + sum(Cprioreps) #E(logp(phi)) I think this is correct but double check ElogphiL
  
  logrnk <- log(rnk)
  logrnk[logrnk == -Inf] <- 0
  Exp5 <- sum(rnk * logrnk) #E(q(Z)) this should be correct - element-wise multiplication
  
  Exp6 <- sum((alpha - 1)*Elogpi) + Cpostalpha #E(logq(pi))
  
  Exp7 <- sum((epsminusone)*ElogphiL) + sum(Cposteps) #Elogq(phi)
  
  ELBO <- Exp1 + Exp2 + Exp3 + Exp4 - Exp5 - Exp6 - Exp7
  
}

#' @keywords internal
#' 

.check_convergence<- function(ELBO, iter, maxiter, tol){
  if (iter > 1 && abs(ELBO[iter * 2] - ELBO[iter * 2-1]) < tol && abs(ELBO[iter * 2-1] - ELBO[iter * 2-2] < tol) && 
      abs(ELBO[iter * 2-2] - ELBO[iter * 2-3] < tol)){
    cat("Stopped after iteration",iter, "\n", sep = " ") #make sure the last 3 ELBOs close to each other
    return(TRUE)
  }
  if (iter == maxiter){
    cat("Not converged after maximum number of iterations")
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}

