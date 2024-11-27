#' runVICatMixVarSel
#'
#' Perform a run of the VICatMixVarSel model on a data-frame including variable
#' selection. Includes an option to include an outcome variable for
#' semi-supervised profile regression.
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
#'      \item{c}{A P-length vector of expected values for the variable selection 
#'      parameter, gamma.}
#'      \item{labels}{A numeric vector listing the cluster assignments for the 
#'      observations.}
#'      \item{nullphi}{A P x maxNCat matrix of maximum likelihood parameters 
#'      for irrelevant variables.}
#'      \item{rnk}{A N x K matrix of responsibilities for the latent variables 
#'      Z.}
#'    }
#'  }
#'   \item{factor_labels}{A data frame showing how variable categories
#'   correspond to numeric factor labels in the model.}
#'
#' @examples
#' # example code
#'
#' \donttest{set.seed(12)
#' generatedData <- generateSampleDataBin(500, 4, c(0.1, 0.2, 0.3, 0.4), 90, 10)
#' result <- runVICatMixVarSel(generatedData$data, 10, 0.01)
#'
#' print(result$labels)
#' #clustering labels
#'
#' print(result$c)
#' #expected values for variable selection parameter; 1 (or close to 1) indicates variable is relevant}
#'
#'
#'
#' @importFrom klaR kmodes
#' @export
runVICatMixVarSel <- function(data, K, alpha, a = 2, maxiter = 2000, tol = 0.00000005, outcome = NA, verbose = FALSE){
  
  if (!is.numeric(K) | K <= 1 | K != round(K)) {
    stop("Error: K must be an integer greater than 1.")}
  if (!is.numeric(alpha) | alpha <= 0) {
    stop("Error: alpha must be positive.")}
  if (any(is.na(data))) {
    stop("Error: The data frame contains NA values. Please remove or handle them before proceeding.")
  }
  
  #Process dataset
  if (is.na(outcome)){
    X <- as.data.frame(data)
  } else if (!is.na(outcome) & (is.vector(outcome) | is.data.frame(outcome))){
    X <- as.data.frame(cbind(data, outcome))
  } else{
    stop("y is not a vector or a column of a data frame. Please correct")
  }
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
  if (is.na(outcome)){
    X <- data.matrix(X)
  } else if (!is.na(outcome)){
    X <- data.matrix(X)
    y <- X[,dim(X)[2]]
    X <- X[,1:(dim(X)[2] - 1)]
  }
  
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
  prior$a = a
  
  #Prior for profile regression
  if (!is.na(outcome)){
    J <- length(unique(y))
    prior$beta = rep(1/J,each=J)
  }
  
  #Initialising cluster labels
  clusterInit <- klaR::kmodes(X, modes = K)$cluster #k modes: analogous to k means
  
  EPSreshape = t(prior$eps) 
  dim(EPSreshape) = c(1,maxNCat,D) 
  model = list(alpha = rep(prior$alpha, K),
               eps = EPSreshape[rep(1,K),,],
               c = rep(1, D), #initialise c_i - all variables initially included
               labels = clusterInit) 
  
  #Setting null phi to the rate of the parameter value in the dataset - max likelihood
  model$nullphi <- .nullphiCalc(X, nCat, maxNCat, D, N)
  
  model$eps <- .firstepsCalc(K, maxNCat, D, N, prior$eps, X, clusterInit)
  
  if (!is.na(outcome)){
    model$beta <- .firstbetaCalc(y, prior$beta, K, J, N, clusterInit)
  }
  
  if (is.na(outcome)){
    for (iter in 1:maxiter){
      model = .expectStepVarSel(X, model) #Expectation step
      model = .maxStepVarSel(X, model, prior) #Maximisation step
      ELBO[iter] = .ELBOCalcVarSelMStep(X, model, prior) #ELBO
      Cl[iter] = length(unique(model$labels)) #Counts number of non-empty clusters
      
      if(verbose){
        cat("Iteration number ", iter, " ELBO : ", ELBO[iter], "\n", sep = "")
      }
      if(.check_convergence(ELBO, iter, maxiter, tol)) break
      
    }
  }
  
  if (!is.na(outcome)){
    for (iter in 1:maxiter){
      model = .expectStepProfCat(X, model) #Expectation step
      model = .maxStepProfCat(X, model, prior) #Maximisation step
      ELBO[iter] = .ELBOCalcProfCatMStep(X, model, prior) #ELBO
      Cl[iter] = length(unique(model$labels)) #Counts number of non-empty clusters
      
      if(verbose){
        cat("Iteration number ", iter, " ELBO : ", ELBO[iter], "\n", sep = "")
      }
      if(.check_convergence(ELBO, iter, maxiter, tol)) break
      
    }
  }
  
  output <- list(labels = model$labels, ELBO = ELBO[1:iter], Cl = Cl[1:iter], model = model, factor_labels = factor_labels)
  return(output)
  
}

#' @keywords internal

.expectStepVarSel <- function(X, model){
  #Model should contain all current parameters: parameters alpha, epsilon, labels
  
  alpha <- model$alpha
  eps <- model$eps
  c <- model$c
  nullphi <- model$nullphi
  
  N = dim(X)[1]
  D = dim(X)[2]
  K = length(model$alpha)
  maxNCat = dim(model$eps)[[2]]
  
  Elogpi <- digamma(alpha) - digamma(sum(alpha)) #k dim vector, expectation of log pi k
  
  Elogphi <- .ElogphiCalc(eps, K, D, N, maxNCat, X)
  
  carray <- replicate(N, matrix(rep(c, K), ncol = D, byrow = TRUE), simplify="array") * Elogphi
  #array of c_i * Elogphi
  cmatrix <- .cmatrixCalc(nullphi, X, c, N, D) #1 - c_i * phi_0ixni
  
  logrhonk <- .logrhonkCalcVarSel(Elogpi, carray, cmatrix, K, D, N)
  lse <- matrixStats::rowLogSumExps(logrhonk)
  rnk <- .rnkCalc(logrhonk, lse, N, K)
  
  labels <- apply(rnk, 1, which.max) #k with the highest responsibility is the cluster that zn is assigned to
  
  model$rnk <- rnk #update responsibilities in model
  model$labels <- labels #update labels in model
  
  return(model)
  
}

#' @keywords internal
#' 

.maxStepVarSel <- function(X, model, prior){
  #Model should contain all current parameters: parameters alpha, epsilon, rnk, c, labels, need prior
  prioralpha <- prior$alpha
  prioreps <- prior$eps
  a <- prior$a
  
  rnk <- model$rnk
  c <- model$c
  nullphi <- model$nullphi
  
  N = dim(X)[1]
  D = dim(X)[2]
  K = length(model$alpha)
  maxNCat = dim(model$eps)[[2]]
  
  #Parameters for pi update - Dirichlet
  alpha <- prioralpha + colSums(rnk)
  
  #Parameters for phi update - Dirichlet
  eps <- .epsCalcVarSel(K, maxNCat, D, N, prioreps, X, rnk, c)
  
  #First calculate c_i
  Elogphi <- .ElogphiCalc(eps, K, D, N, maxNCat, X)
  lognullphi <- .lognullphiCalc(nullphi, X, K, D, N) #phi_0ixni
  
  Elogdelta <- digamma(c + a) - digamma(2*a + 1)
  Elogminusdelta <- digamma(1 - c + a) - digamma(2*a + 1)
  
  logeta1 <- as.vector(.logeta1Calc(Elogphi, rnk, Elogdelta, K, D, N))
  logeta2 <- as.vector(.logeta2Calc(lognullphi, rnk, Elogminusdelta, K, D, N))
  logetas <- matrix(c(logeta1, logeta2), nrow = D, ncol = 2)
  
  clse <- matrixStats::rowLogSumExps(logetas)
  c <- as.vector(.cCalc(logeta1, clse, D))
  
  model$alpha <- alpha #update alpha* in model
  model$eps <- eps #update epsilon* in model
  model$c <- c #update c in model
  return(model)
}

#' @keywords internal
#' 

.expectStepProfCat <- function(X, y, model){
  #Model should contain all current parameters: parameters alpha, epsilon, c, labels
  #Add parameter rnk (responsibilities) in this step; this is the first step before maximisation
  
  alpha <- model$alpha
  eps <- model$eps
  c <- model$c
  nullphi <- model$nullphi
  beta <- model$beta
  
  N = dim(X)[1]
  D = dim(X)[2]
  K = length(model$alpha)
  J = dim(beta)[2]
  maxNCat = dim(model$eps)[[2]]
  
  Elogpi <- digamma(alpha) - digamma(sum(alpha)) #k dim vector, expectation of log pi k
  Elogphi <- .ElogphiCalc(eps, K, D, N, maxNCat, X)
  Elogtheta <- .ElogthetaCalcCat(beta, K, J)
  
  carray <- replicate(N, matrix(rep(c, K), ncol = D, byrow = TRUE), simplify="array") * Elogphi
  cmatrix <- .cmatrixCalc(nullphi, X, c, N, D)
  
  logrhonk <- .logrhonkCalcProfCat(Elogpi, Elogtheta, y, carray, cmatrix, K, D, N)
  lse <- matrixStats::rowLogSumExps(logrhonk)
  rnk <- .rnkCalc(logrhonk, lse, N, K)
  
  labels <- apply(rnk, 1, which.max) 
  
  model$rnk <- rnk 
  model$labels <- labels 
  
  return(model)
  
  
}

#' @keywords internal
#' 

.maxStepProfCat <- function(X, y, model, prior){
  prioralpha <- prior$alpha
  prioreps <- prior$eps
  priorbeta <- prior$beta
  a <- prior$a
  
  rnk <- model$rnk
  c <- model$c
  nullphi <- model$nullphi
  
  N = dim(X)[1]
  D = dim(X)[2]
  K = length(model$alpha)
  J = dim(model$beta)[2]
  maxNCat = dim(model$eps)[[2]]
  
  #Parameters for pi update - Dirichlet
  alpha <- prioralpha + colSums(rnk)
  
  #Parameters for phi and theta update - Dirichlet
  eps <- .epsCalcVarSel(K, maxNCat, D, N, prioreps, X, rnk, c)
  beta <- .betaCalc(priorbeta, y, K, J, N, rnk)
  
  #First calculate c_i
  Elogphi <- .ElogphiCalc(eps, K, D, N, maxNCat, X)
  lognullphi <- .lognullphiCalc(nullphi, X, K, D, N) #phi_0ixni
  
  Elogdelta <- digamma(c + a) - digamma(2*a + 1)
  Elogminusdelta <- digamma(1 - c + a) - digamma(2*a + 1)
  
  logeta1 <- as.vector(.logeta1Calc(Elogphi, rnk, Elogdelta, K, D, N))
  logeta2 <- as.vector(.logeta2Calc(lognullphi, rnk, Elogminusdelta, K, D, N))
  logetas <- matrix(c(logeta1, logeta2), nrow = D, ncol = 2)
  
  clse <- matrixStats::rowLogSumExps(logetas)
  c <- as.vector(.cCalc(logeta1, clse, D))
  
  model$alpha <- alpha #update alpha* in model
  model$eps <- eps #update epsilon* in model
  model$c <- c #update c in model
  model$beta <- beta #update beta* in model
  return(model)
}

#' @keywords internal
#' 

.ELBOCalcVarSelMStep <- function(X, model, prior){
  N = dim(X)[1]
  D = dim(X)[2]
  K = length(model$alpha)
  maxNCat = dim(model$eps)[[2]]
  
  prior2 = list(alpha = rep(prior$alpha, K),
                eps = t(prior$eps))
  prioralpha <- prior2$alpha
  prioreps <- prior2$eps
  a <- prior$a
  
  alpha <- model$alpha
  eps <- model$eps
  rnk <- model$rnk
  c <- model$c
  nullphi <- model$nullphi
  
  Elogpi <- digamma(alpha) - digamma(sum(alpha)) 
  Elogphi <- .ElogphiCalc(eps, K, D, N, maxNCat, X)
  ElogphiL <- .ElogphiLCalc(eps, K, D, maxNCat)
  
  Elogdelta <- digamma(c + a) - digamma(2*a + 1)
  Elogminusdelta <- digamma(1 - c + a) - digamma(2*a + 1)
  
  Tk <- alpha - prioralpha
  
  #(log) normalising constants of Dirichlet
  Cprioralpha <- lgamma(sum(prioralpha)) - sum(lgamma(prioralpha))
  Cpostalpha <- lgamma(sum(alpha)) - sum(lgamma(alpha))
  Cprioreps <- .CpriorepsCalc(prioreps, K, D, maxNCat)
  Cposteps <- .CpostepsCalc(eps, K, D, maxNCat)
  Cpriordelta <- lgamma(a + a) - 2 * lgamma(a)
  Cpostdelta <- as.vector(.CpostdeltaCalc(c, a, D))
  
  #Calculations
  
  carray <- replicate(N, matrix(rep(c, K), ncol = D, byrow = TRUE), simplify="array") * Elogphi
  #array of c_i * Elogphi
  cmatrix <- .cmatrixCalc(nullphi, X, c, N, D)
  sumDElogphi <- .sumDElogphiCalcVarSel(carray, cmatrix, K, D, N)
  
  priorepsminusone <- .priorepsminusoneCalc(prioreps, K, D, maxNCat)
  epsminusone <- .epsminusoneCalc(eps, K, D, maxNCat)
  
  matExp1 <- rnk * sumDElogphi
  
  Exp1 <- sum(matExp1) #E(logp(X|Z,phi)) 
  
  Exp2 <- sum(Tk * Elogpi) #E(logp(Z|pi))
  
  Exp3 <- sum((prioralpha - 1)*Elogpi) + Cprioralpha #E(logp(pi)) 
  
  Exp4 <- sum((priorepsminusone)*ElogphiL) + sum(Cprioreps) #E(logp(phi)) 
  
  Exp5 <- sum((c * Elogdelta) + (1-c)*Elogminusdelta) #E(logp(gamma|delta)) 
  
  Exp6 <- sum((a - 1) * Elogdelta + (a - 1) * Elogminusdelta + Cpriordelta) #E(logp(delta)) 
  
  logrnk <- log(rnk)
  logrnk[logrnk == -Inf] <- 0
  Exp7 <- sum(rnk * logrnk) #E(q(Z)) 
  
  Exp8 <- sum((alpha - 1)*Elogpi) + Cpostalpha #E(logq(pi)) 
  
  Exp9 <- sum((epsminusone)*ElogphiL) + sum(Cposteps) #Elogq(phi) 
  
  matExp10 <- (c * log(c)) + ((1-c) * log(1-c))
  matExp10[is.na(matExp10)] <- 0
  
  Exp10 <- sum(matExp10) #Elogq(gamma) 
  
  Exp11 <- sum((c + a - 1) * Elogdelta + (a - c) * Elogminusdelta + Cpostdelta) #Elogq(delta) 
  
  ELBO <- Exp1 + Exp2 + Exp3 + Exp4 + Exp5 + Exp6 - Exp7 - Exp8 - Exp9 - Exp10 - Exp11
  
}

#' @keywords internal
#' 
.ELBOCalcProfCatMStep <- function(X, y, model, prior){
  N = dim(X)[1]
  D = dim(X)[2]
  J = dim(model$beta)[2]
  K = length(model$alpha)
  maxNCat = dim(model$eps)[[2]]
  
  prior2 = list(alpha = rep(prior$alpha, K),
                eps = t(prior$eps))
  prioralpha <- prior2$alpha
  prioreps <- prior2$eps
  a <- prior$a
  priorbeta <- prior$beta
  
  alpha <- model$alpha
  eps <- model$eps
  rnk <- model$rnk
  c <- model$c
  nullphi <- model$nullphi
  beta <- model$beta
  
  Elogpi <- digamma(alpha) - digamma(sum(alpha)) 
  Elogphi <- .ElogphiCalc(eps, K, D, N, maxNCat, X)
  ElogphiL <- .ElogphiLCalc(eps, K, D, maxNCat)
  
  Tk <- alpha - prioralpha
  
  Elogdelta <- digamma(c + a) - digamma(2*a + 1)
  Elogminusdelta <- digamma(1 - c + a) - digamma(2*a + 1)
  Elogtheta <- .ElogthetaCalcCat(beta, K, J)
  
  Cprioralpha <- lgamma(sum(prioralpha)) - sum(lgamma(prioralpha))
  Cpostalpha <- lgamma(sum(alpha)) - sum(lgamma(alpha))
  Cprioreps <- .CpriorepsCalc(prioreps, K, D, maxNCat)
  Cposteps <- .CpostepsCalc(eps, K, D, maxNCat)
  Cpriordelta <- lgamma(a + a) - 2 * lgamma(a)
  Cpostdelta <- .CpostdeltaCalc(c, a, D)
  Cpriortheta <- .CpriorbetaCalc(priorbeta, K, J)
  Cposttheta <- .CpostbetaCalc(beta, K, J)
  
  carray <- replicate(N, matrix(rep(c, K), ncol = D, byrow = TRUE), simplify="array") * Elogphi
  #c_i * Elogphi
  cmatrix <- .cmatrixCalc(nullphi, X, c, N, D) #c_i * phi_0ixni
  
  sumDElogphi <- .sumDElogphiCalcVarSel(carray, cmatrix, K, D, N)
  priorepsminusone <- .priorepsminusoneCalc(prioreps, K, D, maxNCat)
  epsminusone <- .epsminusoneCalc(eps, K, D, maxNCat)
  
  resptheta <- .respthetaCalc(Elogtheta, rnk, y, N, K)
  matExp1 <- rnk * sumDElogphi
  
  Exp1 <- sum(resptheta) + sum(matExp1) #E(logp(X, Y|Z,phi, theta, gamma)) 
  
  Exp2 <- sum(Tk * Elogpi) #E(logp(Z|pi))
  
  Exp3 <- sum((prioralpha - 1)*Elogpi) + Cprioralpha #E(logp(pi)) 
  
  Exp4 <- sum((priorepsminusone)*ElogphiL) + sum(Cprioreps) #E(logp(phi)) 
  
  Exp5 <- sum((c * Elogdelta) + (1-c)*Elogminusdelta) #E(logp(gamma|delta)) 
  
  Exp6 <- sum((a - 1) * Elogdelta + (a - 1) * Elogminusdelta + Cpriordelta) #E(logp(delta)) 
  
  Exp7 <- sum((priorbeta - 1) * Elogtheta) + sum(Cpriortheta) #Elogptheta 
  
  logrnk <- log(rnk)
  logrnk[logrnk == -Inf] <- 0
  Exp8 <- sum(rnk * logrnk) #E(q(Z)) 
  
  Exp9 <- sum((alpha - 1)*Elogpi) + Cpostalpha #E(logq(pi))
  
  Exp10 <- sum((epsminusone)*ElogphiL) + sum(Cposteps) #Elogq(phi) 
  
  matExp11 <- (c * log(c)) + ((1-c) * log(1-c))
  matExp11[is.na(matExp11)] <- 0
  
  Exp11 <- sum(matExp11) #Elogq(gamma) 
  
  Exp12 <- sum((c + a - 1) * Elogdelta + (a - c) * Elogminusdelta + Cpostdelta) 
  
  Exp13 <- sum((beta - 1) * Elogtheta) + sum(Cposttheta) 
  
  ELBO <- Exp1 + Exp2 + Exp3 + Exp4 + Exp5 + Exp6 + Exp7 - Exp8 - Exp9 - Exp10 - Exp11 - Exp12 - Exp13
  
}
