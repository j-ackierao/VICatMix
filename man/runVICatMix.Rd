% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runVICatMix.R
\name{runVICatMix}
\alias{runVICatMix}
\title{runVICatMix}
\usage{
runVICatMix(data, K, alpha, maxiter = 2000, tol = 5e-08)
}
\arguments{
\item{data}{A data frame or data matrix with N rows of observations, and P
columns of covariates.}

\item{K}{Maximum number of clusters desired.}

\item{alpha}{The Dirichlet prior parameter. Recommended to set this to a
number < 1.}

\item{maxiter}{The maximum number of iterations for the algorithm. Default is
2000.}

\item{tol}{A convergence parameter. Default is 5x10^-8.}
}
\value{
A list with the following components: (maxNCat refers to the maximum
  number of categories for any covariate in the data) 
  \item{labels}{A numeric vector listing the cluster assignments for the 
  observations.} 
  \item{ELBO}{A numeric vector tracking the value of the ELBO in every 
  iteration.}
  \item{Cl}{A numeric vector tracking the number of clusters in every 
  iteration.}
  \item{model}{A list containing all variational model parameters and the 
  cluster labels:
 \describe{
     \item{alpha}{A K-length vector of Dirichlet parameters for alpha.}
     \item{eps}{A K x maxNCat x P array of Dirichlet parameters for epsilon.}
     \item{labels}{A numeric vector listing the cluster assignments for the 
     observations.}
     \item{rnk}{A N x K matrix of responsibilities for the latent variables 
     Z.}
   }
 }
  \item{factor_labels}{A data frame showing how variable categories
  correspond to numeric factor labels in the model.}
}
\description{
Perform a run of the VICatMix model on a data-frame with no variable
selection imposed.
}
\examples{
# example code
result <- runVICatMix(zoo, 10, 0.01)

print(result$labels)



}
