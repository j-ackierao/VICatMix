\name{minVI}
\alias{minVI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Minimize the posterior expected Variation of Information
}
\description{
Finds a representative partition of the posterior by minimizing the lower bound to the posterior
expected Variation of Information from Jensen's Inequality.
}
\usage{
minVI(psm, cls.draw=NULL, method=c("avg","comp","draws","all"), 
      max.k=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{psm}{a posterior similarity matrix, which can be obtained from MCMC samples of clusterings through a call to \code{comp.psm}.}
 \item{cls.draw}{a matrix of the samples of clusterings of the \code{ncol(cls)} data points that have been used to compute \code{psm}. Note: \code{cls.draw} has to be provided if 
  \code{method="draw"} or \code{"all"}.}
 \item{method}{the optimization method used. Should be one of \code{"avg"}, \code{"comp"}, \code{"draws"}, or \code{"all"}. Defaults to \code{"avg"}.}
 \item{max.k}{ integer, if \code{method="avg"} or \code{"comp"} the maximum number of clusters up to which the hierarchical clustering is cut.
  Defaults to \code{ceiling(nrow(psm)/4)}. }
}
\details{
The Variation of Information between two clusterings is defined as the sum of the entropies minus two times the mutual information. Computation of the posterior expected Variation
of Information can be expensive, as it requires a Monte Carlo estimate. We consider a modified posterior expected Variation of Information, obtained by swapping the log and expectation, which is much more computationally efficient as it only depends on the posterior through the posterior similarity matrix. From Jensen's inequality, the problem can be viewed as minimizing a lower bound to the posterior expected loss. 

We provide several optimization methods. For \code{method="avg"} and \code{"comp"}, the search is restricted to the clusterings obtained from a hierarchical clustering with average/complete linkage and \code{1-psm} as a distance matrix (the clusterings with number of clusters \code{1:max.k} are considered).\cr
 Method \code{"draws"} restricts the search to the clusterings sampled. 
If \code{method="all"} all minimization methods are applied by default. 
}
\value{
 \item{cl}{clustering with minimal value of expected loss. If \code{method="all"} a matrix containing the clustering with the smallest value
   of the expected loss over all methods in the first row and the clusterings of the individual methods in the next rows.}
  \item{value}{value of posterior expected loss. A vector corresponding to the rows of \code{cl} if \code{method="all"}.}
  \item{method}{the optimization method used.}
}
\references{
Meila, M. (2007) Bayesian model based clustering
procedures, \emph{Journal of Multivariate Analysis} \bold{98}, 873--895.

Wade, S. and Ghahramani, Z. (2015) Bayesian cluster analysis: Point estimation and credible balls. 
Submitted. arXiv:1505.03339
}
\author{
Sara Wade, \email{sara.wade@ed.ac.uk}
}
\examples{

\donttest{set.seed(15)
generatedData <- generateSampleDataBin(100, 4, c(0.1, 0.2, 0.3, 0.4), 25, 0)
resultforpsm <- list()
for (i in 1:5){ #use 5 initialisations
  mix <- runVICatMix(generatedData$data, 10, 0.01, tol = 0.005)
  resultforpsm[[i]] <- mix$model$labels
}
p1 <- t(matrix(unlist(resultforpsm), 100, 5))
psm <- mcclust::comp.psm(p1)

labels_avg <- minVI(psm, method = 'avg',max.k = 10)$cl

print(labels_avg)}

}
