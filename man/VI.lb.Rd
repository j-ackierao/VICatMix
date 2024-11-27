\name{VI.lb}
\alias{VI.lb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the modified Variation of Information from swapping log and expectation.
}
\description{
Based on samples of partitions (eg. from MCMC or different clustering initialisations), computes the modified Variation of Information which switches the log and expectation in the usual Variation of Information.
}
\usage{

VI.lb(cls, psm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cls}{a matrix of partitions where the posterior expected (modified) Variation of Information is to be evaluated. Each row corresponds to a clustering of \code{ncol(cls)} data points.
}
\item{psm}{a posterior similarity matrix, which can be obtained from clusterings through a call to \code{comp.psm}.}
}
\details{
The Variation of Information (VoI) between two clusterings is defined as the sum of the entropies minus two times the mutual information. Computation of the posterior expected VoI can be expensive, as it requires a Monte Carlo estimate. The modified posterior expected VoI, obtained by swapping the log and expectation, is much more computationally efficient as it only depends on the posterior through the posterior similarity matrix. From Jensen's inequality, the problem of finding the optimal partition which minimizing the posterior expected modified VoI can be viewed as minimizing a lower bound to the posterior expected VoI. 
}
\value{
vector of length \code{nrow(cls)} of the posterior expected (modified) VoI.
}
\references{
Meila, M. (2007) Bayesian model based clustering
procedures, \emph{Journal of Multivariate Analysis} \bold{98}, 873--895.

Wade, S. and Ghahramani, Z. (2015) Bayesian cluster analysis: Point estimation and credible balls. 
Submitted. arXiv:1505.03339.
}
\author{
Sara Wade, \email{sara.wade@ed.ac.uk}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{minVI}} which locates the partition that minimizes the posterior expected modified VoI.
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

# Compute modified Variation of Information for each partition from VICatMix runs
VI.lb(p1, psm)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
