% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{zoo}
\alias{zoo}
\title{Zoo data}
\format{
## `zoo`
A data frame with 100 animal observations (second instance of "frog" has 
been removed) and 17 columns:
\describe{
  \item{hair, feathers, eggs, milk, airborne, aquatic, predator, toothed, 
  backbone, breathes, venomous, fins, tail, domestic, catsize}{Boolean attributes}
  \item{legs}{Categorical attribute}
  \item{type}{Animal type, assumed to be the true clustering structure of 7 clusters}
}
}
\source{
<https://archive.ics.uci.edu/dataset/111/zoo>
}
\usage{
zoo
}
\description{
A subset of data from the UC Irvine Machine Learning Repository, created by 
Richard Forsyth. Each animal has 15 boolean attributes and 1 discrete 
attribute (legs). The animals have been partitioned into 7 "types", assumed 
to be the "true" clustering structure.
}
\keyword{datasets}
