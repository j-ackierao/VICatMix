# VICatMix: Variational Inference for Categorical Mixture Models

`VICatMix` is a variational Bayesian finite mixture model designed for the clustering of categorical data, implemented as an R package incorporating C++ (via `Rcpp` and `RcppArmadillo`) for faster computation. The package provides options to include variable selection to enhance its performance on high-dimensional or noisy data, and to incorporate model averaging and summarisation over multiple different initialisations for improved accuracy. The package additionally contains functions to generate sample clustered binary/categorical data for testing. For more details on the model, please refer to the [arXiv preprint](https://arxiv.org/abs/2406.16227).

## Installation
To install the `VICatMix` package, you can use the `devtools` package to install directly from GitHub:
```R
install.packages("devtools")
devtools::install_github("j-ackierao/VICatMix")
library(VICatMix)
```

Note VICatMix depends on the `Rcpp` and `RcppArmadillo` packages, which both require an appropriate C++ compiler.

## Examples
An example of generating sample binary data for clustering, with 'true' cluster labels.
```R
generatedData <- generateSampleDataBin(1000, 4, c(0.1, 0.2, 0.3, 0.4), 100, 0)
```
An example of running one initialisation of `VICatMix` on sample data without variable selection:
```R
result <- runVICatMix(generatedData[[1]], 10, 0.01) 
```
An example of implementing model averaging over 30 initialisations of `VICatMix` on sample data with variable selection:
```R
result <- runVICatMixVarSelAvg(generatedData[[1]], 10, 0.01, inits = 30)
```
