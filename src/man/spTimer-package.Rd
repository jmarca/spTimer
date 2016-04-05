\name{spTimer-package}
\alias{spTimer-package}
\alias{spTimer}
\docType{package}
\title{
Spatio-Temporal Bayesian Modelling using R
}
\description{
This package uses different hierarchical Bayesian spatio-temporal modelling strategies, namely: \cr
(1) Gaussian processes (GP) models, \cr
(2) Autoregressive (AR) models, \cr
(3) Gaussian predictive processes (GPP) based autoregressive models for big-n problem.
}
\details{
\tabular{ll}{
Package: \tab spTimer\cr
Type: \tab Package\cr
}
The back-end code of this package is built under c language.\cr
Main functions used: \cr
\code{> \link{spT.Gibbs}} \cr
\code{> \link{predict.spT}} \cr
Some other important functions: \cr
\code{> \link{spT.priors}} \cr
\code{> \link{spT.initials}} \cr
\code{> \link{spT.decay}} \cr
\code{> \link{spT.time}} \cr
Data descriptions: \cr
\code{> \link{NYdata}} \cr
}
\author{
K.S. Bakar & S.K. Sahu \cr
Maintainer: K.S. Bakar <shuvo.bakar@gmail.com>
}
\references{
1. Bakar, K.S. & Sahu, S.K. (2015). spTimer: Spatio-Temporal Bayesian Modelling Using R. Journal of Statistical Software. 63(15). 1-32.\cr
2. Sahu, S.K. & Bakar, K.S. (2012). Hierarchical Bayesian Autoregressive Models for Large Space Time Data with Applications to Ozone Concentration Modelling. Applied Stochastic Models in Business and Industry, 28, 395-415.\cr
3. Sahu, S.K., Gelfand, A.E., & Holland, D.M. (2007). High-Resolution Space-Time Ozone Modelling for Assessing Trends. Journal of the American Statistical Association, 102, 1221-1234.\cr
4. Bakar, K.S. (2012). Bayesian Analysis of Daily Maximum Ozone Levels. PhD Thesis, University of Southampton, Southampton, United Kingdom.
}
\keyword{ package }
\seealso{
Packages 'spacetime', 'forecast'; 'spBayes'; 'maps'; 'MBA'; 'coda'; website: \code{http://www.r-project.org/}.
}
