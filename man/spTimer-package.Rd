\name{spTimer-package}
\alias{spTimer-package}
\alias{spTimer}
\docType{package}
\title{
Spatio-Temporal Bayesian Modelling using R
}
\description{
This package uses different hierarchical Bayesian spatio-temporal modelling strategies, namely the gaussian processes (GP) models, the autoregressive (AR) models, and models using Gaussian predictive processes (GPP) approximation 
to analyse space-time observations.
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
\code{> \link{NYdata}} 
}
\author{
K.S. Bakar & S.K. Sahu \cr
Maintainer: K.S. Bakar <shuvo.bakar@csiro.au>
}
\references{

 Bakar, K. S. and Sahu, S. K. (2013) spTimer: Spatio-Temporal Bayesian Modelling Using R. Technical Report, University of Southampton, UK. \cr
 Sahu, S. K. and Bakar, K. S. (2012) A comparison of Bayesian Models for Daily Ozone Concentration Levels Statistical Methodology , 9, 144-157. \cr
 Sahu, S. K. and Bakar, K. S. (2012) Hierarchical Bayesian auto-regressive models for large space time data with applications to ozone concentration modelling. Applied Stochastic Models in Business and Industry, 28, 395-415. \cr
 Sahu, S. K., Bakar, K. S. and Awang, N. (2013) Bayesian Forecasting Using Hierarchical Spatio-temporal Models with Applications to Ozone Levels in the Eastern United States. Technical Report, University of Southampton. \cr
 Sahu, S.K., Gelfand, A.E., & Holland, D.M. (2007). High-Resolution Space-Time Ozone Modelling for Assessing Trends. Journal of the American Statistical Association, 102, 1221-1234.

}
\keyword{ package }
\seealso{
Packages 'forecast'; 'spBayes'; 'maps'; 'MBA'; 'coda'; website: \code{http://www.r-project.org/}.
}
%\examples{
%~~ simple examples of the most important functions ~~
%}
