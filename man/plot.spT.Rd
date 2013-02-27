\name{plot.spT}
\alias{plot.spT}
\title{Plots for spTimer output.}
\description{
This function is used to obtain MCMC summary and residual plots.
}
\usage{
## S3 method for class 'spT'
\method{plot}{spT}(x, residuals=FALSE, ...)

## 
}
\arguments{
\item{x}{Object of class inheriting from "spT".}
\item{residuals}{If TRUE then plot residual vs. fitted and normal qqplot of the residuals. If FALSE then plot MCMC samples of the parameters using coda package. Defaults value is FALSE.}
\item{...}{Other arguments.}
}
\seealso{
\code{\link{spT.Gibbs}}.
}
\examples{
\dontrun{
##

plot(out) # where out is the output from spT class
plot(out, residuals=TRUE) # where out is the output from spT class

##
}
}
\keyword{spT}
