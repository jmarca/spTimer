\name{plot.spT}
\alias{plot.spT}
\title{Plots for spTimer output.}
\description{
This function is used to obtain MCMC summary, residual and fitted surface plots.
}
\usage{
\method{plot}{spT}(x, residuals=FALSE, ...)

## 
}
\arguments{
\item{x}{Object of class inheriting from "spT".}
\item{residuals}{If TRUE then plot residual vs. fitted and normal qqplot of the residuals. If FALSE then plot MCMC samples of the parameters using coda package. Defaults value is FALSE.}
%\item{surface}{To obtain surface plot for the fitted mean values and corresponding standard deviations using linear interpolation, see \code{akima}. Takes values "Mean" and "SD".}
%\item{time}{Used when surface is not \code{NULL}. This argument is a vector of time point. Only one value of time should be given, e.g., time=c(5).} %First one represents the time segment one, e.g., days and second one is for the time segment two, e.g., year.}
%\item{a3d}{Used when surface is not \code{NULL}. This logical argument is to obtain a 3-dimensional plot. Default is FALSE.} 
%\item{points}{Used when surface is not \code{NULL}. This is a logical argument, when TRUE plots the coordinate points. Default is FALSE.}
%\item{title}{Used when surface is not \code{NULL}. This logical arguments is to show title for the interpolated fitted surface plot. Default is TRUE.}
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
