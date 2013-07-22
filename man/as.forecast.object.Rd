\name{as.forecast.object}
\alias{as.forecast.object}
\title{Conversion of spT object into forecast object}
\description{
This function is used to convert the spT object from the temporal prediction output of "spTimer" into forecast object of the package "forecast".
}
\usage{

as.forecast.object(object, site=1, level=c(80,95), ...)

}
\arguments{
\item{object}{Object of class inheriting from "spT".}
\item{site}{Selection of location/site for which the time-series will convert into forecast object.}
\item{level}{Confidence level for temporal prediction intervals, see details in \code{\link{forecast}}.}
\item{...}{Other arguments, see details in \code{\link{forecast}}.}
}
\value{
An object class \code{"forecast"}
\item{model}{Name of the fitted model.} 
\item{method}{Name of the forecasting method.} 
\item{mean}{MCMC mean for the temporal predictions.} 
\item{lower}{Lower limits for the temporal prediction intervals obtained from the MCMC samples.} 
\item{upper}{Upper limits for the temporal prediction intervals obtained from the MCMC samples.} 
\item{level}{Prediction interval limits.} 
\item{x}{The data used for model fitting in time-series format, see \code{\link{ts}}.} 
\item{residuals}{Residuals of the fiited model.} 
\item{fitted}{Fitted MCMC mean.} 
}
\seealso{
\code{\link{spT.Gibbs}, \link{predict.spT}}.
}
\examples{
\dontrun{
##

# 'out' is the output of spT class obtained from temporal prediction
fobj<-as.forecast.object(out) 
class(fobj)
summary(fobj) 
plot(fobj) 

##
}
}
\keyword{spT}
