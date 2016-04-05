\name{fitted.spT}
\alias{fitted.spT}
\title{Extract model fitted values.}
\description{
Extract average fitted values and corresponding standard deviations from model.
}
\usage{
\method{fitted}{spT}(object, ...)

## 
}
\arguments{
\item{object}{Object of class inheriting from "spT".}
\item{...}{Other arguments.}
}
\value{
\item{Mean}{Fitted mean values obtained from the MCMC samples.}
\item{SD}{Corresponding standard deviations.}
}
\seealso{
\code{\link{spT.Gibbs}}.
}
\examples{
\dontrun{
##

fitted(out) # where out is the output from spT class

##
}
}
\keyword{spT}
