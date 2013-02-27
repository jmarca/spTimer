\name{summary.spT}
\alias{summary.spT}
\title{Summary statistics of the parameters.}
\description{
This function is used to obtain MCMC summary statistics.
}
\usage{
## S3 method for class 'spT'
\method{summary}{spT}(object, pack="spTimer", ...)

## 
}
\arguments{
\item{object}{Object of class inheriting from "spT".}
\item{pack}{If "coda" then summary statistics are given using coda package. Defaults value is "spTimer".}
\item{...}{Other arguments.}
}
\seealso{
\code{\link{spT.Gibbs}}.
}
\examples{
\dontrun{
##

summary(out) # where out is the output from spT class
summary(out, pack="coda") # where out is the output from spT class

##
}
}
\keyword{spT}
