\name{confint.spT}
\alias{confint.spT}
\title{Credible intervals for model parameters.}
\description{
This function is used to obtain credible intervals for model parameters from the MCMC samples.
}
\usage{
\method{confint}{spT}(object, parm, level = 0.95, ...)

## 
}
\arguments{
\item{object}{Object of class inheriting from "spT".}
\item{parm}{a specification of which parameters are to be given credible intervals, a vector of names. If missing, all parameters are considered.}
\item{level}{The required credible interval.}
\item{...}{other arguments.}
}
\seealso{
\code{\link{spT.Gibbs}}.
}
\examples{
\dontrun{
##

confint(out) # where out is the output from spT class

##
}
}
\keyword{spT}
