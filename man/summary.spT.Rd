\name{summary.spT}
\alias{summary.spT}
\title{Summary statistics of the parameters.}
\description{
This function is used to obtain MCMC summary statistics.
}
\usage{
\method{summary}{spT}(object, digits=4, package="spTimer", coefficient=NULL, ...)

## 
}
\arguments{
\item{object}{Object of class inheriting from "spT".}
\item{digits}{Rounds the specified number of decimal places (default 4).}
\item{package}{If "coda" then summary statistics are given using coda package. Defaults value is "spTimer". One can also use "spTDyn" for obtaining spatially varying and temporal dynamic models (see details: https://cran.r-project.org/web/packages/spTDyn/index.html)}
\item{coefficient}{Only applicable for package "spTDyn".}
\item{...}{Other arguments.}
}
\value{
\item{sig2eps}{Summary statistics for \eqn{\sigma_\epsilon^2}.} 
\item{sig2eta}{Summary statistics for \eqn{\sigma_\eta^2}.} 
\item{phi}{Summary statistics for spatial decay parameter \eqn{\phi}, if estimated using \code{spT.decay}.} 
\item{...}{Summary statistics for other parameters used in the models.} 
}
\seealso{
\code{\link{spT.Gibbs}}.
}
\examples{
\dontrun{
##

summary(out) # where out is the output from spT class
summary(out, digit=2) # where out is the output from spT class
summary(out, pack="coda") # where out is the output from spT class

##
}
}
\keyword{spT}
