\name{spT.decay}
\alias{spT.decay}
\title{Choice for sampling spatial decay parameter \eqn{\phi}.}
\description{
This function initialises the sampling method for the spatial decay parameter \eqn{\phi}.
}
\usage{
spT.decay(distribution=Gamm(a=2,b=1), tuning=NULL, npoints=NULL, value=NULL)
}
\arguments{
\item{distribution}{Prior distribution for \eqn{\phi}. Currently available methods are, Gamm(a,b) and Unif(low,up). One can also used "FIXED" value for \eqn{\phi} parameter.}
\item{tuning}{If the Gamma prior distribution is used then we need to define the tuning parameter for sampling \eqn{\phi}. The \code{tuning} is the standard deviation for the normal proposal distribution of the random-walk Metropolis algorithm used to sample \eqn{\phi} on the log-scale.  }
\item{npoints}{If Unif distribution is used then need to define the number of segments for the range of limits by npoints. Default value is 5.}
\item{value}{If distribution="FIXED" type is used then need to define the value for \eqn{\phi}. The default value is 3/dmax where dmax is the maximum distance between the fitting sites provided by \code{coords}. }
}
\seealso{
\code{\link{spT.Gibbs}}.
}
\examples{
\dontrun{
## 

# input for random-walk Metropolis within Gibbs 
# sampling for phi parameter
spatial.decay<-spT.decay(distribution=Gamm(2,1), tuning=0.08)

# input for discrete sampling of phi parameter 
# with uniform prior distribution
spatial.decay<-spT.decay(distribution=Unif(0.01,0.02),npoints=5)

# input for spatial decay if FIXED is used
spatial.decay<-spT.decay(distribution="FIXED", value=0.01)

##
}
}
\keyword{spT}
