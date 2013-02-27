\name{spT.decay}
\alias{spT.decay}
\title{Choice for sampling spatial decay parameter \eqn{\phi}.}
\description{
This function initialises the sampling method for the spatial decay parameter \eqn{\phi}.
}
\usage{
spT.decay(type="MH", tuning=NULL, limit=NULL, segments=NULL, value=NULL)
}
\arguments{
\item{type}{The sampling method, currently available methods are, DISCRETE and MH. One can also used FIXED value for \eqn{\phi} parameter.}
\item{tuning}{If MH type is used then need to define the tuning parameter for \eqn{\phi}.}
\item{limit}{If DISCRETE type is used then need to define the lower and upper limits.}
\item{segments}{If DISCRETE type is used then need to define the number of segments for the range of limits.}
\item{value}{If FIXED type is used then need to define the value for \eqn{\phi}.}
}
\seealso{
\code{\link{spT.Gibbs}}.
}
\examples{
\dontrun{
## 

# input for random-walk Metropolis within Gibbs 
# sampling for phi parameter
spatial.decay<-spT.decay(type="MH", tuning=0.08)

# input for discrete sampling of phi parameter
spatial.decay<-spT.decay(type="DISCRETE",
           limit=c(0.01,0.02),segments=5)

# input for spatial decay if FIXED is used
spatial.decay<-spT.decay(type="FIXED", value=0.01)

##
}
}
\keyword{spT}
