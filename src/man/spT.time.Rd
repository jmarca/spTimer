\name{spT.time}
\alias{spT.time}
\title{Timer series information.}
\description{
This function defines the time series in the spatio-temporal data.
}
\usage{
spT.time(t.series, segments=1)
}
\arguments{
\item{t.series}{Number of times within each segment in each series. It could be either a scalar or a vector. It should be a scalar if the segments are of equal length and should be a vector of length \code{segments} whose entries give the length of the segments.}  
\item{segments}{Number of segments in each time series. This should  be a scalar. }
}
\seealso{
\code{\link{spT.Gibbs}}.
}
\examples{
\dontrun{
## 

# Equal length time-series in each of 3 years
time.data<-spT.time(t.series=365,segments=3)

# Un-equal length time-series in 5 years
time.data<-spT.time(t.series=c(366, 365, 365, 365, 366),segments=5)

##
}
}
\keyword{spT}