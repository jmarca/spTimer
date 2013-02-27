\name{spT.time}
\alias{spT.time}
\title{Timer series information.}
\description{
This function defines the time series in the spatio-temporal data.
}
\usage{
spT.time(t.series, segment=1)
}
\arguments{
\item{t.series}{Number of times within each segment in each series. This should be a constant. Currently it is not possible to make it a variable? Like its 30 for April, and 31 for May etc.}
\item{segment}{Number of segments in each time series. This should  be a constant. }
}
\seealso{
\code{\link{spT.Gibbs}}.
}
\examples{
\dontrun{
## 

time.data<-spT.time(t.series=31,segment=2)

##
}
}
\keyword{spT}