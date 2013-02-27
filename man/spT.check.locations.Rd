\name{spT.check.locations}
\alias{spT.check.locations}
\alias{spT.check.sites.inside}
\title{Distance Monitoring Function}
\description{
This function is used to check the minimum distance between two locations.
}
\usage{
spT.check.locations(fit.locations, pred.locations, 
   method="geodetic:km", tol = 5) 

spT.check.sites.inside(coords, method)
}
\arguments{
\item{fit.locations}{The locations for the fitted observations.}
\item{pred.locations}{The locations for the predicted observations.}
\item{method}{The distance measurement. The available methods are "geodetic:km", "geodetic:mile", "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".}
\item{tol}{The tolerance limit for the distance.}
\item{coords}{The longitude and latitude positions in a matrix format.}
}
\details{
\code{spT.check.locations} is used to check minimum distance between two locations, e.g., fitted and prediction.
\code{spT.check.sites.inside} is used to check distance within the location points. Here, the tol. limit is 0.01. If output shows nothing then we can say distances are alright.
These functions are used to avoid of occuring non-positive definite correlation matrix.
}
\seealso{
\code{\link{spT.geodist}, \link{dist}, \link{spT.data.selection}}.
}
\examples{
\dontrun{
##

data(NYdata)	
head(NYdata)
NYsite<-unique(NYdata[,1:3])	

# Sample 4 sites randomly from the data NYdata.

r4<-spT.data.selection(data=NYsite, random=TRUE, num.rs=4)

# Choose purposively sites numbered as 2, 8, and 12, 15.

p4<-spT.data.selection(data=NYsite, random=FALSE, s=c(2,8,12,15))
   
# Check locations of datasets r4 and p4

spT.check.locations(fit.locations=r4, pred.locations=p4,
    method="geodetic:km", tol=5)

#
spT.check.sites.inside(NYsite[,2:3],method="geodetic:km")

# if nothing appears then distances are alright

##   
}
}
\keyword{utility}