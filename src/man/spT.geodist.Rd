\name{spT.geodist}
\alias{spT.geodist}
\alias{spT.geo.dist}
\alias{spT.geo_dist}
\title{Geodetic/geodesic Distance}
\description{
This geodetic distance provides the distance between the locations in Kilometers (k.m.) and Miles, using spherical law of Cosines.
}
\usage{
spT.geodist(Lon, Lat, KM = TRUE)

spT.geo.dist(point1, point2)
spT.geo_dist(points)
}
\arguments{
\item{Lon}{The longitude position.}
\item{Lat}{The latitude position.}
\item{KM}{A logical value, if 'TRUE' then output is in `kilometers', otherwise in `miles'.}
\item{point1}{In the form of (longitude, latitude) position.}
\item{point2}{In the form of (longitude, latitude) position.}
\item{points}{In the form of points 1:(longitude, latitude) 2:(longitude, latitude) positions.}
}
\details{
\code{spT.geodist} is used to get geodetic distance in both miles and kilometers.
\code{spT.geo.dist} is only used to get geodetic distance in kilometers with a different format.
\code{spT.geo_dist} is only used to get geodetic distance in kilometers with a different format.
}
\seealso{
\code{\link{NYdata}}, \code{\link{spT.grid.coords}}.
}
\examples{
\dontrun{
##

# Load 28 ozone monitoring locations of New York.
data(NYdata)	
head(NYdata)
NYsite<-unique(NYdata[,1:3])	

# Find the geodetic distance in km
spT.geodist(Lon=NYsite$Longitude, Lat=NYsite$Latitude, KM=TRUE)
   
# Find the geodetic distance in miles
spT.geodist(Lon=NYsite$Longitude, Lat=NYsite$Latitude, KM=FALSE)

##

# using spT.geo.dist
point1<-c(-73.757,42.681)
point2<-c(-73.881,40.866)
spT.geo.dist(point1,point2)

# using spT.geo_dist
points<-c(point1,point2)
spT.geo_dist(points)

##
}
}
\keyword{utility}