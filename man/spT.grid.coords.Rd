\name{spT.grid.coords}
\alias{spT.grid.coords}
\title{Grid Coordinates}
\description{
This function is used to obtain Longitude/x and Latitude/y coordinates in a grid set.
}
\usage{
spT.grid.coords(Longitude = c(max, min),
   Latitude = c(max, min), by = c(NA,NA))
}
\arguments{
\item{Longitude}{ The maximum and minimum longitude position.}
\item{Latitude}{ The maximum and minimum latitude position.}
\item{by}{The number of x and y points in each axis.}
}
\seealso{
\code{\link{spT.geodist}}.
}
\examples{
\dontrun{
##

# Load 29 ozone monitoring locations in New York.

data(NYsite)	
NYsite
coords <- as.matrix(NYsite[,c(3,2)])

# Find the knots coordinates

knots.coords <- spT.grid.coords(Longitude=c(max(coords[,1]),
          min(coords[,1])), Latitude=c(max(coords[,2]),
          min(coords[,2])),by=c(4,4))      
knots.coords

##
}
}
\keyword{utility}