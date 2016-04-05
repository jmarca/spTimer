\name{spT.keep.morethan.dist}
\alias{spT.keep.morethan.dist}
\alias{spT.keep.morethan.dist}
\title{Present one coordinate in a defined area for presentation}
\description{
This function is used to present one coordinate in a defined area to avoid cutter.
}
\usage{
spT.keep.morethan.dist(coords, tol.dist=100)
}
\arguments{
\item{coords}{X and Y axes/ longitude and latitude values.}
\item{tol.dist}{The tolerance limit for the distance.}
}
\seealso{
\code{\link{spT.geodist}, \link{dist}}.
}
\examples{
\dontrun{
##

data(NYdata)	
head(NYdata)
NYsite<-unique(NYdata[,2:3])	
head(NYsite)
spT.keep.morethan.dist(NYsite,tol.dist=100)

# Including values 
dat<-cbind(NYsite,value=rnorm(dim(NYsite)[[1]]))
head(dat)
spT.keep.morethan.dist(dat,tol.dist=100)

##   
}
}
\keyword{utility}