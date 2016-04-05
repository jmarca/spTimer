\name{spT.subset}
\alias{spT.subset}
\title{Select a subset of Spatial data.}
\description{
This command selects a subset of the dataset using the site numbers.
}
\usage{
spT.subset(data, var.name, s = NULL, reverse = FALSE) 
}
\arguments{
\item{data}{The dataset.}
\item{var.name}{The name of the variable for which data will be sub-setted, e.g., "s.index".}
\item{s}{The site numbers to be selected/deselected based on the argument \code{reverse}, e.g., c(2,8,12).}
\item{reverse}{Logical value: if TRUE then \code{num.rs} will be discarded from the data.}
}
\seealso{
\code{\link{NYdata}}.
}
\examples{
\dontrun{
##

# Load ozone concentration data for New York.
data(NYdata)	
NYdata	
# Choose sites 2, 8, and 12.
subdata<-spT.subset(data=NYdata, var.name=c("s.index"), s=c(2,8,12))
# Do not choose purposively defined sites numbered as 2, 8, and 12.
subdata<-spT.subset(data=NYdata, var.name=c("s.index"), s=c(2,8,12), reverse=TRUE)
   
##
}
}
\keyword{utility}