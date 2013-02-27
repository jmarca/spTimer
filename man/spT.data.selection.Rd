\name{spT.data.selection}
\alias{spT.data.selection}
\title{Selection of Spatial data from a big dataset.}
\description{
This command selects a part of the big spatial dataset using the site numbers.
}
\usage{
spT.data.selection(data, random = TRUE, num.rs = NULL, 
                     s = NULL, reverse=FALSE)
}
\arguments{
\item{data}{The dataset.}
\item{random}{Logical value: if TRUE then the \code{num.rs} sites are randomly sampled, if FALSE then we need to provide the value for \code{s}.}
\item{num.rs}{The number of sites to be selected, e.g., 3.}
\item{s}{The site numbers to be selected, e.g., c(2,8,12).}
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

# Sample 4 sites randomly from the data NYdata.

r4<-spT.data.selection(data=NYdata, random=TRUE, num.rs=4)

# Choose purposively defined sites numbered as 2, 8, and 12.

p4<-spT.data.selection(data=NYdata, random=FALSE, s=c(2,8,12))
   
# Donot choose purposively defined sites numbered as 2, 8, and 12.

p4<-spT.data.selection(data=NYdata, random=FALSE, s=c(2,8,12), reverse=TRUE)
   
##
}
}
\keyword{utility}