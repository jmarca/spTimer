\name{spT.pCOVER}
\alias{spT.pCOVER}
\title{Nominal Coverage}
\description{
This function is used to obtain nominal coverage.}
\usage{
spT.pCOVER(z=NULL,zup=NULL,zlow=NULL,zsample=NULL,level=95)
}
\arguments{
\item{z}{The original values (matrix or vector).}
\item{zup}{The predicted values for upper interval (matrix or vector).}
\item{zlow}{The predicted values for lower interval (matrix or vector).}
\item{zsample}{Predicted MCMC samples.}
\item{level}{Level of coverages.}
}
\seealso{
\code{\link{spT.validation}}.
}
\examples{
\dontrun{
##

# Create `x': the true values.
# Create `yup': the upper interval.
# Create `ylow': the lower interval.

x <- rnorm(1000,5,0.1)
yup <- rnorm(1000,7,2)
ylow <- rnorm(1000,3,2)
	
# The pCOVER is:

spT.pCOVER(z=x, zup=yup, zlow=ylow)

# create predicted MCMC samples

y <- matrix(rnorm(1000*5000,5,1),1000,5000)

# The pCOVER is:

spT.pCOVER(z=x, zsample=y)
spT.pCOVER(z=x, zsample=y, level=50)

##  
}
}
\keyword{utility}