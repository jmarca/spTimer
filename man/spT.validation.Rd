\name{spT.validation}
\alias{spT.validation}
\title{Validation Commands}
\description{
The following function is used to validate the predicted observations with the actual values.
}
\usage{
spT.validation(z, zhat, names=FALSE)
}
\arguments{
\item{z}{The original values (matrix or vector).}
\item{zhat}{The predicted values (matrix or vector).}
\item{names}{Logical, if TRUE then print the names of the validation statistics.}
}
\value{
\item{MSE}{Mean Squared Error.}
\item{RMSE}{Root Mean Squared Error.}
\item{MAE}{Mean Absolute Error.}
\item{MAPE}{Mean Absolute Percentage Error.}
\item{BIAS}{Bias.}
\item{rBIAS}{Relative Bias.}
\item{rMSEP}{Relative Mean Separation.}
}
\seealso{
\code{\link{spT.pCOVER}}.
}
\examples{
%\dontrun{
##

# Create `x', which is the true values.
# Create `y', which is the predicted values.

x <- rnorm(10,5,0.1)
y <- rnorm(10,5,1)
spT.validation(x, y)

##
%}
}
\keyword{utility}