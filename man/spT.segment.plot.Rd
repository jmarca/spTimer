\name{spT.segment.plot}
\alias{spT.segment.plot}
\title{Utility plot for prediction/forecast}
\description{
This function is used to obtain scatter plots with 95 percent CI for predictions/forecasts. 
}
\usage{
spT.segment.plot(obs, est, up, low, limit = NULL)
}
\arguments{
\item{obs}{Observed values.}
\item{est}{Estimated values.}
\item{up}{Upper limit of the estimated values.}
\item{low}{Lower limit of the estimated values.}
\item{limit}{x-axis and y-axis limits.}
}
\seealso{
\code{\link{summary.spT}, \link{plot.spT}}.
}
\examples{
%\dontrun{
##

obs<-rnorm(10,15,1)
est<-rnorm(10,15,1.5)
up<-rnorm(10,25,0.5)
low<-rnorm(10,5,0.5)
spT.segment.plot(obs,est,up,low,limit=c(0,30)) 

##
%}
}
\keyword{utility}