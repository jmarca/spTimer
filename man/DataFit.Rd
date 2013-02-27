\name{DataFit}
\docType{data}
\alias{DataFit}
\title{Ozone concentration data used for Model fitting.}
\description{
       A sub-set of the NYdata used for model fitting.
}
\usage{DataFit}
\format{Columns: each contains 1200 observations. 
	\itemize{
	\item   1st col = Site index (s.index), 
	\item   2nd col = Longitude,
	\item   3rd col = Latitude,  
	\item	4th col = Year (year),
	\item	5th col = Month (month),
	\item	6th col = Day (day), 
	\item	7th col = Ozone (o8hrmax), 
	\item	8th col = Maximum temperature (cMAXTMP), 
	\item	9th col = Wind speed (WDSP).
	\item	10th col = Relative humidity (RH).
        }
}
\seealso{
\code{\link{NYdata}}.
}
\examples{
##

data(DataFit)
head(DataFit)

##
}
\source{US EPA}
\keyword{datasets}
