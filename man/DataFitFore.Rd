\name{DataFitFore}
\docType{data}
\alias{DataFitFore}
\title{Ozone concentration data used in Model validations for temporal forecast in the ovserved spatial locations.}
\description{
       A sub-set of the NYdata used for model validation.
}
\usage{DataFitFore}
\format{Columns: each contains 40 observations. 
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

data(DataFitFore)
head(DataFitFore)

##
}
\source{US EPA}
\keyword{datasets}
