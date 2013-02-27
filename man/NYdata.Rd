\name{NYdata}
\docType{data}
\alias{NYdata}
\title{Observations of ozone concentration levels, maximum temperature and wind speed.}
\description{
       This data set provides the simulated values of 
	 daily 8-hour maximum average ozone concentrations (parts per billion (ppb)),
	 maximum  temperature (in degree Celsius),
	 wind speed (knots), and relative humidity, 
	 obtained from 29 monitoring sites of New York, USA.
}
\usage{NYdata}
\format{Columns: each contains 1798 observations. 
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
\code{\link{DataFit}, \link{DataFitFore}, \link{DataValFore}, \link{DataValPred}}.
}
\examples{
\dontrun{
##

data(NYdata)
head(NYdata)

# plots in NY map

NYsite<-unique(cbind(NYdata[,1:3]))
head(NYsite)

# map
library(maps)
map(database="state",regions="new york")
points(NYsite[,2:3],pch=19)


##
}}
\source{US EPA}
%\references{
%http://www.epa.gov/
%We need a publication reference here, for example Sahu and Bakar (2011). Why are the data simulated? Can you explain which fields were simulated and which were kriged? Why don't you take the actual observations? We need see also NYsite. 
%}
\keyword{datasets}
