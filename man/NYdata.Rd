\name{NYdata}
\docType{data}
\alias{NYdata}
\alias{NYgrid}
\title{Observations of ozone concentration levels, maximum temperature and wind speed.}
\description{
       This data set contains values of daily 
	 8-hour maximum average ozone concentrations (parts per billion (ppb)),
	 maximum  temperature (in degree Celsius),
	 wind speed (knots), and 
	 relative humidity, 
	 obtained from 28 monitoring sites of New York, USA.

NYgrid: This dataset contains total 6200 rows for 62 days of observations for 10x10 = 100 grid points. 

}
\usage{NYdata}
\format{Columns for NYdata: each contains 1798 observations. 
	\itemize{
	\item   1st col = Site index (s.index), 
	\item   2nd col = Longitude,
	\item   3rd col = Latitude,  
	\item	4th col = Year,
	\item	5th col = Month,
	\item	6th col = Day, 
	\item	7th col = Ozone (o8hrmax), 
	\item	8th col = Maximum temperature (cMAXTMP), 
	\item	9th col = Wind speed (WDSP).
	\item	10th col = Relative humidity (RH).
        }
}
\seealso{
\code{\link{NYgrid}, \link{spT.Gibbs}, \link{spT.subset}}.
}
\examples{
%\dontrun{
##
  library("spTimer")
# NY data
  data(NYdata)
  head(NYdata)
# plots in NY map
  NYsite<-unique(cbind(NYdata[,1:3]))
  head(NYsite)
# map
  #library(maps)
  #map(database="state",regions="new york")
  #points(NYsite[,2:3],pch=19)

  # Grid data
  data(NYgrid)
  head(NYgrid)
  grid.coords<-unique(cbind(NYgrid[,8:9]))
  #library(maps)
  plot(grid.coords,pch=19,col=1)
  #map(database="state",regions="new york",add=TRUE)

##
%}
}
\source{US EPA}
%\references{
%http://www.epa.gov/
%We need a publication reference here, for example Sahu and Bakar (2011). Why are the data simulated? Can you explain which fields were simulated and which were kriged? Why don't you take the actual observations? We need see also NYsite. 
%}
\keyword{datasets}