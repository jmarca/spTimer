\name{NYdata}
\docType{data}
\alias{NYdata}
\alias{DataFit}
\alias{DataValPred}
\alias{DataValFore}
\alias{DataFitFore}
\alias{NYgrid}
\title{Observations of ozone concentration levels, maximum temperature and wind speed.}
\description{
       This data set contains values of daily 
	 8-hour maximum average ozone concentrations (parts per billion (ppb)),
	 maximum  temperature (in degree Celsius),
	 wind speed (knots), and 
	 relative humidity, 
	 obtained from 28 monitoring sites of New York, USA.

Total 4 seperate subset datasets are created from NYdata that contains same variables as NYdata.

1. DataFit: This data is used for model fitting, contains 20 monitoring locations with 62 days observations; total 1200 rows.

2. DataValPred: This data is used for model validation for spatial prediction, contains 8 monitoring locations with 60 days observations; total 480 rows.

3. DataValFore: This data is used for model validation for temporal prediction in the unobserved locations, contains 8 monitoring locations with 2 days observations; total 8 rows.

4. DataFitFore: This data is used for model validation for temporal prediction in the observed locations, contains 20 monitoring locations with 2 days observations; total 40 rows.

We also include a dataset that contains gridded observations over NY state.

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
\code{\link{DataFit}, \link{DataFitFore}, \link{DataValFore}, \link{DataValPred}, \link{NYgrid}}.
}
\examples{
\dontrun{
##
  library("spTimer")
# NY data
  data(NYdata)
  head(NYdata)
# plots in NY map
  NYsite<-unique(cbind(NYdata[,1:3]))
  head(NYsite)
# map
  library(maps)
  map(database="state",regions="new york")
  points(NYsite[,2:3],pch=19)

# DataFit
  data(DataFit)
  head(DataFit)
# DataValPred
  data(DataValPred)
  head(DataValPred)
# DataValFore
  data(DataValFore)
  head(DataValFore)
# DataFitFore
  data(DataFitFore)
  head(DataFitFore)

# Plot fitted and validation locations in map
  fit.coords<-unique(cbind(DataFit[,1:3]))
  val.coords<-unique(cbind(DataValPred[,1:3]))

  library(maps)
  map(database="state",regions="new york")
  points(fit.coords[,2:3],pch=19,col=2)
  points(val.coords[,2:3],pch=7,col=1)
  legend(x=-78,y=41.5,pch=c(19,7),col=c(2,1),bty="n",
  legend=c("Fitted locations", "Validation locations"))

# Grid data
  data(NYgrid)
  head(NYgrid)
  grid.coords<-unique(cbind(NYgrid[,8:9]))
  library(maps)
  plot(grid.coords,pch=19,col=1)
  map(database="state",regions="new york",add=TRUE)

##
}}
\source{US EPA}
%\references{
%http://www.epa.gov/
%We need a publication reference here, for example Sahu and Bakar (2011). Why are the data simulated? Can you explain which fields were simulated and which were kriged? Why don't you take the actual observations? We need see also NYsite. 
%}
\keyword{datasets}
