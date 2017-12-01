\name{spT.Gibbs}
\alias{spT.Gibbs}
\title{MCMC sampling for the spatio-temporal models.}
\description{
This function is used to draw MCMC samples using the Gibbs sampler. 
}
\usage{
spT.Gibbs(formula, data = parent.frame(), model = "GP", time.data = NULL, 
	coords, knots.coords, newcoords = NULL, newdata = NULL, priors = NULL, 
	initials = NULL, nItr = 5000, nBurn = 1000, report = 1, tol.dist = 0.05, 
	distance.method = "geodetic:km", cov.fnc = "exponential", 
	scale.transform = "NONE", spatial.decay = spT.decay(distribution = "FIXED"), 
	truncation.para = list(at = 0,lambda = 2), annual.aggrn = "NONE",
	fitted.values="TRANSFORMED")
}
\arguments{
\item{formula}{The symnbolic description of the model equation of the regression part of the space-time model.} 
\item{data}{An optional data frame containing the variables in the model. If omitted, the variables are taken from environment(formula), typically the environment from which spT.Gibbs is called. The data should be ordered first by the time and then by the sites specified by the \code{coords} below. One can also supply coordinates through this argument, where coordinate names should be "Latitude" and "Longitude".}
\item{model}{The spatio-temporal models to be fitted, current choices are: "GP", "truncatedGP", "AR", "GPP", and "truncatedGPP", with the first one as the default.}
\item{time.data}{Defining the segments of the time-series set up using the function \code{\link{spT.time}}.}
\item{coords}{The n by 2 matrix or data frame defining the locations (e.g., longitude/easting, latitude/northing) of the fitting sites, where n is the number of fitting sites. One can also supply coordinates through a formula argument such as ~Longitude+Latitude.}
\item{knots.coords}{The locations of the knots in similar format to coords above, only required if \code{model="GPP"}.}
\item{newcoords}{The locations of the prediction sites in similar format to coords above, only required if fit and predictions are to be performed simultaneously. If omitted, no predictions will be performed.}
\item{newdata}{The covariate values at the prediction sites specified by \code{newcoords}. This should have same space-time structure as the original data frame.}
\item{priors}{The prior distributions for the parameters. Default distributions are specified if these are not provided. If priors=NULL a flat prior distribution will be used with large variance. See details in \code{\link{spT.priors}}.}
\item{initials}{The preferred initial values for the parameters. If omitted, default values are provided automatically. Further  details are provided in  \code{\link{spT.initials}}.}
\item{nItr}{Number of MCMC iterations. Default value is 5000.}
\item{nBurn}{Number of burn-in samples. This number of samples will be discarded before making any inference. Default value is 1000.}
\item{report}{Number of reports to display while running the Gibbs sampler. Defaults to number of iterations. }
\item{distance.method}{The preferred method to calculate the distance between any two locations. The available options are "geodetic:km", "geodetic:mile", "euclidean", "maximum", "manhattan", and "canberra". See details in \code{\link{dist}}. The default is "geodetic:km".}
\item{tol.dist}{Minimum separation distance between any two locations out of those  specified by coords, knots.coords and pred.coords. The default is 0.005. The programme will exit if the minimum distance is less than the non-zero specified value. This will ensure non-singularity of the covariance matrices.}
\item{cov.fnc}{Covariance function for the spatial effects. The available options  are "exponential", "gaussian", "spherical" and "matern". If "matern" is used then by default the smooth parameter (\eqn{\nu}) is estimated from (0,1) uniform distribution using discrete samples.}
\item{scale.transform}{The transformation method for the response variable. Currently implemented options are: "NONE", "SQRT", and "LOG" with  "NONE" as the deault.}
\item{spatial.decay}{Provides the prior distribution for the spatial decay parameter \eqn{\phi}. Currently implemented options  are  "FIXED", "Unif", or "Gamm".  Further details for each of these are specified by \code{\link{spT.decay}}.}
\item{truncation.para}{Provides truncation parameter \eqn{\lambda} and truncation point "at" using list.}
\item{annual.aggrn}{This provides the options for calculating annual summary statistics by aggregating different time segments (e.g., annual mean). Currently implemented options are: "NONE", "ave" and "an4th", where "ave" = annual average, "an4th"= annual 4th highest. Only applicable if \code{\link{spT.time}} inputs more than one segment and when fit and predict are done simultaneously.}
\item{fitted.values}{This option provides calculating fitted values and corresponding sd in the original scale. Currently implemented options are: "ORIGINAL" and "TRANSFORMED". Only applicable if \code{scale.transform} inputs "SQRT" or "LOG". Note that the PMCC (model validation criteria) values will be changed accordingly.}
%\item{X.out}{Logical value, if TRUE then the covariates are present in the output.}
%\item{Y.out}{Logical value, if TRUE then the dependent variable is present in the output.}
}
\value{
%\item{X}{The covariate values, if X.out=TRUE.}
%\item{Y}{The observed values, if Y.out=TRUE.}
\item{accept}{The acceptance rate for the \eqn{\phi} parameter if the "MH" method of sampling is chosen.}
\item{phip}{MCMC samples for the parameter \eqn{\phi}.}
\item{nup}{MCMC samples for the parameter \eqn{\nu}. Only available if "matern" covariance function is used.}
\item{sig2eps}{MCMC samples for the parameter \eqn{\sigma^2_\epsilon}.}
\item{sig2etap}{MCMC samples for the parameter \eqn{\sigma^2_\eta}.}
%\item{sig2betap}{MCMC samples for the parameter \eqn{\sigma^2_\beta}, only applicable for spatially varying model.}
%\item{sig2deltap}{MCMC samples for the parameter \eqn{\sigma^2_\delta}, only applicable for temporally varying dynamic model.}
%\item{sig2op}{MCMC samples for the parameter \eqn{\sigma^2_o}, only applicable for temporally varying dynamic model.}
\item{betap}{MCMC samples for the parameter \eqn{\beta}.}
\item{rhop}{MCMC samples for \eqn{\rho} for the AR or GPP model.}
\item{op}{MCMC samples for the true observations.}
\item{fitted}{MCMC summary (mean and sd) for the fitted values.}
\item{tol.dist}{Minimum tolerance distance limit between the locations.}
\item{distance.method}{Name of the distance calculation method.}
\item{cov.fnc}{Name of the covariance function used in model fitting.}
\item{scale.transform}{Name of the scale.transformation method.}
\item{sampling.sp.decay}{The method of sampling for the spatial decay parameter \eqn{\phi}.}
\item{covariate.names}{Name of the covariates used in the model.}
\item{Distance.matrix}{The distance matrix.}
\item{coords}{The coordinate values.}
\item{n}{Total number of sites.}
\item{r}{Total number of segments in time, e.g., years.}
\item{T}{Total points of time, e.g., days within each year.}
\item{p}{Total number of model coefficients, i.e., \eqn{\beta}'s including the intercept.}
\item{initials}{The initial values used in the model.}
\item{priors}{The prior distributions used in the model.}
\item{PMCC}{The predictive model choice criteria obtained by minimising the expected value of a loss function, see Gelfand and Ghosh (1998). Results for both goodness of fit and penalty are given.}
\item{iterations}{The number of samples for the MCMC chain, without burn-in.}
\item{nBurn}{The number of burn-in period for the MCMC chain.}
\item{computation.time}{The computation time required for the fitted model.}
\item{model}{The spatio-temporal model used for analyse the data.}

\item{Text Output}{This option is only applicable when fit and predictions are done simultaneously.\cr

For GP models:\cr
OutGP_Values_Parameter.txt: (nItr x parameters matrix) has the MCMC samples for the parameters, ordered as: beta's, sig2eps, sig2eta, and phi. \cr
OutGP_Stats_FittedValue.txt: (N x 2) matrix of fitted summary, with 1st column as mean and 2nd column as standard deviations, where N=nrT. \cr
OutGP_Stats_PredValue.txt: ((predsites*r*T) x 2) matrix of prediction summary, with 1st column as mean and 2nd column as standard deviations. \cr
OutGP_Values_Prediction.txt: (nItr x (predsites*r*T)) matrix of MCMC predicted values in the predicted sites. \cr
If \code{annual.aggregation="ave"} then we get text output as: \cr
OutGP_Annual_Average_Prediction.txt: (nItr x (predsites*r)) matrix. \cr
If \code{annual.aggregation="an4th"} then we get text output as: \cr
OutGP_Annual_4th_Highest_Prediction.txt: (nItr x (predsites*r)) matrix. \cr


For AR models: \cr
OutAR_Values_Parameter.txt: (nItr x parameters matrix) has the MCMC samples for the parameters, ordered as: beta's, rho, sig2eps, sig2eta, mu_l's, sig2l's and phi. \cr
OutAR_Stats_TrueValue.txt: (N x 2) matrix of true summary values, with 1st column as mean and 2nd column as standard deviations. \cr
OutAR_Stats_FittedValue.txt: (N x 2) matrix of fitted summary, with 1st column as mean and 2nd column as standard deviations. \cr
OutAR_Stats_PredValue.txt: ((predsites*r*T) x 2) matrix of prediction summary, with 1st column as mean and 2nd column as standard deviations. \cr
OutAR_Values_Prediction.txt: (nItr x (predsites*r*T)) matrix of MCMC predicted values in the predicted sites. \cr
If \code{annual.aggregation="ave"} then we get text output as: \cr
OutAR_Annual_Average_Prediction.txt: (nItr x (predsites*r)) matrix. \cr
If \code{annual.aggregation="an4th"} then we get text output as: \cr
OutAR_Annual_4th_Highest_Prediction.txt: (nItr x (predsites*r)) matrix. \cr


For models using GPP approximations: \cr
OutGPP_Values_Parameter.txt: (nItr x parameters matrix) has the MCMC samples for the parameters, ordered as: beta's, rho, sig2eps, sig2eta, and phi. \cr
OutGPP_Stats_FittedValue.txt: (N x 2) matrix of fitted summary, with 1st column as mean and 2nd column as standard deviations. \cr
OutGPP_Stats_PredValue.txt: ((predsites*r*T) x 2) matrix of prediction summary, with 1st column as mean and 2nd column as standard deviations. \cr
OutGPP_Values_Prediction.txt: (nItr x (predsites*r*T)) matrix of MCMC predicted values in the predicted sites. \cr
If \code{annual.aggregation="ave"} then we get text output as: \cr
OutGPP_Annual_Average_Prediction.txt: (nItr x (predsites*r)) matrix. \cr
If \code{annual.aggregation="an4th"} then we get text output as: \cr
OutGPP_Annual_4th_Highest_Prediction.txt: (nItr x (predsites*r)) matrix. \cr

}
}
\references{

 Bakar, K. S. and Sahu, S. K. (2015) spTimer: Spatio-Temporal Bayesian Modelling Using R. Journal of Statistical Software, 63(15). 1--32. 

 Sahu, S. K. and Bakar, K. S. (2012a) A comparison of Bayesian Models for Daily Ozone Concentration Levels Statistical Methodology, 9, 144-157.

 Sahu, S. K. and Bakar, K. S. (2012b) Hierarchical Bayesian auto-regressive models for large space time data with applications to ozone concentration modelling. Applied Stochastic Models in Business and Industry, 28, 395-415. 

}
%\note{
%}
\seealso{
\code{\link{spT.priors}, \link{spT.initials}, \link{spT.geodist}, \link{dist}, \link{summary.spT}, \link{plot.spT}, \link{predict.spT}}.
}
\examples{
\dontshow{
# Read data 
data(NYdata)
#
nItr <- 100
nBurn <- 0
#
# GP: MCMC via Gibbs using default choices
set.seed(11)
post.gp <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,   
        nItr=nItr, nBurn=nBurn,
        data=NYdata, model="GP", coords=~Longitude+Latitude, 
        scale.transform="SQRT")
print(post.gp)
summary(post.gp)
# AR: MCMC via Gibbs using default choices
set.seed(11)
post.ar <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,   
        nItr=nItr,  nBurn=nBurn, 
        data=NYdata, model="AR", coords=~Longitude+Latitude, 
        scale.transform="SQRT")
print(post.ar)
summary(post.ar)
# GPP: MCMC via Gibbs using default choices
# Define the coordinates
coords<-as.matrix(unique(cbind(NYdata[,2:3])))
knots<-spT.grid.coords(Longitude=c(max(coords[,1]),
              min(coords[,1])),Latitude=c(max(coords[,2]),
              min(coords[,2])), by=c(4,4))
set.seed(11)
post.gpp <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,   
        nItr=nItr,  nBurn=nBurn,
        data=NYdata, model="GPP", coords=coords, 
        knots.coords=knots, scale.transform="SQRT")
print(post.gpp)
summary(post.gpp)
# Truncated GP: Truncation at 30 (say)
NYdata$o8hrmax[NYdata$o8hrmax<=30] <- 30
set.seed(11)
out <- spT.Gibbs(formula=o8hrmax~cMAXTMP+WDSP+RH,
  data=NYdata,
  model="truncatedGP",coords=~Longitude+Latitude,
  distance.method="geodetic:km",
  nItr=nItr,  nBurn=nBurn,
  truncation.para = list(at = 30,lambda = 4))
print(out)
summary(out)
# Truncated GPP: Truncation at 30 (say)
coords<-as.matrix(unique(cbind(NYdata[,2:3])))
knots<-spT.grid.coords(Longitude=c(max(coords[,1]),
              min(coords[,1])),Latitude=c(max(coords[,2]),
              min(coords[,2])), by=c(4,4))
set.seed(11)			  
out <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,   
         data=NYdata, 
		 model="truncatedGPP",coords=~Longitude+Latitude,
         knots.coords=knots, distance.method="geodetic:km",
         nItr=nItr,  nBurn=nBurn,
         truncation.para = list(at = 30,lambda = 4))
print(out)
summary(out)
}
\donttest{
##

###########################
## Attach library spTimer
###########################

library(spTimer)

###########################
## The GP models:
###########################

##
## Model fitting
##

# Read data 
data(NYdata)

# MCMC via Gibbs using default choices
set.seed(11)
post.gp <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,   
        data=NYdata, model="GP", coords=~Longitude+Latitude, 
        scale.transform="SQRT")
print(post.gp)

# MCMC via Gibbs not using default choices
# Read data 
s<-c(8,11,12,14,18,21,24,28)
DataFit<-spT.subset(data=NYdata, var.name=c("s.index"), s=s, reverse=TRUE) 
DataFit<-subset(DataFit, with(DataFit, !(Day \%in\% c(30, 31) & Month == 8)))
DataValPred<-spT.subset(data=NYdata, var.name=c("s.index"), s=s) 
DataValPred<-subset(DataValPred, with(DataValPred, !(Day \%in\% c(30, 31) & Month == 8)))

# define the time-series 
time.data<-spT.time(t.series=60,segment=1)

# hyper-parameters for the prior distributions
priors<-spT.priors(model="GP",inv.var.prior=Gamm(2,1),
        beta.prior=Norm(0,10^4))

# initial values for the model parameters
initials<-spT.initials(model="GP", sig2eps=0.01, 
            sig2eta=0.5, beta=NULL, phi=0.001)

# input for spatial decay, any one approach from below
#spatial.decay<-spT.decay(distribution="FIXED", value=0.01)
spatial.decay<-spT.decay(distribution=Gamm(2,1), tuning=0.08)
#spatial.decay<-spT.decay(distribution=Unif(0.01,0.02),npoints=5)

# Iterations for the MCMC algorithms
nItr<-5000

# MCMC via Gibbs
set.seed(11)
post.gp <- spT.Gibbs(formula=o8hrmax ~ cMAXTMP+WDSP+RH, 
         data=DataFit, model="GP", time.data=time.data, 
         coords=~Longitude+Latitude, priors=priors, initials=initials, 
         nItr=nItr, nBurn=0, report=nItr, 
         tol.dist=2, distance.method="geodetic:km", 
         cov.fnc="exponential", scale.transform="SQRT", 
         spatial.decay=spatial.decay)
print(post.gp)

# Summary and plots
summary(post.gp)
summary(post.gp,pack="coda")
plot(post.gp)
plot(post.gp,residuals=TRUE)

coef(post.gp)
confint(post.gp)
terms(post.gp)
formula(post.gp)
model.frame(post.gp)
model.matrix(post.gp)

# Model selection criteria
post.gp$PMCC 

##
## Fit and spatially prediction simultaneously
##

# Read data 
s<-c(8,11,12,14,18,21,24,28)
DataFit<-spT.subset(data=NYdata, var.name=c("s.index"), s=s, reverse=TRUE) 
DataFit<-subset(DataFit, with(DataFit, !(Day \%in\% c(30, 31) & Month == 8)))
DataValPred<-spT.subset(data=NYdata, var.name=c("s.index"), s=s) 
DataValPred<-subset(DataValPred, with(DataValPred, !(Day \%in\% c(30, 31) & Month == 8)))
# Define the coordinates
coords<-as.matrix(unique(cbind(DataFit[,2:3])))
pred.coords<-as.matrix(unique(cbind(DataValPred[,2:3])))
# MCMC via Gibbs will provide output in *.txt format  
# from C routine to avoide large data problem in R
set.seed(11)
post.gp.fitpred <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,   
         data=DataFit, model="GP", coords=coords, 
         newcoords=pred.coords, newdata=DataValPred,
         scale.transform="SQRT")
print(post.gp.fitpred)
summary(post.gp.fitpred)
coef(post.gp.fitpred)
plot(post.gp.fitpred)
names(post.gp.fitpred)
# validation criteria
spT.validation(DataValPred$o8hrmax,c(post.gp.fitpred$prediction[,1]))  


######################################
## The GP model for sp class data
######################################

# Creating sp class data
library(sp)
data(meuse)
summary(meuse)
coordinates(meuse) <- ~x+y
class(meuse)
out<-spT.Gibbs(formula=zinc~sqrt(dist),data=meuse,
               model="GP", scale.transform="LOG")
summary(out)

# Create a dataset with spacetime class
library(spTimer)
site<-unique(NYdata[,c("Longitude","Latitude")])
library(spacetime)
row.names(site)<-paste("point",1:nrow(site),sep="")
site <- SpatialPoints(site)
ymd<-as.POSIXct(seq(as.Date("2006-07-01"),as.Date("2006-08-31"),by=1))
# introduce class STFDF
newNYdata<-STFDF(sp=site, time=ymd, data=NYdata) # full lattice
class(newNYdata)
out <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,   
      data=newNYdata, model="GP", scale.transform="SQRT")
summary(out)


###########################
## The AR models:
###########################

##
## Model fitting
##

# Read data 
data(NYdata)

# Define the coordinates
coords<-as.matrix(unique(cbind(NYdata[,2:3])))

# MCMC via Gibbs using default choices
set.seed(11)
post.ar <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,   
         data=NYdata, model="AR", coords=coords, 
         scale.transform="SQRT")
print(post.ar)

# MCMC via Gibbs not using default choices
# define the time-series 
time.data<-spT.time(t.series=62,segment=1)

# hyper-parameters for the prior distributions
priors<-spT.priors(model="AR",inv.var.prior=Gamm(2,1),
        beta.prior=Norm(0,10^4))

# initial values for the model parameters
initials<-spT.initials(model="AR", sig2eps=0.01, 
            sig2eta=0.5, beta=NULL, phi=0.001)

# Input for spatial decay
#spatial.decay<-spT.decay(distribution="FIXED", value=0.01)
spatial.decay<-spT.decay(distribution=Gamm(2,1), tuning=0.08)
#spatial.decay<-spT.decay(distribution=Unif(0.01,0.02),npoints=5)

# Iterations for the MCMC algorithms
nItr<-5000

# MCMC via Gibbs
set.seed(11)
post.ar <- spT.Gibbs(formula=o8hrmax~cMAXTMP+WDSP+RH, 
         data=NYdata, model="AR", time.data=time.data, 
         coords=coords, priors=priors, initials=initials, 
         nItr=nItr, nBurn=0, report=nItr, 
         tol.dist=2, distance.method="geodetic:km", 
         cov.fnc="exponential", scale.transform="SQRT", 
         spatial.decay=spatial.decay)
print(post.ar)

# Summary and plots
summary(post.ar)
plot(post.ar)

# Model selection criteria
post.ar$PMCC 

##
## Fit and spatially prediction simultaneously
##

# Read data 
s<-c(8,11,12,14,18,21,24,28)
DataFit<-spT.subset(data=NYdata, var.name=c("s.index"), s=s, reverse=TRUE) 
DataFit<-subset(DataFit, with(DataFit, !(Day \%in\% c(30, 31) & Month == 8)))
DataValPred<-spT.subset(data=NYdata, var.name=c("s.index"), s=s) 
DataValPred<-subset(DataValPred, with(DataValPred, !(Day \%in\% c(30, 31) & Month == 8)))
# Define the coordinates
coords<-as.matrix(unique(cbind(DataFit[,2:3])))
pred.coords<-as.matrix(unique(cbind(DataValPred[,2:3])))
# MCMC via Gibbs will provide output in *.txt format  
# from C routine to avoide large data problem in R
set.seed(11)
post.ar.fitpred <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,   
         data=DataFit, model="AR", coords=coords, 
         newcoords=pred.coords, newdata=DataValPred,
         scale.transform="SQRT")
print(post.ar.fitpred)
summary(post.ar.fitpred)
names(post.ar.fitpred)
# validation criteria
spT.validation(DataValPred$o8hrmax,c(post.ar.fitpred$prediction[,1]))  

#################################
## The GPP approximation models:
#################################

##
## Model fitting
##

# Read data 
data(NYdata); 

# Define the coordinates
coords<-as.matrix(unique(cbind(NYdata[,2:3])))
# Define knots
knots<-spT.grid.coords(Longitude=c(max(coords[,1]),
              min(coords[,1])),Latitude=c(max(coords[,2]),
              min(coords[,2])), by=c(4,4))

# MCMC via Gibbs using default choices
set.seed(11)
post.gpp <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,   
         data=NYdata, model="GPP", coords=coords, 
         knots.coords=knots, scale.transform="SQRT")
print(post.gpp)

# MCMC via Gibbs not using default choices
# define the time-series 
time.data<-spT.time(t.series=62,segment=1)

# hyper-parameters for the prior distributions
priors<-spT.priors(model="GPP",inv.var.prior=Gamm(2,1),
        beta.prior=Norm(0,10^4))

# initial values for the model parameters
initials<-spT.initials(model="GPP", sig2eps=0.01, 
            sig2eta=0.5, beta=NULL, phi=0.001)

# input for spatial decay
#spatial.decay<-spT.decay(distribution="FIXED", value=0.001)
spatial.decay<-spT.decay(distribution=Gamm(2,1), tuning=0.05)
#spatial.decay<-spT.decay(distribution=Unif(0.001,0.009),npoints=10)

# Iterations for the MCMC algorithms
nItr<-5000 

# MCMC via Gibbs
set.seed(11)
post.gpp <- spT.Gibbs(formula=o8hrmax~cMAXTMP+WDSP+RH, 
         data=NYdata, model="GPP", time.data=time.data, 
         coords=coords, knots.coords=knots,
         priors=priors, initials=initials, 
         nItr=nItr, nBurn=0, report=nItr, 
         tol.dist=2, distance.method="geodetic:km", 
         cov.fnc="exponential", scale.transform="SQRT", 
         spatial.decay=spatial.decay)
print(post.gpp)

# Summary and plots
summary(post.gpp)
plot(post.gpp)

# Model selection criteria
post.gpp$PMCC 

##
## Fit and spatially prediction simultaneously
##

# Read data 
s<-c(8,11,12,14,18,21,24,28)
DataFit<-spT.subset(data=NYdata, var.name=c("s.index"), s=s, reverse=TRUE) 
DataFit<-subset(DataFit, with(DataFit, !(Day \%in\% c(30, 31) & Month == 8)))
DataValPred<-spT.subset(data=NYdata, var.name=c("s.index"), s=s) 
DataValPred<-subset(DataValPred, with(DataValPred, !(Day \%in\% c(30, 31) & Month == 8)))
# Define the coordinates
coords<-as.matrix(unique(cbind(DataFit[,2:3])))
pred.coords<-as.matrix(unique(cbind(DataValPred[,2:3])))
knots<-spT.grid.coords(Longitude=c(max(coords[,1]),
              min(coords[,1])),Latitude=c(max(coords[,2]),
              min(coords[,2])), by=c(4,4))
# MCMC via Gibbs will provide output in *.txt format  
# from C routine to avoide large data problem in R
set.seed(11)
post.gpp.fitpred <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,   
         data=DataFit, model="GP", coords=coords, knots.coords=knots,
         newcoords=pred.coords, newdata=DataValPred,
         scale.transform="SQRT")
print(post.gpp.fitpred)
summary(post.gpp.fitpred)
plot(post.gpp.fitpred)
names(post.gpp.fitpred)
# validation criteria
spT.validation(DataValPred$o8hrmax,c(post.gpp.fitpred$prediction[,1]))  


######################################################
## The Truncated/Censored GP models:
######################################################

##
## Model fitting
##

data(NYdata)

# Truncation at 30 (say)
NYdata$o8hrmax[NYdata$o8hrmax<=30] <- 30

# Read data 
s<-c(8,11,12,14,18,21,24,28)
DataFit<-spT.subset(data=NYdata, var.name=c("s.index"), s=s, reverse=TRUE) 
DataFit<-subset(DataFit, with(DataFit, !(Day \%in\% c(30, 31) & Month == 8)))
DataValPred<-spT.subset(data=NYdata, var.name=c("s.index"), s=s) 
DataValPred<-subset(DataValPred, with(DataValPred, !(Day \%in\% c(30, 31) & Month == 8)))
DataValFore<-spT.subset(data=NYdata, var.name=c("s.index"), s=c(8,11,12,14,18,21,24,28)) 
DataValFore<-subset(DataValFore, with(DataValFore, (Day \%in\% c(30, 31) & Month == 8)))
DataFitFore<-spT.subset(data=NYdata, var.name=c("s.index"), s=c(8,11,12,14,18,21,24,28),
			reverse=TRUE) 
DataFitFore<-subset(DataFitFore, with(DataFitFore, (Day \%in\% c(30, 31) & Month == 8)))

#
nItr <- 5000 # number of MCMC samples for each model
nBurn <- 1000 # number of burn-in from the MCMC samples
# Truncation at 30 
# fit truncated GP model 
out <- spT.Gibbs(formula=o8hrmax~cMAXTMP+WDSP+RH,data=DataFit,
  model="truncatedGP",coords=~Longitude+Latitude,
  distance.method="geodetic:km",nItr=nItr,nBurn=nBurn,report=5,
  truncation.para = list(at = 30,lambda = 4),
  fitted.values="ORIGINAL")
#  
summary(out)
head(fitted(out))
plot(out,density=FALSE)
#
head(cbind(DataFit$o8hrmax,fitted(out)[,1]))
plot(DataFit$o8hrmax,fitted(out)[,1])
spT.validation(DataFit$o8hrmax,fitted(out)[,1])

##
## prediction (spatial)
##

pred <- predict(out,newdata=DataValPred, newcoords=~Longitude+Latitude, tol=0.05)
names(pred)
plot(DataValPred$o8hrmax,c(pred$Mean)) 
spT.validation(DataValPred$o8hrmax,c(pred$Mean)) 
#pred$prob.below.threshold

##
## forecast (temporal)
##

# unobserved locations
fore <- predict(out,newdata=DataValFore, newcoords=~Longitude+Latitude,
   type="temporal", foreStep=2, tol=0.05)
spT.validation(DataValFore$o8hrmax,c(fore$Mean)) 
plot(DataValFore$o8hrmax,c(fore$Mean)) 
#fore$prob.below.threshold

# observed locations 
fore <- predict(out,newdata=DataFitFore, newcoords=~Longitude+Latitude,
   type="temporal", foreStep=2, tol=0.05)
spT.validation(DataFitFore$o8hrmax,c(fore$Mean)) 
plot(DataFitFore$o8hrmax,c(fore$Mean)) 
#fore$prob.below.threshold


######################################################
## The Truncated/Censored GPP models:
######################################################

##
## Model fitting
##

data(NYdata)

# Define the coordinates
coords<-as.matrix(unique(cbind(NYdata[,2:3])))
# Define knots
knots<-spT.grid.coords(Longitude=c(max(coords[,1]),
              min(coords[,1])),Latitude=c(max(coords[,2]),
              min(coords[,2])), by=c(4,4))

# Truncation at 30 (say)
NYdata$o8hrmax[NYdata$o8hrmax<=30] <- 30

# Read data 
s<-c(8,11,12,14,18,21,24,28)
DataFit<-spT.subset(data=NYdata, var.name=c("s.index"), s=s, reverse=TRUE) 
DataFit<-subset(DataFit, with(DataFit, !(Day \%in\% c(30, 31) & Month == 8)))
DataValPred<-spT.subset(data=NYdata, var.name=c("s.index"), s=s) 
DataValPred<-subset(DataValPred, with(DataValPred, !(Day \%in\% c(30, 31) & Month == 8)))
DataValFore<-spT.subset(data=NYdata, var.name=c("s.index"), s=c(8,11,12,14,18,21,24,28)) 
DataValFore<-subset(DataValFore, with(DataValFore, (Day \%in\% c(30, 31) & Month == 8)))
DataFitFore<-spT.subset(data=NYdata, var.name=c("s.index"), s=c(8,11,12,14,18,21,24,28),
			reverse=TRUE) 
DataFitFore<-subset(DataFitFore, with(DataFitFore, (Day \%in\% c(30, 31) & Month == 8)))

#
nItr <- 5000 # number of MCMC samples for each model
nBurn <- 1000 # number of burn-in from the MCMC samples
# Truncation at 30 
# fit truncated GPP model 
out <- spT.Gibbs(formula=o8hrmax ~cMAXTMP+WDSP+RH,   
         data=DataFit, model="truncatedGPP",coords=~Longitude+Latitude,
         knots.coords=knots, distance.method="geodetic:km",
         nItr=nItr,nBurn=nBurn,report=5,fitted="ORIGINAL",
         truncation.para = list(at = 30,lambda = 4))
#  
summary(out)
head(fitted(out))
plot(out,density=FALSE)
#
head(cbind(DataFit$o8hrmax,fitted(out)[,1]))
plot(DataFit$o8hrmax,fitted(out)[,1])
spT.validation(DataFit$o8hrmax,fitted(out)[,1])

##
## prediction (spatial)
##

pred <- predict(out,newdata=DataValPred, newcoords=~Longitude+Latitude, tol=0.05)
names(pred)
plot(DataValPred$o8hrmax,c(pred$Mean)) 
spT.validation(DataValPred$o8hrmax,c(pred$Mean)) 
#pred$prob.below.threshold

##
## forecast (temporal)
##

# unobserved locations
fore <- predict(out,newdata=DataValFore, newcoords=~Longitude+Latitude,
   type="temporal", foreStep=2, tol=0.05)
spT.validation(DataValFore$o8hrmax,c(fore$Mean)) 
plot(DataValFore$o8hrmax,c(fore$Mean)) 
#fore$prob.below.threshold

# observed locations 
fore <- predict(out,newdata=DataFitFore, newcoords=~Longitude+Latitude,
   type="temporal", foreStep=2, tol=0.05)
spT.validation(DataFitFore$o8hrmax,c(fore$Mean)) 
plot(DataFitFore$o8hrmax,c(fore$Mean)) 
#fore$prob.below.threshold


######################################################
######################################################
##
}
}
\keyword{spT}