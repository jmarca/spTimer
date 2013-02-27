##
## Summary Statistics another built in
##
spT.MCMC.stat<-function(x, nBurn=0)
{
  options(warn=-1)
  if(class(x) != "spT"){
    stop("\n# Error: provide valid posterior output \n")
  }
  model<-x$model
  nItr<-x$iterations
  if(!is.null(nBurn)){
    if(nBurn == 0){
    nBurn<-0
    nItr<-nItr-x$nBurn
    }
    if(nBurn != 0){ 
    nBurn<-nBurn
    nItr<-nItr-x$nBurn
    }
  }
  if((nBurn+x$nBurn) >= (nItr+x$nBurn)){
   cat("# Number of Iterations:            ", nItr+x$nBurn, "\n")
   cat("# Number of Burn-in (fitted model):", x$nBurn, "\n")
   cat("# More Burn-in:                    ", nBurn, "\n")
   cat("# Total Number of Burn-in:         ", nBurn+x$nBurn, "\n")
   stop("\n# Error: iteration (",nItr+x$nBurn,") is less than or equal to total burn-in (",nBurn+x$nBurn,") \n")
  }
  cat("\n")
  cat("# Model:", model, "\n")
  if(is.null(model)==TRUE){
   stop("\n# Error: need to define the model")
  }
  else if(model=="AR"){
  cat("\n")
   cat("# Number of Iterations:   ", nItr+x$nBurn, "\n")
   cat("# Total number of Burn-in:", nBurn+x$nBurn, "\n")
  cat("\n")
     if(nItr <= nBurn){
     stop("\n# Error: iteration (",nItr,") is less than or equal to nBurn (",nBurn,") \n")
     }
  r<-dim(x$mu_lp)[[1]]
  p<-dim(x$betap)[[1]]
  #para<-rbind((x$betap[,(nBurn+1):nItr]),
  #            t(x$rhop[(nBurn+1):nItr]),
  #            t(x$sig2ep[(nBurn+1):nItr]),
  #            t(x$sig2etap[(nBurn+1):nItr]),
  #            (x$sig2lp[,(nBurn+1):nItr]),
  #            (x$mu_lp[,(nBurn+1):nItr]),
  #            t(x$phip[(nBurn+1):nItr]))
  #para<-spT.Summary.Stat(para)
  #dimnames(para)[[1]][1:(1+p+2)]<-c(dimnames(x$X)[[2]],"rho",
  #   "sig2eps","sig2eta")
  #dimnames(para)[[1]][(1+p+2+1):(1+p+2+r)] <- paste("sig2:", 1:r)  
  #dimnames(para)[[1]][(1+p+2+r+1):(1+p+2+r+r)] <- paste("mu:", 1:r)  
  #dimnames(para)[[1]][(1+p+2+r+r+1)] <-c("phi")
  #
  if(x$cov.fnc=="matern"){
  para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$rhop[(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]),t(x$nup[(nBurn+1):nItr]))
  }
  else{
  para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$rhop[(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]))
  }
  #
  para<-spT.Summary.Stat(para)
  dimnames(para)[[1]][1:(1+p+3)]<-c(dimnames(x$X)[[2]],"rho",
     "sig2eps","sig2eta","phi")
  round(para,4)
  }
  else if(model == "GPP"){
  cat("\n")
   cat("# Number of Iterations:   ", nItr+x$nBurn, "\n")
   cat("# Total number of Burn-in:", nBurn+x$nBurn, "\n")
  cat("\n")
     if(nItr <= nBurn){
     stop("\n# Error: iteration (",nItr,") is less than or equal to nBurn (",nBurn,") \n")
     }
  r<-x$r
  p<-x$p
  #
  if(x$cov.fnc=="matern"){
  para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$rhop[(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]),t(x$nup[(nBurn+1):nItr]))
  }
  else{
  para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$rhop[(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]))
  }
  #
  para<-spT.Summary.Stat(para)
  dimnames(para)[[1]][1:(1+p+2+1)]<-c(dimnames(x$X)[[2]],"rho",
     "sig2eps","sig2eta","phi")
  round(para,4)
  }
  else if(model == "GP"){
  cat("\n")
   cat("# Number of Iterations:   ", nItr+x$nBurn, "\n")
   cat("# Total number of Burn-in:", nBurn+x$nBurn, "\n")
  cat("\n")
     if(nItr <= nBurn){
     stop("\n# Error: iteration (",nItr,") is less than or equal to nBurn (",nBurn,") \n")
     }
  r<-x$r
  p<-x$p
  #
  if(x$cov.fnc=="matern"){
  para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]),t(x$nup[(nBurn+1):nItr]))
  }
  else{
  para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]))
  }
  #
  para<-spT.Summary.Stat(para)
  dimnames(para)[[1]]<-c(dimnames(x$X)[[2]],
     "sig2eps","sig2eta","phi")
  round(para,4)
  }
  else{

  }
}
##
## MCMC plot another built in 
##
spT.MCMC.plot<-function(x, nBurn=0, ACF="FALSE", PARTIAL.acf="FALSE")
{
  options(warn=-1)
  if(class(x) != "spT"){
    stop("\n# Error: provide valid posterior output \n")
  }
  model<-x$model
  nItr<-x$iterations
  if(!is.null(nBurn)){
    if(nBurn == 0){
    nBurn<-0
    nItr<-nItr-x$nBurn
    }
    if(nBurn != 0){ 
    nBurn<-nBurn
    nItr<-nItr-x$nBurn
    }
  }
  if((nBurn+x$nBurn) >= (nItr+x$nBurn)){
   cat("# Number of Iterations:            ", nItr+x$nBurn, "\n")
   cat("# Number of Burn-in (fitted model):", x$nBurn, "\n")
   cat("# More Burn-in:                    ", nBurn, "\n")
   cat("# Total Number of Burn-in:         ", nBurn+x$nBurn, "\n")
   stop("\n# Error: iteration (",nItr+x$nBurn,") is less than or equal to total burn-in (",nBurn+x$nBurn,") \n")
  }
  cat("\n")
  cat("# Model:", model, "\n")
  if(is.null(model)==TRUE){
   stop("\n# Error: need to define the model")
  }
  else if(model=="AR"){
  cat("\n")
   cat("# Number of Iterations:   ", nItr+x$nBurn, "\n")
   cat("# Total number of Burn-in:", nBurn+x$nBurn, "\n")
  cat("\n")
     if(nItr <= nBurn){
     stop("\n# Error: iteration (",nItr,") is less than or equal to nBurn (",nBurn,") \n")
     }
  r<-x$r
  p<-x$p
  #para<-rbind((x$betap[,(nBurn+1):nItr]),
  #            t(x$rhop[(nBurn+1):nItr]),
  #            t(x$sig2ep[(nBurn+1):nItr]),
  #            t(x$sig2etap[(nBurn+1):nItr]),
  #            (x$sig2lp[,(nBurn+1):nItr]),
  #            (x$mu_lp[,(nBurn+1):nItr]),
  #            t(x$phip[(nBurn+1):nItr]))
  #dimnames(para)[[1]][1:(1+p+2+r+r+1)]<-c(dimnames(x$X)[[2]],"rho",
  #   "sig2eps","sig2eta",paste("sig2:", 1:r),paste("mu:", 1:r),"phi")
  #initial_val<-c(x$initials[[7]],x$initials[6],x$initials[2],
  #           x$initials[3],x$initials[5],x$initials[4],x$initials[1])
  #
  if(x$cov.fnc=="matern"){
  para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$rhop[(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]),t(x$nup[(nBurn+1):nItr]))
  dimnames(para)[[1]][1:(1+p+2+1+1)]<-c(dimnames(x$X)[[2]],"rho",
     "sig2eps","sig2eta","phi","nu")
  }
  else{
  para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$rhop[(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]))
  dimnames(para)[[1]][1:(1+p+2+1)]<-c(dimnames(x$X)[[2]],"rho",
     "sig2eps","sig2eta","phi")
  }
  #
  #
  }
  else if(model == "GPP"){
  cat("\n")
   cat("# Number of Iterations:   ", nItr+x$nBurn, "\n")
   cat("# Total number of Burn-in:", nBurn+x$nBurn, "\n")
  cat("\n")
     if(nItr <= nBurn){
     stop("\n# Error: iteration (",nItr,") is less than or equal to nBurn (",nBurn,") \n")
     }
  r<-dim(x$mu_lp)[[1]]
  p<-dim(x$betap)[[1]]
  #para<-rbind((x$betap[,(nBurn+1):nItr]),
  #            t(x$rhop[(nBurn+1):nItr]),
  #            t(x$sig2ep[(nBurn+1):nItr]),
  #            t(x$sig2etap[(nBurn+1):nItr]),
  #            t(x$phip[(nBurn+1):nItr]))
  #dimnames(para)[[1]][1:(1+p+2+1)]<-c(dimnames(x$X)[[2]],"rho",
  #   "sig2eps","sig2eta","phi")
  #initial_val<-c(x$initials[[7]],x$initials[4],x$initials[2],
  #           x$initials[3],x$initials[1])
  #
  if(x$cov.fnc=="matern"){
  para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$rhop[(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]),t(x$nup[(nBurn+1):nItr]))
  dimnames(para)[[1]][1:(1+p+2+1+1)]<-c(dimnames(x$X)[[2]],"rho",
     "sig2eps","sig2eta","phi","nu")
  }
  else{
  para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$rhop[(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]))
  dimnames(para)[[1]][1:(1+p+2+1)]<-c(dimnames(x$X)[[2]],"rho",
     "sig2eps","sig2eta","phi")
  }
  #
  }
  else if(model == "GP"){
  cat("\n")
   cat("# Number of Iterations:   ", nItr+x$nBurn, "\n")
   cat("# Total number of Burn-in:", nBurn+x$nBurn, "\n")
  cat("\n")
     if(nItr <= nBurn){
     stop("\n# Error: iteration (",nItr,") is less than or equal to nBurn (",nBurn,") \n")
     }
  r<-x$r
  p<-x$p
  #para<-rbind((x$betap[,(nBurn+1):nItr]),
  #            t(x$sig2ep[(nBurn+1):nItr]),
  #            t(x$sig2etap[(nBurn+1):nItr]),
  #            t(x$phip[(nBurn+1):nItr]))
  #dimnames(para)[[1]]<-c(dimnames(x$X)[[2]],
  #   "sig2eps","sig2eta","phi")
  #initial_val<-c(x$initials[[4]],x$initials[2],
  #             x$initials[3],x$initials[1])
  #
  if(x$cov.fnc=="matern"){
  para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]),t(x$nup[(nBurn+1):nItr]))
  dimnames(para)[[1]][1:(1+p+2+1+1)]<-c(dimnames(x$X)[[2]],
     "sig2eps","sig2eta","phi","nu")
  }
  else{
  para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]))
  dimnames(para)[[1]][1:(1+p+2+1)]<-c(dimnames(x$X)[[2]],
     "sig2eps","sig2eta","phi")
  }
  #
  #
  }
  else{

  }
  #
  ##
   x11()
   for(i in 1:dim(para)[[1]]){
    MCMC.plot.obj(para[i,],name=c(dimnames(para)[[1]][i]),ACF=ACF,PARTIAL.acf=PARTIAL.acf)
    par(ask=TRUE)
   }
}
##
## segment plot for upper and lower limit
##  
spT.segment.plot<-function(obs, est, up, low, limit=NULL){
  #
  tmp<-cbind(obs,est,up,low)
  tmp<-na.omit(tmp)
  if(is.null(limit)==TRUE){
  plot(tmp[,1],tmp[,2],xlab="Observations",ylab="Predictions", pch="*",
  xlim=c(min(c(tmp),na.rm=TRUE),max(c(tmp),na.rm=TRUE)),
  ylim=c(min(c(tmp),na.rm=TRUE),max(c(tmp),na.rm=TRUE)))
  }
  #
  else{
  plot(tmp[,1],tmp[,2],xlab="Observations",ylab="Predictions",
   xlim=c(limit[1],limit[2]),ylim=c(limit[1],limit[2]),pch="*")
  }
  #
  segments(tmp[,1],tmp[,2],tmp[,1],tmp[,3])
  segments(tmp[,1],tmp[,2],tmp[,1],tmp[,4])
  #  
}
##
## hit and false alarm function for forecast
##
 spT.hit.false<-function(obs,fore,tol){
   #
   tmp<-cbind(obs,fore)
   tmp<-na.omit(tmp)
   #
   c11<-tmp[tmp[,1]<=tol & tmp[,2]<=tol,]
   c11<-length(c11)/2
   c12<-tmp[tmp[,1]<=tol & tmp[,2]>tol,]
   c12<-length(c12)/2  
   c21<-tmp[tmp[,1]>tol & tmp[,2]<=tol,]
   c21<-length(c21)/2  
   c22<-tmp[tmp[,1]>tol & tmp[,2]>tol,]
   c22<-length(c22)/2
   mat<-matrix(c(c11,c21,c12,c22),2,2)
   dimnames(mat)[[1]]<-c(paste("[Obs:<=",tol,"]"),paste("[Obs:> ",tol,"]")) 
   dimnames(mat)[[2]]<-c(paste("[Forecast:<=",tol,"]"),paste("[Forecast:>",tol,"]")) 
   POD<-round(mat[1,1]/sum(diag(mat)),4)
   FAR<-round(mat[1,2]/sum(mat[1,]),4)
   HAR<-round(sum(diag(mat))/sum(mat),4)
   top<-2*(mat[1,1]*mat[2,2]-mat[1,2]*mat[2,1])
   bot<-mat[1,2]^2+mat[2,1]^2+2*mat[1,1]*mat[2,2]+(mat[1,2]+mat[2,1])*sum(diag(mat))
   S<-round(top/bot,4)
   x<-list(False.Alarm=FAR,Hit.Rate=HAR,Probability.of.Detection=POD,
      Heidke.Skill=S,cross.table=mat,tolerance.limit=tol) 
   x
}
##
## For data split
##
 spT.data.selection<-function(data, random = TRUE, num.rs = NULL, 
                     s = NULL, reverse=FALSE) 
{
#
# This function is to select and deduct the sites used to
# fit or valid the model
# Input: 	data
#		s = the site numbers to be selected, e.g., c(4,7,10)
#		num.rs = the number of sites to be selected , e.g., 3
  if(reverse==FALSE){
	if(random != TRUE){
		num.rs<-NULL
		num.rs<-length(s)
		dat<-NULL
			for(i in 1:num.rs){
			dat<-rbind(dat,data[data[,1]==s[i], ])
			}	
		dat	
	}	
	else{
		s<-NULL
		a.s<-unique(data[,1])
		s<-sort(sample(a.s, num.rs))
		cat('Randomly selected sites =', s, '\n')
		dat<-NULL
			for(i in 1:num.rs){
			dat<-rbind(dat,data[data[,1]==s[i], ])
			}	
		dat


	}
  }
  else{
	if(random != TRUE){
		num.rs<-NULL
		num.rs<-length(s)
			for(i in 1:num.rs){
			data<-data[data[,1] != s[i], ]
			}	
		data	
	}	
	else{
		s<-NULL
		a.s<-unique(data[,1])
		s<-sort(sample(a.s, num.rs))
		cat('Randomly deducted sites =', s, '\n')
			for(i in 1:num.rs){
			data<-data[data[,1] != s[i], ]
			}	
		data
	}
  }
 }
##
## grid coordinates
##
 spT.grid.coords<-function(Longitude=c(max,min),Latitude=c(max,min),by=c(NA,NA))
{
      max.lon <- Longitude[1]
      min.lon <- Longitude[2]
      max.lat <- Latitude[1]
      min.lat <- Latitude[2]
      if(is.na(by[1]) || is.na(by[2])){
        stop('Error: need to specify grid dimension, n x m')
      }
      knots.lon <- seq(max.lon,min.lon,length=by[1])
      knots.lat <- seq(max.lat,min.lat,length=by[2])
      #knots.coords <- cbind(rep(knots.lon,by),c(outer(rep(1,by),knots.lat)))
      knots.coords <- cbind(c(outer(rep(1,by[2]),knots.lon)),rep(knots.lat,by[1]))
      as.matrix(knots.coords)
}
##
## Percentage Coverage
##
 spT.pCOVER<-function(z=NULL,zup=NULL,zlow=NULL,zsample=NULL,level=95)
{
    if(is.null(z)){
      stop('\n # Provide observed values\n')
    }
    if(!is.null(zsample)){
      z<-c(z)
      low<-(1-level/100)/2
      up<-1-low
      x<-apply(zsample,1,quantile,prob=c(low,up))
      x<-rbind(x,z)
      x<-t(x)
      x<-na.exclude(x)
      y<-x[x[, 1] <= x[, 3] & x[, 3] <= x[, 2], ]
      round(length(y[,1])/length(x[,1])*100,2)
    }
    else if(is.null(zsample)){
      z<-as.matrix(z)
      zlow<-as.matrix(zlow)
      zup<-as.matrix(zup)
      x<-cbind(zlow,z,zup)
      x<-na.exclude(x)
      y<- x[x[,1] <= x[,2] & x[,2] <= x[,3],]
      round(length(y[,1])/length(x[,1])*100,2)
    }
}
##
## Geodetic distance in K.M. or Miles
##
 spT.geodist<-function(Lon, Lat, KM=TRUE){
 #
 # This function is for geodistance using C/C++
 #
   n<-length(Lat)
   if(KM==TRUE){	
   ds<-.C("GeoDist_km", as.integer(n),as.double(Lat),
    as.double(Lon),out=matrix(double(n*n),n,n))$out
   }
   else {
   ds<-.C("GeoDist_miles", as.integer(n),as.double(Lat),
    as.double(Lon),out=matrix(double(n*n),n,n))$out
   }
   #dis<-matrix(ds,length(Lat),length(Lat))
   #return(dis)
   ds
 }
##
## Geodetic distance: another approach: two points
##
 spT.geo.dist <- function(point1, point2){
 #
 # The following program computes the distance on the surface 
 # of the earth between two points point1 and point2. 
 # Both the points are of the form (Longitude, Latitude)
 #
   R <- 6371
   p1rad <- point1 * pi/180
   p2rad <- point2 * pi/180
   d <- sin(p1rad[2])*sin(p2rad[2])+cos(p1rad[2])*cos(p2rad[2])*cos(abs(p1rad[1]-p2rad[1]))	
   d <- acos(d)
   R*d
 }
##
## Geodetic distance: using points 1:(lon, lat) 2:(lon, lat)
##
spT.geo_dist <- function(points)
{
        point1 <- c(points[1], points[2])
        point2 <- c(points[3], points[4])
        #The following program computes the distance on the surface 
        #The argument points is of the form (Long, Lat, Long, Lat)
        dist <- 0
        if(sum((point1 - point2)^2) > 1e-05) {
                R <- 6371
                p1rad <- (point1 * pi)/180
                p2rad <- (point2 * pi)/180
                d <- sin(p1rad[2]) * sin(p2rad[2]) + cos(p1rad[2]) * cos(p2rad[
                        2]) * cos(abs(p1rad[1] - p2rad[1]))
                d <- acos(d)
                dist <- R * d
        }
        dist
}
##
## To keep the long lat positions based on tol.dist
##
 spT.keep.morethan.dist <- function(coords, tol.dist=100)  
{
 #
 # coords must have two columns named long and lat 
 #
   a <- as.data.frame(coords)
   names(a) <- c("long","lat")
   n <- nrow(coords)
   c1 <- rep(1:n, each=n)
   c2 <- rep(1:n, n)
   b <- matrix(NA, nrow=n, ncol=n)
   w <- as.vector(upper.tri(b))
   bigmat <- matrix(0, nrow=n*n, ncol=7)
   bigmat[, 1] <- c1
   bigmat[, 2] <- c2
   bigmat[, 3] <- a$long[c1]
   bigmat[, 4] <- a$lat[c1]
   bigmat[, 5] <- a$long[c2]
   bigmat[, 6] <- a$lat[c2]
   ubmat <- bigmat[w, ]
   ubmat[,7] <- as.vector(apply(ubmat[,3:6], 1, spT.geo_dist)) 

   v <- ubmat[,7]
   w <- ubmat[v<tol.dist, ]

   z <- unique(w[,1])
   a <- coords[-z,  ]
   a
}
##
##
##
# Sigma_inv<-function(n, phi, D){
#
# Exponential covariance function
# phi = decay parameter
# D = distance matrix
#
#	Si_eta<-.C("MInv", as.double(exp(-phi*D)), inv=double(n*n), 
#		as.integer(n), dt=as.double(1))$inv
#	Si_eta
# }
##
##
##
# Mean<-function(x){
#
# Mean function for missing data vector 'x'
#
#	x<-mean(x,na.rm=TRUE)
#	x
# }
##
##
##
# Median<-function(x){
#
# Median function for missing data vector 'x'
#
#	x<-median(x,na.rm=TRUE)
#	x
# }
##
## This function extracts the data using formula
##
 Formula.matrix <-  function(formula, data, na.rm=FALSE)
{
      if(na.rm == FALSE){
           mt <- terms(formula, data=data)
           if(missing(data)) data <- sys.frame(sys.parent())
      
           mf <- model.frame(mt, data=data, na.action=na.pass)
	
      # null model support
           X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)
           X <- as.matrix(X)         # X matrix
           xvars <- dimnames(X)[[2]] # X variable names
  
           Y <- as.matrix(model.response(mf, "numeric")) # Y matrix
      return(list(Y, X, xvars))
      }
      else{
           mt <- terms(formula, data=data)
           if(missing(data)) data <- sys.frame(sys.parent())
      
           mf <- model.frame(mt, data=data)
	
      # null model support
           X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)
           X <- as.matrix(X)         # X matrix
           xvars <- dimnames(X)[[2]] # X variable names
  
           Y <- as.matrix(model.response(mf, "numeric")) # Y matrix
      return(list(Y, X, xvars))
      }

}
##
## Summary Statistics
##
 spT.Summary.Stat <- function(y) 
{

     ## define Upper and lower limit function
       up.low.limit<-function(y,limit)
       {
         #
          y<-sort(y)
          y<-y[limit]
          y
         #
       }
       #
         if(is.vector(y)==TRUE){
            y <- matrix(y[!is.na(y)])
         }
         else{
            y <- t(y)
         } 
         N <- length(y[1,  ])
         nItr <- length(y[, 1])
         z <- matrix(nrow = N, ncol = 5)
         dimnames(z) <- list(dimnames(y)[[2]], c("Mean","Median","SD","Low2.5p","Up97.5p"))
        #
         if (nItr < 40) {
           stop("\n##\n# Error: number of samples must be >= 40\n##")
         }
         z[, 1] <- apply(y,2,mean)
         z[, 2] <- apply(y,2,median)
         z[, 3] <- apply(y,2,sd)
         #z[, 4] <- apply(y,2,quantile,0.025)
         #z[, 5] <- apply(y,2,quantile,0.975)
        nl <- as.integer(nItr * 0.025)
        nu <- as.integer(nItr * 0.975)
         z[, 4] <- apply(y,2,up.low.limit,limit=nl)
         z[, 5] <- apply(y,2,up.low.limit,limit=nu)
         
        #
        # out <- .C('sum_stat',as.integer(nItr),as.integer(N),
        #    as.double(y),as.integer(1),Mean=double(N),Median=double(N),
        #    Var=double(N),Low=double(N),Up=double(N))[5:9]
        #
        # z[, 1] <- out$Mean
        # z[, 2] <- out$Median
        # z[, 3] <- sqrt(out$Var)
        # z[, 4] <- out$Low
        # z[, 5] <- out$Up
        #
        #nl <- as.integer(nItr * 0.025)
        #nu <- as.integer(nItr * 0.975)
        #for(i in 1:N) {
        #        tmp <- sort(y[, i])
        #        z[i, 1] <- mean(tmp)
        #        z[i, 2] <- median(tmp)
        #        z[i, 3] <- sqrt(var(tmp))
        #        z[i, 4] <- tmp[nl]
        #        z[i, 5] <- tmp[nu]
        #}
        #
        as.data.frame(z)
        #z
}
##
## Predictive Model Choice Criteria
##
 PMCC<-function(z=NULL, z.mean=NULL, z.sd=NULL, z.samples=NULL)
{
    #
    # Predictive Model Choice Criteria
    #
    if(is.null(z)){
      stop("Error: need to provide z values.")
    }
    #
    if(is.null(z.mean) | is.null(z.sd)){
      if(!is.null(z.mean)){
       stop("Error: need to provide z.sd value.")
      }   
      if(!is.null(z.sd)){
       stop("Error: need to provide z.mean value.")
      }   
      if(is.null(z.samples)){
       stop("Error: need to provide z.samples value.")
      }
    }
    #
    if(!is.null(z.samples)){
     if ( !is.matrix(z.samples)) {
         stop("Error: z.samples must be a (N x nItr) matrix")
     }
     if (dim(z.samples)[1] != length(z)) {
         stop("Error: observations in z.samples in each iteration must be equal to length of z")
     }
     if ( dim(z.samples)[2] < 40) {
         stop("Error: samples are too small to obtain summary statistics")
     }
     sum.stat = matrix(NA,length(c(z)),6)
     sum.stat[,1:5] = as.matrix(spT.Summary.Stat(z.samples))
     sum.stat[,6] = c(z)
     sum.stat = sum.stat[!is.na(sum.stat[,6]),]
     goodness.of.fit = round(sum((sum.stat[,1]-sum.stat[,6])^2),2)
     penalty = round(sum(sum.stat[,3]^2),2)
     pmcc =  round(goodness.of.fit + penalty,2)
     out = NULL
     out$pmcc = pmcc; 
     out$goodness.of.fit = goodness.of.fit
     out$penalty = penalty
     #class(out) <- "PMCC"
     out
    }
    else{
     if(is.null(z.mean) | is.null(z.sd)){
       stop("Error: need to provide z.mean and/or z.sd values.")
     }
     if(length(c(z)) != length(c(z.mean))){
       stop("Error: z and z.mean should be in same length.")
     }
     if(length(c(z)) != length(c(z.sd))){
       stop("Error: z and z.sd should be in same length.")
     }
     sum.stat = matrix(NA,length(c(z)),3)
     sum.stat[,1] = c(z)
     sum.stat[,2] = c(z.mean)
     sum.stat[,3] = c(z.sd)
     sum.stat = sum.stat[!is.na(sum.stat[,1]),]
     goodness.of.fit = round(sum((sum.stat[,1]-sum.stat[,2])^2),2)
     penalty = round(sum(sum.stat[,3]^2),2)
     pmcc =  round(goodness.of.fit + penalty,2)
     out = NULL
     out$pmcc = pmcc; 
     out$goodness.of.fit = goodness.of.fit
     out$penalty = penalty
     #class(out) <- "PMCC"
     out
    }
}

##
## To identify the Integer Number
##

#is.wholenumber<-function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

##
## To check sites that has <tol km of distances
## 
 spT.check.locations<-function(fit.locations, pred.locations,
              method="geodetic:km", tol=5){
  #
      #
      if(!method %in% c("geodetic:km", "geodetic:mile", "euclidean",
        "maximum", "manhattan", "canberra")){
        stop("\n# Error: correctly define distance.method \n")
      }
      #
           coords.all <- rbind(fit.locations,pred.locations)
           tn.fitsites <- length(fit.locations[, 1])
           nfit.sites <- 1:tn.fitsites
           tn.predsites <- length(coords.all[, 1]) - tn.fitsites
           npred.sites <- (tn.fitsites + 1):(length(coords.all[, 1]))
      #
      if(method=="geodetic:km"){
         coords.D <- as.matrix(spT.geodist(Lon=coords.all[,1],Lat=coords.all[,2], KM=TRUE))
      }
      else if(method=="geodetic:mile"){
         coords.D <- as.matrix(spT.geodist(Lon=coords.all[,1],Lat=coords.all[,2], KM=FALSE))
      }
      else{
       coords.D <- as.matrix(dist(coords.all, method, diag = T, upper = T))
      }  
           coords.D[is.na(coords.D)]<-0
      #
           diag(coords.D)<-NA
      # 
           fdmis<-coords.D[nfit.sites, npred.sites]
     if(is.matrix(fdmis)==TRUE){
           fdmis<-cbind(c(t(fdmis)),1:dim(fdmis)[[2]],sort(rep(1:dim(fdmis)[[1]],dim(fdmis)[[2]]))) # add pred sites and fitted sites
           fdmis<-fdmis[fdmis[,1] < tol,]
      #
           if(!is.na(fdmis[1])==TRUE){
            cat("#\n# Tolerance Limit:", paste(tol))
            cat("\n# There are some Prediction locations very close to the Fitted locations.\n#\n")
            fdmis<-matrix(fdmis,(length(fdmis)/3),3) 
            for(i in 1:dim(fdmis)[[1]]){
            print(paste("Distance:", round(fdmis[i,1],2)," Predicted location:",fdmis[i,2]," Fitted location:", fdmis[i,3],""))
            }
            cat("#\n# Romove the locations and run again. \n#\n")
            dimnames(fdmis)[[2]]<-c('distance','pred_location','fit_location')
            fdmis
           }
           else{
            cat("#\n# Tolerance Limit (unit):", paste(tol))
            cat("\n# Fitted and Predicted location distances are alright \n#\n")  
           }
      }
      else{
           fdmis<-cbind(c(fdmis),1:length(fdmis)) # 
           fdmis<-fdmis[fdmis[,1] < tol,]
           if(!is.na(fdmis[1])==TRUE){
            cat("#\n# Tolerance Limit:", paste(tol))
            cat("\n# There are some Prediction locations very close to the Fitted locations.\n#\n")
            fdmis<-matrix(fdmis) 
            for(i in 1:dim(fdmis)[[1]]){
            print(paste("Distance:", round(fdmis[i,1],2)," Predicted location:",1," Fitted location:", fdmis[i,2],""))
            }
            cat("#\n# Romove the locations and run again. \n#\n")
            dimnames(fdmis)[[2]]<-c('distance','fit_location')
            fdmis
           }
           else{
            cat("#\n# Tolerance Limit (unit):", paste(tol))
            cat("\n# Fitted and Predicted location distances are alright \n#\n")  
           }
       }
}
##
## To check sites inside codes
## 
spT.check.sites.inside<-function(coords, method){
      #
      #
      if(!method %in% c("geodetic:km", "geodetic:mile", "euclidean",
        "maximum", "manhattan", "canberra")){
        stop("\n# Error: correctly define distance.method \n")
      }
      #
      #
      if(method=="geodetic:km"){
          fdm<- as.matrix(spT.geodist(Lon=coords[,1],Lat=coords[,2], KM=TRUE))
      }
      else if(method=="geodetic:mile"){
          fdm<- as.matrix(spT.geodist(Lon=coords[,1],Lat=coords[,2], KM=FALSE))
      }
      else{
           fdm<- as.matrix(dist(coords, method, diag = TRUE, upper = TRUE))
      } 
      #
           diag(fdm)<-NA
      # 
           fdm<-cbind(c(fdm),1:dim(fdm)[[2]],sort(rep(1:dim(fdm)[[1]],dim(fdm)[[2]]))) 
      #
           fdm<-fdm[!is.na(fdm[,1]),]
      #
           tol <- 0.01
           fdmis<-fdm[fdm[,1] < tol,]
      #
           if(!is.na(fdmis[1])==TRUE){
            cat("\n# There are some locations very close ( < 0.01 unit) to each other.\n#\n")
            fdmis<-matrix(fdmis,(length(fdmis)/3),3) 
            for(i in 1:dim(fdmis)[[1]]){
            print(paste("Distance (unit):", round(fdmis[i,1],2)," site:",fdmis[i,2]," site:", fdmis[i,3],""))
            }
            dimnames(fdmis)[[2]]<-c('dist_km','pred_site','fit_site')
            fdmis
            cat("#\n# Romove the sites and run again. \n#\n")
            stop("Error: Termination.")
           }
       #
           tol <- 0.5
           fdmis<-fdm[fdm[,1] < tol,]
      #
           if(!is.na(fdmis[1])==TRUE){
            cat("#\n# Warnings: There are some locations very close ( < 0.5 unit) to each other.\n#\n")
           }
      #
}
##
##
#check.sites.inside.pred<-spT.check.locations
##
## convert seconds into min. hour. and day
##
 fnc.time<-function(t)
{
     #
     if(t < 60){
      t <- round(t,2)
      tt <- paste(t," - Sec.")
      cat(paste("##\n# Elapsed time:",t,"Sec.\n##\n"))
     } 
     #
     if(t < (60*60) && t >= 60){
      t1 <- as.integer(t/60)
      t <- round(t-t1*60,2) 
      tt <- paste(t1," - Mins.",t," - Sec.")
      cat(paste("##\n# Elapsed time:",t1,"Min.",t,"Sec.\n##\n"))
     }
     #
     if(t < (60*60*24) && t >= (60*60)){
      t2 <- as.integer(t/(60*60))
      t <- t-t2*60*60
      t1 <- as.integer(t/60)
      t <- round(t-t1*60,2) 
      tt <- paste(t2," - Hour/s.",t1," - Mins.",t," - Sec.")
      cat(paste("##\n# Elapsed time:",t2,"Hour/s.",t1,"Min.",t,"Sec.\n##\n"))
     }
     #
     if(t >= (60*60*24)){
      t3 <- as.integer(t/(60*60*24))
      t <- t-t3*60*60*24
      t2 <- as.integer(t/(60*60))
      t <- t-t2*60*60
      t1 <- as.integer(t/60)
      t <- round(t-t1*60,2)
      tt <- paste(t3," - Day/s.",t2," - Hour/s.",t1," - Mins.",t," - Sec.")
      cat(paste("##\n# Elapsed time:",t3,"Day/s.",t2,"Hour/s.",t1,"Mins.",t,"Sec.\n##\n"))
     }
     #
     tt
}
##
## validation criteria
##
spT.validation <- function(z, zhat)
{
 ##
 ## Validation Mean Squared Error (VMSE) 
 ##
 VMSE <- function(z, zhat)
 {
       z<-as.matrix(z)
       zhat<-as.matrix(zhat)
       x <- c(z-zhat)
       u <- x[!is.na(x)]
       round(sum(u^2)/length(u), 4)
 }
 ##
 ## Root Mean Squared Error (RMSE) 
 ##
 RMSE<-function (z,zhat) 
 {
       z<-as.matrix(z)
       zhat<-as.matrix(zhat)
       x <- c(z-zhat)
       u <- x[!is.na(x)]
       round(sqrt(sum(u^2)/length(u)), 4)
 }
 ##
 ## Mean Absolute Error (MAE) 
 ##
 MAE<-function (z,zhat) 
 {
    z<-as.matrix(z)
    zhat<-as.matrix(zhat)
    x <- abs(c(zhat-z))
    u <- x[!is.na(x)]
    round(sum(u)/length(u), 4)
 }
 ##
 ## Mean Absolute Percentage Error (MAPE) 
 ##
 MAPE<-function (z,zhat) 
 {
    z<-as.matrix(z)
    zhat<-as.matrix(zhat)
    x <- abs(c(zhat-z))/z
    u <- x[!is.na(x)]
    u <- u[!is.infinite(u)]
    round(sum(u)/length(u)*100, 4)
 }
 ##
 ## Bias (BIAS) 
 ##
 BIAS<-function (z,zhat) 
 {
    z<-as.matrix(z)
    zhat<-as.matrix(zhat)
    x <- c(zhat-z)
    u <- x[!is.na(x)]
    round(sum(u)/length(u), 4)
 }
 ##
 ## Relative Bias (rBIAS) 
 ##
 rBIAS<-function (z,zhat) 
 {
    z<-as.matrix(z)
    zhat<-as.matrix(zhat)
    x <- c(zhat-z)
    u <- x[!is.na(x)]
    round(sum(u)/(length(u)*mean(z,na.rm=TRUE)), 4)
 }
 ##
 ## Relative Mean Separation (rMSEP) 
 ##
 rMSEP<-function (z,zhat) 
 {
    z<-as.matrix(z)
    zhat<-as.matrix(zhat)
    x <- c(zhat-z)
    u <- x[!is.na(x)]
    y <- c(mean(zhat,na.rm=TRUE)-z)
    v <- y[!is.na(y)]
    round(sum(u^2)/sum(v^2), 4)
 }
 ##
 cat("##\n Mean Squared Error (MSE) \n Root Mean Squared Error (RMSE) \n Mean Absolute Error (MAE) \n Mean Absolute Percentage Error (MAPE) \n Bias (BIAS) \n Relative Bias (rBIAS) \n Relative Mean Separation (rMSEP)\n##\n") 
 ##
   out<-NULL
   out$VMSE<-VMSE(z, zhat)
   out$RMSE<-RMSE(z, zhat)
   out$MAE<-MAE(z, zhat)
   #out$MAPE<-MAPE(z, zhat)
   #out$BIAS<-BIAS(z, zhat)
   out$rBIAS<-rBIAS(z, zhat)
   out$rMSEP<-rMSEP(z, zhat)
   unlist(out)
}
##
## MCMC plots for individual values (trace, density, acf)
##
MCMC.plot.obj<-function(post_val, nItr, nBurn=0, 
                name=c('....'), ACF="FALSE", PARTIAL.acf="FALSE") 
  {
      nBurn=nBurn+1; x<-post_val;
      #if(missing(initial_val)){ 
      #   y<-0
      #}
      #else{
      #   y<-initial_val
      #} 
      if(missing(nItr)){ 
         nItr<-length(post_val)
      }
      if(ACF==FALSE & PARTIAL.acf==FALSE){
        #windows()
        #x11()
        par(mfrow=c(1,2))
        #plot(nBurn:nItr, x[nBurn:nItr], xlab = "Iterations", 
        #   ylab = paste("Values of  (", name, ")", sep=''), type = "l",
        #   ylim=c(min(y,x[nBurn:nItr]),max(y,x[nBurn:nItr])),main='Trace plot')
        #abline(h = y, lty = 2,col="blue")
        #plot(density(x[nBurn:nItr]),main='Density plot',xlim=c(min(y,x[nBurn:nItr]),max(y,x[nBurn:nItr])))
        #abline(v = y, lty = 2,col="blue")
        plot(nBurn:nItr, x[nBurn:nItr], xlab = "Iterations", 
           ylab = paste("Values of  (", name, ")", sep=''), type = "l",
           ylim=c(min(x[nBurn:nItr]),max(x[nBurn:nItr])),main='Trace plot')
        plot(density(x[nBurn:nItr]),main='Density plot',xlim=c(min(x[nBurn:nItr]),max(x[nBurn:nItr])))
     }	
     else if(ACF==TRUE & PARTIAL.acf==TRUE){
      #windows()
      #x11() 
      par(mfrow=c(1,4))
      plot(nBurn:nItr, x[nBurn:nItr], xlab = "Iterations", 
        ylab = paste("Values of  (", name, ")", sep=''), type = "l", 
        ylim=c(min(x[nBurn:nItr]), max(x[nBurn:nItr])),main='Trace plot')
      #abline(h = y, lty = 2,col="blue")
      plot(density(x[nBurn:nItr]),main='Density plot',xlim=c(min(x[nBurn:nItr]), max(x[nBurn:nItr])))
      #abline(v = y, lty = 2,col="blue")
      acf(x[nBurn:nItr],main='ACF plot',col="red")
      pacf(x[nBurn:nItr],main='Partial ACF plot',col="red")
     }
     else if(ACF=="TRUE" & PARTIAL.acf=="FALSE"){
      #x11() 
      par(mfrow=c(1,3))
      plot(nBurn:nItr, x[nBurn:nItr], xlab = "Iterations", 
        ylab = paste("Values of  (", name, ")", sep=''), type = "l", 
        ylim=c(min(x[nBurn:nItr]), max(x[nBurn:nItr])),main='Trace plot')
      #abline(h = y, lty = 2, col="blue")
      plot(density(x[nBurn:nItr]),main='Density plot',xlim=c(min(x[nBurn:nItr]), max(x[nBurn:nItr])))
      #abline(v = y, lty = 2, col="blue")
      acf(x[nBurn:nItr],main='ACF plot',col="red")
     }
     else {

     }
  }
##
## Multiple imputation for initial values using Amelia with m=1
##
#imputation.z<-function(x){
  #
  # x is the r*T by n matrix
  #
#  x<-amelia(x,m=1)$imputation[[1]]
#  x
#}
##
## To use in coda package
##
as.mcmc.spT<-function(x, ...){

    model <- x$model
    if (is.null(model) == TRUE) {
        stop("\n# Error: need to define the model")
    }
    else if (model == "AR") {
        r <- x$r
        p <- x$p
        #para <- rbind((x$betap), t(x$rhop), t(x$sig2ep), t(x$sig2etap), (x$sig2lp), (x$mu_lp), 
        #        t(x$phip))
        #dimnames(para)[[1]] <- c(dimnames(x$X)[[2]],"rho", "sig2eps", "sig2eta",
        #        paste("sig2:",1:r),paste("mu:", 1:r),"phi")
      #          
      if(x$cov.fnc=="matern"){
        para <- rbind((x$betap), t(x$rhop), t(x$sig2ep), t(x$sig2etap), t(x$phip), t(x$nup))
        dimnames(para)[[1]] <- c(dimnames(x$X)[[2]],"rho", "sig2eps", "sig2eta", "phi", "nu")
      }
      else {
        para <- rbind((x$betap), t(x$rhop), t(x$sig2ep), t(x$sig2etap), t(x$phip))
        dimnames(para)[[1]] <- c(dimnames(x$X)[[2]],"rho", "sig2eps", "sig2eta", "phi")
      }
      #
        para<-t(para)
        para<-mcmc(para)
        para
    }
    else if (model == "GPP") {
        r <- x$r
        p <- x$p
      #          
      if(x$cov.fnc=="matern"){
        para <- rbind((x$betap), t(x$rhop), t(x$sig2ep), t(x$sig2etap), t(x$phip), t(x$nup))
        dimnames(para)[[1]] <- c(dimnames(x$X)[[2]],"rho", "sig2eps", "sig2eta", "phi", "nu")
      }
      else {
        para <- rbind((x$betap), t(x$rhop), t(x$sig2ep), t(x$sig2etap), t(x$phip))
        dimnames(para)[[1]] <- c(dimnames(x$X)[[2]],"rho", "sig2eps", "sig2eta", "phi")
      }
      #
        para<-t(para)
        para<-mcmc(para)
        para
    }
    else if (model == "GP") {
        r <- x$r
        p <- x$p
      #          
      if(x$cov.fnc=="matern"){
        para <- rbind((x$betap), t(x$sig2ep), t(x$sig2etap), t(x$phip), t(x$nup))
        dimnames(para)[[1]] <- c(dimnames(x$X)[[2]], "sig2eps", "sig2eta", "phi", "nu")
      }
      else {
        para <- rbind((x$betap), t(x$sig2ep), t(x$sig2etap), t(x$phip))
        dimnames(para)[[1]] <- c(dimnames(x$X)[[2]], "sig2eps", "sig2eta", "phi")
      }
      #
        para<-t(para)
        para<-mcmc(para)
        para
    }
    else {
    }
}
##
## use of model.frame
##
model.frame.spT<-function(formula, ...){
#   tmp<-cbind(formula$Y,formula$X[,-1])
#   dimnames(tmp)[[2]]<-dimnames(attr(terms(formula$call),"factors"))[[1]]
#   tmp
   if(formula$combined.fit.pred==TRUE){
      stop("\n# Error: not useful for output with combined fit and predict \n")
    }
   else{
      model.frame(formula$call,formula$data,na.action=na.pass)
   }
}
##
## use of model.matrix
##
model.matrix.spT<-function(object, ...){
   if(object$combined.fit.pred==TRUE){
      stop("\n# Error: not useful for output with combined fit and predict \n")
    }
   else{
      Formula.matrix(object$call,object$data)[[2]]
   }
}
##
## use of summary
##
summary.spT<-function(object, pack="spTimer", ...){
   if(pack=="coda"){
    if(object$combined.fit.pred==TRUE){
      stop("\n# Error: coda package is not useful for output with combined fit and predict \n")
    }
    else{
     cat("\n#### MCMC summary statistics using coda package ####\n")
     tmp<-as.mcmc(object)
     summary(tmp)
    }
   }
   else{
   cat("Parameters:\n")
   print(object$parameter); #cat("\n");
   }
}
##
## use of plot
##
plot.spT<-function(x, residuals=FALSE, ...){
   if(as.logical(residuals)==FALSE){
     if(x$combined.fit.pred==TRUE){
      #cat("\n# Error: not useful to get MCMC plots for output with combined fit and predict \n#        use residuals=TRUE \n")
      #cat("# ")
      cat("\n## Only residual plots are available for output with combined fit and predict option \n\n")
      plot(x$fitted[,1],residuals(x),ylab="Residuals",xlab="Fitted values");abline(h=0,lty=2);title("Residuals vs Fitted")
      par(ask=TRUE)
      qqnorm(residuals(x));qqline(residuals(x),lty=2)
     }
     else{
      tmp<-as.mcmc(x)
      plot(tmp)
     } 
   }
   else{
   plot(x$fitted[,1],residuals(x),ylab="Residuals",xlab="Fitted values");abline(h=0,lty=2);title("Residuals vs Fitted")
   par(ask=TRUE)
   qqnorm(residuals(x));qqline(residuals(x),lty=2)
   }
}
##
## use of coefficients
##
coef.spT<-function(object, ...){
   t(object$parameter)[1,]
}
##
## use of residuals
##
residuals.spT<-function(object, ...){
   #if(object$combined.fit.pred==TRUE){
   #   stop("\n# Error: not useful for output with combined fit and predict")
   #}
   #else{
     if(object$scale.transform=="NONE"){
     tmp<-object$Y-object$fitted[,1]
     tmp
     }
     else if(object$scale.transform=="SQRT"){
     tmp<-sqrt(object$Y)-object$fitted[,1]
     tmp
     }
     else if(object$scale.transform=="LOG"){
     tmp<-log(object$Y)-object$fitted[,1]
     tmp
     }
     else{
     }
   #}
}
##
## use of formula
##
formula.spT<-function(x, ...){
  x$call
}
##
## use of terms
##
terms.spT<-function(x, ...){
  terms(x$call)
}
##
## use of display
##


##
## use of offset
##
#offset.spT<-function(object){
#  Call<-object$call
  #Coefficients<-coef(object)
  #return(Call,Coefficients)
#}
##
##
##