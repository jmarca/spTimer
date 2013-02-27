##
## prior function for spTimer
##
spT.priors<-function(model="GP",var.prior=Gam(a = 2,b = 1), 
           beta.prior=Nor(0,10^10), rho.prior=Nor(0,10^10),
           phi.prior=Gam(a = 2, b = 1))
{
   #
   #if(!model %in% c("GP", "AR", "GPP")){
   # stop("\n# Error: correctly define model \n")
   #}
   #
   if(model=="GP"){
     #
     if(is.na(var.prior[,1])){
       prior_a<-NULL
     }
     else{
       prior_a<-var.prior[,1]
     }
     if(is.na(var.prior[,2])){
       prior_b<-NULL
     }
     else{
       prior_b<-var.prior[,2]
     }
     if(is.na(beta.prior[,1])){
       prior_mubeta<-NULL
     }
     else{
       prior_mubeta<-beta.prior[,1]
     }
     if(is.na(beta.prior[,2])){
       prior_sigbeta<-NULL
     }
     else{
       prior_sigbeta<-beta.prior[,2]
     }
     if(is.na(rho.prior[,1])){
       prior_murho<-NULL
     }
     else{
       prior_murho<-rho.prior[,1]
     }
     if(is.na(rho.prior[,2])){
       prior_sigrho<-NULL
     }
     else{
       prior_sigrho<-rho.prior[,2]
     }
     ##
     ## Priors for the GP models
     ##
      priors.gp<-function(prior_a=NULL, prior_b=NULL, prior_mubeta=NULL,  
            prior_sigbeta=NULL, prior_omu=NULL, prior_osig=NULL)
     {
      out <- NULL
      out$prior_a <- prior_a
      out$prior_b <- prior_b
      out$prior_mubeta <- prior_mubeta
      out$prior_sigbeta <- prior_sigbeta
      out$prior_omu <- prior_omu
      out$prior_osig <- prior_osig
      class(out) <- "spGP"
      out
     }
     #
       priors.gp(prior_a, prior_b, prior_mubeta, prior_sigbeta,
       prior_omu=NULL, prior_osig=NULL)
     #
   }
   #
   else if(model=="AR"){
     #
     if(is.na(var.prior[,1])){
       prior_a<-NULL
     }
     else{
       prior_a<-var.prior[,1]
     }
     if(is.na(var.prior[,2])){
       prior_b<-NULL
     }
     else{
       prior_b<-var.prior[,2]
     }
     if(is.na(beta.prior[,1])){
       prior_mubeta<-NULL
     }
     else{
       prior_mubeta<-beta.prior[,1]
     }
     if(is.na(beta.prior[,2])){
       prior_sigbeta<-NULL
     }
     else{
       prior_sigbeta<-beta.prior[,2]
     }
     if(is.na(rho.prior[,1])){
       prior_murho<-NULL
     }
     else{
       prior_murho<-rho.prior[,1]
     }
     if(is.na(rho.prior[,2])){
       prior_sigrho<-NULL
     }
     else{
       prior_sigrho<-rho.prior[,2]
     }
     ##
     ## Priors for the AR models
     ##
      priors.ar<-function(prior_a=NULL,prior_b=NULL,prior_sig=NULL)
     {
      out <- NULL
      out$prior_a <- prior_a
      out$prior_b <- prior_b
      out$prior_sig <- prior_sig
      class(out) <- "spAR"
      out
     }
     #
       priors.ar(prior_a,prior_b,prior_sigbeta)
     #
   }
   #
   else if(model=="GPP"){
     #
     if(is.na(var.prior[,1])){
       prior_a<-NULL
     }
     else{
       prior_a<-var.prior[,1]
     }
     if(is.na(var.prior[,2])){
       prior_b<-NULL
     }
     else{
       prior_b<-var.prior[,2]
     }
     if(is.na(beta.prior[,1])){
       prior_mubeta<-NULL
     }
     else{
       prior_mubeta<-beta.prior[,1]
     }
     if(is.na(beta.prior[,2])){
       prior_sigbeta<-NULL
     }
     else{
       prior_sigbeta<-beta.prior[,2]
     }
     if(is.na(rho.prior[,1])){
       prior_murho<-NULL
     }
     else{
       prior_murho<-rho.prior[,1]
     }
     if(is.na(rho.prior[,2])){
       prior_sigrho<-NULL
     }
     else{
       prior_sigrho<-rho.prior[,2]
     }
     ##
     ## Priors for the GPP models
     ##
     priors.gpp<-function(prior_a=NULL,prior_b=NULL,
          mu_beta=NULL,delta2_beta=NULL,
          mu_rho=NULL,delta2_rho=NULL,
          alpha_l=NULL,delta2_l=NULL)
     {
        out <- NULL
        out$prior_a <- prior_a
        out$prior_b <- prior_b
        out$mu_beta <- mu_beta
        out$delta2_beta <- matrix(0,length(delta2_beta),length(delta2_beta))
        diag(out$delta2_beta)<-delta2_beta
        out$mu_rho <- mu_rho
        out$delta2_rho <- delta2_rho
        out$alpha_l <- alpha_l
        out$delta2_l <- delta2_l
        class(out) <- "spGPP"
        out
     }
     #
       priors.gpp(prior_a,prior_b,mu_beta=prior_mubeta,
       delta2_beta=prior_sigbeta,mu_rho=prior_murho,
       delta2_rho=prior_sigrho,alpha_l=NULL,delta2_l=NULL)
     #
   }
   #
   else {
        stop("\n# Error: correctly define model  \n")
   }
   #
}
##
## Gamma prior
##
Gam<-function(a=NA,b=NA){
   matrix(c(a,b),1,2)
}
##
## Normal prior
##
Nor<-function(mu=NA,sig=NA){
   matrix(c(mu,sig),1,2)
}
##
## initial function for spTimer
##
spT.initials<-function(model="GP", sig2eps=0.01, 
            sig2eta=NULL, rho=NULL, beta=NULL, phi=NULL)
{
   #
   #if(!model %in% c("GP", "AR", "GPP")){
   # stop("\n# Error: correctly define model \n")
   #}
   #
   #
   if(model=="GP"){
   ##
   ## Initial values for the GPP models
   ##
    initials.gp<-function(phi=NULL,sig2eps=NULL,sig2eta=NULL,beta=NULL)
   {
      out <- NULL
      out$phi=phi; out$sig2eps=sig2eps; out$sig2eta=sig2eta; 
      out$beta=beta; 
      class(out) <- "spGP"
      out
   }
   #
      initials.gp(phi,sig2eps,sig2eta,beta)
   }
   #
   else if(model=="AR"){
   ##
   ## Initial values for the AR models
   ##
    initials.ar<-function(phi=NULL,sig2eps=NULL,sig2eta=NULL,   
            sig2l=NULL,rho=NULL,beta=NULL,mu_l=NULL)
   {
      out <- NULL
      out$phi=phi; out$sig2eps=sig2eps; out$sig2eta=sig2eta; out$sig_l0=sig2l;
      out$rho=rho; out$beta=beta; out$mu_l=mu_l; 
      class(out) <- "spAR"
      out
   }
   #
      initials.ar(phi, sig2eps, sig2eta, sig2l=NULL,
      rho,beta,mu_l=NULL)
   }
   #
   else if(model=="GPP"){
   ##
   ## Initial values for the GPP models
   ##
    initials.gpp<-function(phi=NULL, sig2eps=NULL, sig2eta=NULL, 
            sig2l=NULL, rho=NULL, beta=NULL, mu_l=NULL)
   {
      out <- NULL
      out$phi=phi; out$sig2eps=sig2eps; out$sig2eta=sig2eta; out$sig2l=sig2l;
      out$rho=rho; out$beta=beta; out$mu_l=mu_l; 
      class(out) <- "spGPP"
      out
   }
   #
      initials.gpp(phi, sig2eps, sig2eta, sig2l=NULL, 
      rho, beta, mu_l=NULL)
   }
   #
   else {
        stop("\n Error: correctly define model \n")
   }
   #
}
##
## Function for time data
##
spT.time<-function(t.series,segment=1)
{
  c(segment,0,t.series)
}
##
## Function for the spatial decay selection
##
 spT.decay<-function(type="MH", tuning=NULL, limit=NULL, segments=NULL, value=NULL)
{
   #
   if(!type %in% c("FIXED", "DISCRETE", "MH")){
    stop("\n# Error: correctly define spatial.decay \n")
   }
   #
   out<-NULL
   if(type == "FIXED"){
    if(is.null(value)){
     stop("\n# Error: correctly define discrete limits \n")
    }
    out$type<-"FIXED"
    out$value<-value
    out
   }
   else if(type == "DISCRETE"){
    if(is.null(limit)){
     stop("\n# Error: correctly define discrete limits \n")
    }
    if(is.null(segments)){
     stop("\n# Error: correctly define discrete segments \n")
    }
    out$type<-"DISCRETE"
    out$values<-seq(from=limit[1],to=limit[2],
        by=(limit[2]-limit[1])/(segments-1))
    out$segments<-segments
    out
   }
   else if(type == "MH"){
    if(is.null(tuning)){
     stop("\n# Error: correctly define tuning parameter \n")
    }
    out$type<-"MH"
    out$tuning<-tuning
    out
   }
   else{
    stop("\n# Error: correctly define spatial.decay \n")
   }
}
##
## Gibbs function for spTimer
##
#spT.Gibbs<-function(formula, ...) UseMethod("spT")
##
spT.Gibbs<-function(formula, data=parent.frame(), model="GP",
         time.data=NULL, coords, knots.coords, 
         newcoords=NULL, newdata=NULL, 
         priors=NULL, initials=NULL, nItr=13000, nBurn=3000, report=1, 
         tol.dist=2, distance.method="geodetic:km", 
         cov.fnc="exponential", scale.transform="NONE", 
         spatial.decay=spT.decay(type="MH",tuning=0.1),
         annual.aggrn="NONE")
{
   ##
   sum.stat<-function(x)
   {
     model<-x$model
     nItr<-x$iterations-x$nBurn
     nBurn<-0
     cat("\n")
     cat("# Model:", model, "\n")
     if(is.null(model)==TRUE){
      stop("\n# Error: need to define the model")
     }
     else if(model=="AR"){
        if(nItr <= nBurn){
        stop("\n# Error: iteration (",nItr,") is less than or equal to nBurn (",nBurn,") \n")
        }
     r<-dim(x$mu_lp)[[1]]
     p<-dim(x$betap)[[1]]
     if(cov.fnc=="matern"){
     para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$rhop[(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]),t(x$nup[(nBurn+1):nItr]))
     para<-spT.Summary.Stat(para)
     dimnames(para)[[1]][1:(1+p+4)]<-c(dimnames(x$X)[[2]],"rho","sig2eps","sig2eta", "phi", "nu")
     }
     else {
     para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$rhop[(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]))
     para<-spT.Summary.Stat(para)
     dimnames(para)[[1]][1:(1+p+3)]<-c(dimnames(x$X)[[2]],"rho","sig2eps","sig2eta", "phi")
     }
     #round(para,4)
     para
     }
     else if(model == "GPP"){
        if(nItr <= nBurn){
        stop("\n# Error: iteration (",nItr,") is less than or equal to nBurn (",nBurn,") \n")
        }
     r<-x$r
     p<-x$p
     if(cov.fnc=="matern"){
     para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$rhop[(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]),t(x$nup[(nBurn+1):nItr]))
     para<-spT.Summary.Stat(para)
     dimnames(para)[[1]][1:(1+p+2+1+1)]<-c(dimnames(x$X)[[2]],"rho","sig2eps","sig2eta","phi","nu")
     }
     else {
     para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$rhop[(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]))
     para<-spT.Summary.Stat(para)
     dimnames(para)[[1]][1:(1+p+2+1)]<-c(dimnames(x$X)[[2]],"rho","sig2eps","sig2eta","phi")
     }
     #round(para,4)
     para
     }
     else if(model == "GP"){
        if(nItr <= nBurn){
        stop("\n# Error: iteration (",nItr,") is less than or equal to nBurn (",nBurn,") \n")
        }
     r<-x$r
     p<-x$p
     if(cov.fnc=="matern"){
     para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]),t(x$nup[(nBurn+1):nItr]))
     para<-spT.Summary.Stat(para)
     dimnames(para)[[1]]<-c(dimnames(x$X)[[2]],"sig2eps","sig2eta","phi","nu")
     }
     else {
     para<-rbind((x$betap[,(nBurn+1):nItr]),
              t(x$sig2ep[(nBurn+1):nItr]),
              t(x$sig2etap[(nBurn+1):nItr]),
              t(x$phip[(nBurn+1):nItr]))
     para<-spT.Summary.Stat(para)
     dimnames(para)[[1]]<-c(dimnames(x$X)[[2]],"sig2eps","sig2eta","phi")
     }
     #round(para,4)
     para
     }
     else{

     }
   }
   # check coords: built-in the data
   if(missing(coords) || is.null(coords)){
     if(!missing(data)){
       sstr <- names(data)
       coords <- data[,sstr[sstr %in% c("Longitude","Latitude","xcoords","ycoords")]]
       if(dim(coords)[[2]]==0){
         stop("\n Error: need to specify the coords using argument coords or\n through the supplied dataset, see help in spT.Gibbs")
       }
       else if(dim(coords)[[2]]>2){
         stop("\n Error: check the dataset for supplied coordinates. \n It contains more than two columns of coordinate information.")
       }
       else{ 
         coords <- as.matrix(unique(coords))
       }
     } 
   }
   ##
   if(is.null(newcoords) | is.null(newdata)) {
   # 
   if (missing(formula)) {
          stop("\n# Error: formula must be specified \n")
   }
   if (class(formula) != "formula") {
          stop("\n# Error: equation must be in formula-class \n ...")
   }
   #
   if(!distance.method %in% c("geodetic:km", "geodetic:mile", "euclidean",
      "maximum", "manhattan", "canberra")){
    stop("\n# Error: correctly define distance.method \n")
   }
   #
   if(!cov.fnc %in% c("exponential", "gaussian", "spherical",
      "matern")){
    stop("\n# Error: correctly define cov.fnc \n")
   }
   #
   if(!scale.transform %in% c("NONE", "SQRT", "LOG")){
    stop("\n# Error: correctly define scale.transform \n")
   }
   #
   if(!spatial.decay$type %in% c("FIXED", "DISCRETE", "MH")){
    stop("\n# Error: correctly define spatial.decay \n")
   }
   #
   if(!model %in% c("GP", "AR", "GPP")){
    stop("\n# Error: correctly define model \n")
   }
   #
   if(report < 1){
    report<-1
   }
   #
   if(model=="GP"){
      cat("\n Output: GP models \n")
      if(class(priors) != "spGP" & class(priors) != "NULL"){
        stop("\n# Error: correctly define the GP models for function spT.priors.")
      }
      if(class(initials) != "spGP" & class(initials) != "NULL"){
        stop("\n# Error: correctly define the GP models for function spT.initials.")
      }
      if(!missing(knots.coords)){
        stop("\n# Error: please check options for the GP model and/or knots.coords \n")
      }
      out<- spGP.Gibbs(formula=formula, data=data, time.data=time.data, coords=coords, 
            priors=priors, initials=initials, nItr=nItr, nBurn=nBurn, report=report, tol.dist=tol.dist,   
            distance.method=distance.method, cov.fnc=cov.fnc, scale.transform=scale.transform, spatial.decay=spatial.decay,
            X.out=TRUE, Y.out=TRUE)
      out$combined.fit.pred<-FALSE
      out$model<-model
      out$parameter<-sum.stat(out)
      out$data<-data
      class(out)<-"spT"
      out
   }
   #
   else if(model=="AR"){
      cat("\n Output: AR models \n")
      if(class(priors) != "spAR" & class(priors) != "NULL"){
        stop("\n# Error: correctly define the AR models for function spT.priors.")
      }
      if(class(initials) != "spAR" & class(initials) != "NULL"){
        stop("\n# Error: correctly define the AR models for function spT.initials.")
      }
      if(!missing(knots.coords)){
        stop("\n# Error: please check options for the AR model and/or knots.coords \n")
      }
      out<- spAR.Gibbs(formula=formula, data=data, time.data=time.data, coords=coords, 
            priors=priors, initials=initials, nItr=nItr, nBurn=nBurn, report=report, tol.dist=tol.dist,   
            distance.method=distance.method, cov.fnc=cov.fnc, scale.transform=scale.transform, spatial.decay=spatial.decay,
            X.out=TRUE, Y.out=TRUE)
      out$combined.fit.pred<-FALSE
      out$model<-model
      out$parameter<-sum.stat(out)
      out$data<-data
      class(out)<-"spT"
      out
   }
   #
   else if(model=="GPP"){
      cat("\n Output: GPP approximation models \n")
      if(missing(knots.coords)){
        stop("\n# Error: define knots.coords \n")
      }
      if(class(priors) != "spGPP" & class(priors) != "NULL"){
        stop("\n# Error: correctly define the GPP models for function spT.priors.")
      }
      if(class(initials) != "spGPP" & class(initials) != "NULL"){
        stop("\n# Error: correctly define the GPP models for function spT.initials.")
      }
      out<- spGPP.Gibbs(formula=formula, data=data, time.data=time.data, 
            knots.coords=knots.coords, coords=coords, priors=priors, initials=initials, 
            nItr=nItr, nBurn=nBurn, report=report, tol.dist=tol.dist, 
            distance.method=distance.method, cov.fnc=cov.fnc,
            scale.transform=scale.transform, spatial.decay=spatial.decay, 
            X.out=TRUE, Y.out=TRUE)
      out$combined.fit.pred<-FALSE
      out$model<-model
      out$parameter<-sum.stat(out)
      out$data<-data
      class(out)<-"spT"
      out
   }
   #
   else{
    stop("\n# Error: correctly define model \n")
   }
   #
   } # end of loop
   ##
   ## for both fit and pred
   ##
   else if (!is.null(newcoords) & !is.null(newdata)) {
     #
     if(!annual.aggrn %in% c("NONE", "ave", "an4th", "w126")){
      stop("\n# Error: correctly define annual.aggrn \n")
     }
     if(is.null(newcoords)){
      stop("\n# Error: correctly define newcoords \n")
     }
     if(is.null(newdata)){
      stop("\n# Error: correctly define newdata \n")
     }
     #
    out<-spT.fit.pred(formula=formula, data=data, model=model, time.data=time.data, coords=coords, 
            knots.coords=knots.coords, pred.coords=newcoords, priors=priors,
            initials=initials, pred.data=newdata, nItr=nItr, nBurn=nBurn, report=report,
            tol.dist=tol.dist, distance.method=distance.method, cov.fnc=cov.fnc, scale.transform=scale.transform,
            spatial.decay=spatial.decay, annual.aggrn=annual.aggrn)
    class(out)<-"spT"
    out
   } 
   else {
    stop("\n# Error: correctly define newcoords and newdata \n")
   }
}
##
## Prediction function for spTimer
##
##
print.spTprd<-function(x, ...) {
    cat("--------------------------------------"); cat('\n');
    cat("Prediction with Models:", x$model, "\n")
    cat("Covariance function:", x$cov.fnc, "\n")
    cat("Distance method:", x$distance.method, "\n")
    cat("Computation time: ",x$computation.time, "\n")
    cat("--------------------------------------"); cat('\n');
}
##
spT.prediction<-function(nBurn=0, pred.data=NULL, pred.coords,
     posteriors, tol.dist=2, Summary=TRUE)
{
   #
   if (missing(posteriors)) {
       stop("Error: need to provide the posterior MCMC samples.")
   }
   #
   if(posteriors$combined.fit.pred == TRUE){
    stop("\n#\n# Termination: Prediction results are already obtained with the fitted models \n#\n")
   }
   #
   model<-posteriors$model
   #
   if(!model %in% c("GP", "AR", "GPP")){
    stop("\n# Error: model is not correctly defined \n")
   }
   #
   else if(model=="GP"){
      cat("\n Prediction: GP models \n")
      out<-spGP.prediction(nBurn=nBurn, pred.data=pred.data, pred.coords=pred.coords, 
           posteriors=posteriors, tol.dist=tol.dist, Summary=Summary)
      out$model<-model
      class(out)<-"spTprd"
      out
   }
   else if(model=="AR"){
      cat("\n Prediction: AR models \n")
      out<-spAR.prediction(nBurn=nBurn, pred.data=pred.data, pred.coords=pred.coords, 
           posteriors=posteriors, tol.dist=tol.dist, Summary=Summary)
      out$model<-model
      class(out)<-"spTprd"
      out
   }
   else if(model=="GPP"){
      cat("\n Prediction: GPP models \n")
      out<-spGPP.prediction(nBurn=nBurn, pred.data=pred.data, pred.coords=pred.coords, 
           posteriors=posteriors, Summary=Summary)
      out$model<-model
      class(out)<-"spTprd"
      out
   }
   else {
    stop("\n# Error: model is not correctly defined \n")
   }
}
##
## Forecast function for spTimer
##
print.spTfore<-function(x, ...) {
    cat("--------------------------------------"); cat('\n');
    cat("Forecast with Models:", x$model, "\n")
    cat("Covariance function:", x$cov.fnc, "\n")
    cat("Distance method:", x$distance.method, "\n")
    cat("Computation time: ",x$computation.time, "\n")
    cat("--------------------------------------"); cat('\n');
}
##
spT.forecast<-function(nBurn=0, K=1, fore.data=NULL, fore.coords, 
            posteriors, pred.samples.ar=NULL, Summary=TRUE)
{
  #
   if (missing(posteriors)) {
       stop("Error: need to provide the posterior MCMC samples.")
   }
   #
   model<-posteriors$model
   #
   if(!model %in% c("GP", "AR", "GPP")){
    stop("\n# Error: model is not correctly defined \n")
   }
   #
   else if(model=="GP"){
      cat("\n Forecast: GP models \n")
      out<-spGP.forecast(nBurn=nBurn, K=K, fore.data=fore.data,
           fore.coords=fore.coords, posteriors=posteriors, 
           Summary=Summary)
      out$model<-model
      class(out)<-"spTfore"
      out
   }
   #
   else if(model=="AR"){
      if(missing(pred.samples.ar)){
        stop("Error: define predAR by value or by NULL.")
      }
      cat("\n Forecast: AR models \n")
      out<-spAR.forecast(nBurn=nBurn, K=K, fore.data=fore.data, fore.coords=fore.coords, 
           posteriors=posteriors, pred.samples.ar=pred.samples.ar, Summary=Summary)
      out$model<-model
      class(out)<-"spTfore"
      out
   }
   #
   else if(model=="GPP"){
      cat("\n Forecast: GPP models \n")
      out<-spGPP.forecast(nBurn=nBurn, K=K, fore.data=fore.data, fore.coords=fore.coords, 
           posteriors=posteriors, Summary=Summary)
      out$model<-model
      class(out)<-"spTfore"
      out
   }
   #
   else {
    stop("\n# Error: model is not correctly defined \n")
   }
}
##
## use of predict function
##
predict.spT<-function(object, newdata=NULL, newcoords, foreStep=NULL, type="spatial", nBurn=0, 
          tol.dist=2, predAR=NULL, ...)
{
     if(type=="spatial"){
        if(is.null(newcoords)){
          stop("Error: define newcoords.")
        }
     out<-spT.prediction(nBurn=nBurn, pred.data=newdata, pred.coords=newcoords,
          posteriors=object, tol.dist=tol.dist, Summary=TRUE)
     out$type<-"spatial"
     class(out)<-"spTpred" 
     out
     }
     else if(type=="temporal"){
        if(is.null(newcoords)){
          stop("Error: define newcoords.")
        }
        if(is.null(foreStep)){
          stop("Error: define foreStep.")
        }
     out<-spT.forecast(nBurn=nBurn, K=foreStep, fore.data=newdata, fore.coords=newcoords, 
            posteriors=object, pred.samples.ar=predAR, Summary=TRUE)
     out$type<-"temporal"
     class(out)<-"spTpred"
     out
     }
     else{
       stop("\n# Error: prediction type is not correctly defined \n")
     }
}
#
print.spTpred<-function(x, ...) {
    if(x$type=="spatial"){
    cat("--------------------------------------"); cat('\n');
    cat("Spatial prediction with Model:", x$model, "\n")
    cat("Covariance function:", x$cov.fnc, "\n")
    cat("Distance method:", x$distance.method, "\n")
    cat("Computation time: ",x$computation.time, "\n")
    cat("--------------------------------------"); cat('\n');
    }
    else if(x$type=="temporal"){
    cat("--------------------------------------"); cat('\n');
    cat("Temporal prediction/forecast with Model:", x$model, "\n")
    cat("Covariance function:", x$cov.fnc, "\n")
    cat("Distance method:", x$distance.method, "\n")
    cat("Computation time: ",x$computation.time, "\n")
    cat("--------------------------------------"); cat('\n');
    }
    else{
    }
}
##
## Fit and prediction function for AR and GPP
##
 spT.fit.pred<-function(formula, data=parent.frame(), model="AR", 
            time.data, coords, knots.coords, pred.coords, priors,
            initials, pred.data=NULL, nItr, nBurn=0, report=1,
            tol.dist=2, distance.method="geodetic:km",  
            cov.fnc="exponential", scale.transform="NONE",
            spatial.decay,annual.aggrn="NONE")
{
   # 
   if (missing(formula)) {
          stop("\n# Error: formula must be specified \n")
   }
   if (class(formula) != "formula") {
          stop("\n# Error: equation must be in formula-class \n ...")
   }
   #
   if(!distance.method %in% c("geodetic:km", "geodetic:mile", "euclidean",
      "maximum", "manhattan", "canberra")){
    stop("\n# Error: correctly define distance.method \n")
   }
   #
   if(!cov.fnc %in% c("exponential", "gaussian", "spherical",
      "matern")){
    stop("\n# Error: correctly define cov.fnc \n")
   }
   #
   if(!scale.transform %in% c("NONE", "SQRT", "LOG")){
    stop("\n# Error: correctly define scale.transform \n")
   }
   #
   if(!spatial.decay$type %in% c("FIXED", "DISCRETE", "MH")){
    stop("\n# Error: correctly define spatial.decay \n")
   }
   #
   if(!model %in% c("GP", "AR", "GPP")){
    stop("\n# Error: correctly define model \n")
   }
   #
   if(model=="GP"){
      #cat("\n Currently unavailable \n")
      if(class(priors) != "spGP" & class(priors) != "NULL"){
        stop("\n# Error: correctly define the GP models for function spT.priors.")
      }
      if(class(initials) != "spGP" & class(initials) != "NULL"){
        stop("\n# Error: correctly define the GP models for function spT.initials.")
      }
      out <- spGP.MCMC.Pred(formula=formula, data=data,
             time.data=time.data, coords=coords, 
             pred.coords=pred.coords, priors=priors,
             initials=initials, pred.data=pred.data, 
             nItr=nItr, nBurn=nBurn, report=report, tol.dist=tol.dist,
             distance.method=distance.method, cov.fnc=cov.fnc, 
             scale.transform=scale.transform, 
             spatial.decay=spatial.decay,
             annual.aggregation=annual.aggrn)
      out$combined.fit.pred<-TRUE
      out$model<-model
      out
   }
   #
   else if(model=="AR"){
      cat("\n Output: AR models \n")
      if(class(priors) != "spAR" & class(priors) != "NULL"){
        stop("\n# Error: correctly define the AR models for function spT.priors.")
      }
      if(class(initials) != "spAR" & class(initials) != "NULL"){
        stop("\n# Error: correctly define the AR models for function spT.initials.")
      }
      out <- spAR.MCMC.Pred(formula=formula, data=data,
             time.data=time.data, coords=coords, 
             pred.coords=pred.coords, priors=priors,
             initials=initials, pred.data=pred.data, 
             nItr=nItr, nBurn=nBurn, report=report, tol.dist=tol.dist,
             distance.method=distance.method, cov.fnc=cov.fnc, 
             scale.transform=scale.transform, 
             spatial.decay=spatial.decay,
             annual.aggregation=annual.aggrn)
      out$combined.fit.pred<-TRUE
      out$model<-model
      out
   }
   #
   else if(model=="GPP"){
      cat("\n Output: GPP approximation models \n")
      if(missing(knots.coords)){
        stop("\n# Error: define knots.coords \n")
      }
      if(class(priors) != "spGPP" & class(priors) != "NULL"){
        stop("\n# Error: correctly define the GPP models for function spT.priors.")
      }
      if(class(initials) != "spGPP" & class(initials) != "NULL"){
        stop("\n# Error: correctly define the GPP models for function spT.initials.")
      }
      out <- spGPP.MCMC.Pred(formula=formula, data=data,
             time.data=time.data, knots.coords=knots.coords, 
             coords=coords, pred.coords=pred.coords, priors=priors,
             initials=initials, pred.data=pred.data, nItr=nItr, nBurn=nBurn, 
             report=report, tol.dist=tol.dist,
             distance.method=distance.method, cov.fnc=cov.fnc, 
             scale.transform=scale.transform, 
             spatial.decay=spatial.decay,
             annual.aggregation=annual.aggrn) 
      out$combined.fit.pred<-TRUE
      out$model<-model
      out
   }
   #
   else{
    stop("\n# Error: correctly define model for both fit.pred fnc. \n#    currently GP model is not available.")
   }
   #
}
##
##
#setMethod(f="plot", signature="spGP", 
#   definition=spT.MCMC.plot
#)
#setMethod(f="plot", signature="spGPP", 
#   definition=spT.MCMC.plot
#)
#setMethod(f="summary", signature(object="spAR"), 
#   definition=spT.MCMC.stat
#)
##
##

##
##
##
