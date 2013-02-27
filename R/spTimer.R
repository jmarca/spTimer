
####################################################################################


library(coda)


.onLoad <-
    function(libname, pkgname)
{
      library.dynam(pkgname, pkgname, lib.loc=libname)
}


.onAttach <-
    function(libname, pkgname)
{
	## figureout this year automatically
	this.year <- substr(as.character(Sys.Date()),1,4)
	## output to screen
        packageStartupMessage("##\n## Spatio-Temporal Bayesian Modelling using R")
        packageStartupMessage("## 2010-", this.year, 
		", License: GPL\n##", sep="")
}


####################################################################################
