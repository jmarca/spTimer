
####################################################################################


library(coda)
library(forecast)

.onLoad <-
    function(libname, pkgname)
{
      library.dynam(pkgname, pkgname, lib.loc=libname)
}


.onAttach <-
    function(libname, pkgname)
{
        ## figureout the version automatically
        library(help=spTimer)$info[[1]] -> version
        version <- version[pmatch("Version",version)]
        um <- strsplit(version," ")[[1]]
        version <- um[nchar(um)>0][2]
	## figureout this year automatically
	this.year <- substr(as.character(Sys.Date()),1,4)
	## output to screen
        packageStartupMessage("##\n## Spatio-Temporal Bayesian Modelling using R")
        packageStartupMessage("## 2010-", this.year, 
		", License: GPL", sep="")
        packageStartupMessage("## Version: ", version," \n##")
}


####################################################################################
