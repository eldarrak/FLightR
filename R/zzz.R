print.FLightR.version <- function()
{ library(help=FLightR)$info[[1]] -> version
	version <- version[pmatch("Version",version)]
	if(!is.null(version))
	{
		um <- strsplit(version," ")[[1]]
  	    version <- um[nchar(um)>0][2]
	}
	hello <- paste("This is FLightR ",version,"\n Note that for use of plotting functions relying on google maps you should get the Google maps api key\n",sep="")
	packageStartupMessage(hello)
}

.onAttach <- function(...) { 
	print.FLightR.version()
}
