print.FLightR.version <- function()
{ library(help=FLightR)$info[[1]] -> version
	version <- version[pmatch("Version",version)]
	if(!is.null(version))
	{
		um <- strsplit(version," ")[[1]]
  	    version <- um[nchar(um)>0][2]
	}
	hello <- paste("This is FLightR ",version,"\n NB: All functions that started from plot.**.**() are renamed to plot_**_** in order to pass CRAN check \n",sep="")
	packageStartupMessage(hello)
}

.onAttach <- function(...) { 
	print.FLightR.version()
}
