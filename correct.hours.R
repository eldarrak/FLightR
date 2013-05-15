# this function will correct time values to make them not crossing 0.
correct.hours<-function(datetime) {
	# this function is supposed to correct hours by adding some amount of them.
   hours <- as.numeric(format(datetime,"%H"))+as.numeric(format(datetime,"%M"))/60
	cor <- rep(NA, 24)
	for(i in 0:23){
		cor[i+1] <- max(abs((c(hours[1],hours)+i)%%24 - 
		            (c(hours,hours[length(hours)])+i)%%24),na.rm=T)
	}
	hours <- (hours + (which.min(round(cor,2)))-1)%%24
	return(hours)
}
