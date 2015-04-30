
convert.lux.to.tags<-function(file, log=T, log.light.borders=c(1,10)) {
	# the function takes the current (2015)
	# .lux format and converts it to .csv format that
	# TAGS service accepts

	Dat<-read.csv(file, skip=20, sep="\t", stringsAsFactors =F)
	names(Dat)<-c("datetime", "light")
	Dat$datetime<-as.POSIXct(Dat$datetime, tz="UTC", format="%d/%m/%Y %H:%M:%S")

	Dat_new<-Dat
	Dat_new$datetime<-format(Dat_new$datetime, format="%Y-%m-%d %H:%M:%S")
	
	if (log) {
	Dat_new$light<-log(Dat_new$light)
	
	if (!is.null(log.light.borders)) {
	Dat_new$light[Dat_new$light<log.light.borders[1]]<-log.light.borders[1]
	Dat_new$light[Dat_new$light>log.light.borders[2]]<-log.light.borders[2]
	}
	}
	# now I need to save as csv...
	write.csv(Dat_new, file=paste(unlist(strsplit(file, ".lux")), "csv", sep="."), quote = F, row.names=F) 
	cat("Success!\n")
	cat("file", paste(unlist(strsplit(file, ".lux")), "csv", sep="."), "\nwas saved to", getwd(), "\n")
	return(NULL)
}