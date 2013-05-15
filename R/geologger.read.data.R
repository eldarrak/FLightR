

geologger.read.data<-function( file) {
	track <- read.csv(file,header=FALSE, stringsAsFactors=F) 
	names(track)<-c('date','time','light')
	track$datetime <- paste(track$date,track$time,sep=' ') #makes a new column called datetime with date and time concatenated together with a space between
	track$gmt <- as.POSIXct(strptime(track$datetime,'%d/%m/%Y %H:%M:%S'),'UTC') #makes a new column called gmt with date and time data in a formate that R can work with
	d <- data.frame(id =1:nrow(track), gmt = track$gmt, light = track$light) 
	Data<-list(d=d)
	return(Data)
}

