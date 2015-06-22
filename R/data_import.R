#######################
# here are the functions that read data


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


read.tags.light.twilight<-function(lig.raw, start.date=NULL, end.date=NULL) {
	lig.raw$datetime<-as.POSIXct(lig.raw$datetime, tz="UTC", format="%Y-%m-%dT%T")

	###
	### I also want to exclude the last days...
	
	if (!is.null(start.date)) lig.raw<-lig.raw[as.numeric(lig.raw$datetime) < as.numeric(as.POSIXct(start.date , tz="UTC")),]
	if (!is.null(end.date)) lig.raw<-lig.raw[as.numeric(lig.raw$datetime) < as.numeric(as.POSIXct(end.date , tz="UTC")),]

	## and also I want to exclude the interpolated and excluded afterwards points
		
	lig<-lig.raw[-which(lig.raw$interp==T),]

##########################################################################
##Convert the datetime/light fields into the format that FLightR works on##
##########################################################################

	lig.new<-data.frame(
	format(lig$datetime, "%d/%m/%Y"),
	format(lig$datetime, "%H:%M:%S"),
 	lig$light)

########################################################
##Save that date/light file so it can be read in later##
########################################################
	tmpname<-tempfile(fileext = ".csv")
	write.table(lig.new, file=tmpname, sep="," , row.names = F, col.names = F,  quote = F)
	Data<-geologger.read.data(file=tmpname)
	unlink(x=tmpname)
#########################################################
##Use TAGS output to define twilight periods################
#########################################################
# 

Filtered_tw <- lig.raw[(lig.raw$twilight)>0 & !lig.raw$excluded, c("datetime", "twilight", "light")]
names(Filtered_tw)[2]="type"
Filtered_tw <- Filtered_tw[!duplicated(Filtered_tw$datetime),]
Filtered_tw <- Filtered_tw[order(Filtered_tw[,1]),]


# now I want to pair data and twilights..		  
Filtered_tw$id<-0
Data$d$type<-0
Data$d<-rbind(Data$d, data.frame(id=Filtered_tw$id, gmt= Filtered_tw$datetime, light=Filtered_tw$light, type=Filtered_tw$type))

All.p<-Data$d[order(Data$d$gmt),]
rownames(All.p)<-1:nrow(All.p)

#############################################################

# now I want to pair data and twilights..		  
Filtered_tw$light<-approx(x=Data$d$gmt, y=Data$d$light, xout=Filtered_tw$datetime)$y
Filtered_tw$id<-0
Data$d$type<-0
Data$d<-rbind(Data$d, data.frame(id=Filtered_tw$id, gmt= Filtered_tw$datetime, light=Filtered_tw$light, type=Filtered_tw$type))

Filtered_tw$excluded=0

All.p<-Data$d[order(Data$d$gmt),]

#All.p<-All.p[!duplicated(All.p[,2:3], fromLast=T),]
rownames(All.p)<-1:nrow(All.p)
Res<-list(Data=All.p, twilights=Filtered_tw)
return(Res)
}


process.geolight.output<-function(datetime, light, gl_twl) {

#filter <- loessFilter(gl_twl[,1],gl_twl[,2],gl_twl[,3],k=10)
Filtered_tw <- data.frame(datetime=as.POSIXct(c(gl_twl$tFirst,gl_twl$tSecond),"UTC"),type=c(gl_twl$type,ifelse(gl_twl$type==1,2,1)))

Filtered_tw <- Filtered_tw[!duplicated(Filtered_tw$datetime),]
Filtered_tw <- Filtered_tw[order(Filtered_tw[,1]),]

# now I want to pair data and twilights..		  
Filtered_tw$light<-approx(x=datetime, y=light, xout=Filtered_tw$datetime)$y
Filtered_tw$id<-0
Data<-data.frame(id=1:length(datetime), gmt=datetime, light=light, type=0)

Data<-rbind(Data, data.frame(id=Filtered_tw$id, gmt= Filtered_tw$datetime, light=Filtered_tw$light, type=Filtered_tw$type))

All.p<-Data[order(Data$gmt),]
#All.p<-All.p[!duplicated(All.p[,2:3], fromLast=T),]
rownames(All.p)<-1:nrow(All.p)
Filtered_tw$excluded=0
FlightR.data=list(Data=All.p, twilights=Filtered_tw)
return(FlightR.data)

}


geologger.read.data<-function( file) {
	track <- read.csv(file,header=FALSE, stringsAsFactors=F) 
	names(track)<-c('date','time','light')
	track$datetime <- paste(track$date,track$time,sep=' ') #makes a new column called datetime with date and time concatenated together with a space between
	track$gmt <- as.POSIXct(strptime(track$datetime,'%d/%m/%Y %H:%M:%S'),'UTC') #makes a new column called gmt with date and time data in a formate that R can work with
	d <- data.frame(id =1:nrow(track), gmt = track$gmt, light = track$light) 
	Data<-list(d=d)
	return(Data)
}

