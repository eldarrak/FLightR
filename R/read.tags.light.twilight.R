read.tags.light.twilight<-function(lig.raw, end.date=NULL) {
	lig.raw$datetime<-as.POSIXct(lig.raw$datetime, tz="UTC", format="%Y-%m-%dT%T")

	###
	### I also want to exclude the last days...
	
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