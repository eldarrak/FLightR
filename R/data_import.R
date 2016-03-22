#######################
# here are the functions that read data


get.tags.data<-function(filename=NULL, start.date=NULL, end.date=NULL, log.light.borders='auto', log.irrad.borders='auto', saves=c("max", "mean"), measurement.period=NULL,  impute.on.boundaries=FALSE) {
 # measurement.period should be set to saving period for swiss tags
   TAGS.twilights<-read.csv(filename, stringsAsFactors =F)

   # now we have to find tag type and figure out whether data were logtransformed..

   detected<-get.tag.type(TAGS.twilights)

   if (is.null(detected)) {
      if (log.light.borders=='auto') stop("Unrecognized tag type, please supply log.light.borders\n")
      if (log.irrad.borders=='auto') stop("Unrecognized tag type, please supply log.irrad.borders\n")
      if (sum(saves==c("max", "mean"))==2) stop("Unrecognized tag type, please supply either tell me what was tag saving - max or mean over saving period\n")
   } else {
      if (detected$log_tranformed) TAGS.twilights$light<-exp(TAGS.twilights$light)
      if(detected$tagtype=="Intigeo_Mode_1") {
	     if (log.light.borders=='auto') log.light.borders<-c(1.5, 9)
	     if (log.irrad.borders=='auto') log.irrad.borders<-c(-3,3)
	     if (sum(saves==c("max", "mean"))==2) saves<-"max"
	  }
      if(detected$tagtype=="Intigeo_Mode_4") {
	     if (log.light.borders=='auto') log.light.borders<-c(1.5, 7)
	     if (log.irrad.borders=='auto') log.irrad.borders<-c(-3,3)
	     if (sum(saves==c("max", "mean"))==2) saves<-"max"
      }
      if(detected$tagtype=="mk") {
	     if (log.light.borders=='auto') log.light.borders<-log(c(2, 63)) 
	     if (log.irrad.borders=='auto') log.irrad.borders<-c(-5.75,1.5)
	     if (sum(saves==c("max", "mean"))==2) saves<-"max"
      }
   }

   FLightR.data<-read.tags.light.twilight(TAGS.twilights,
                      start.date=start.date, end.date=end.date)

   # ok, now we want to specify measurement period
   saving.period<-round(diff(as.numeric(FLightR.data$Data$gmt[1:2]))/10)*10
   if (is.null(measurement.period)) {
      if(saves=="mean") {
	     measurement.period<-saving.period	
      } else {
         measurement.period<-60
      }
   }
   cat("tag saved data every", saving.period, "seconds, and is assumed to measure data every", measurement.period, "seconds, and write down", saves[1], "\n" )
   if (max(TAGS.twilights$light) ==64 & saving.period>500) {
      impute.on.boundaries=TRUE
      cat("saving period was too long for this type of tag, FLightR will impute data\n")
   }
   Proc.data<-process.twilights(FLightR.data$Data,FLightR.data$twilights, 
                             measurement.period=measurement.period, saving.period=saving.period, impute.on.boundaries=impute.on.boundaries)
   if (!is.null(detected)) {
   Proc.data$tagtype<-detected$tagtype
   Proc.data$log_tranformed<-detected$log_tranformed
   Proc.data$log.light.borders=log.light.borders
   Proc.data$log.irrad.borders=log.irrad.borders
   }
   Proc.data$FLightR.data<-FLightR.data
   return(Proc.data)			
}


get.tag.type<-function(TAGS.twilights) {

   Max_light<-max(TAGS.twilights$light)
   recognized<-FALSE
   if(round(Max_light,2)==11.22) {
      tagtype<-"Intigeo_Mode_1"
      log_tranformed<-TRUE
	  recognized<-TRUE
   }
   
   if (round(Max_light/10)==7442) {
     tagtype<-"Intigeo_Mode_1"
	 log_tranformed<-FALSE
     recognized<-TRUE
   }

   if(round(Max_light,2)==7.06) {
      tagtype<-"Intigeo_Mode_4"
      log_tranformed<-TRUE
	  recognized<-TRUE
   }
   if (round(Max_light/10)==116) {
     tagtype<-"Intigeo_Mode_4"
	 log_tranformed<-FALSE
     recognized<-TRUE
   }
   
   if (Max_light == 64) {
      tagtype<-"mk"
	  log_tranformed<-FALSE
	  recognized<-TRUE
   }

   if (Max_light == log(64)) {
      tagtype<-"mk"
	  log_tranformed<-TRUE
	  recognized<-TRUE
   }
   if (recognized==FALSE) { 
   cat("tag type was not recognised!\nmail me details of your tag and I will add them to the list!\n")
   return(NULL)
   } else {
   cat("Detected", tagtype, "tag\n")
   if (log_tranformed) cat("Data found to be logtransformed\n")
   Res<-list(tagtype=tagtype, log_tranformed=log_tranformed)
   return(Res)
  }
}



convert.lux.to.tags<-function(file, log=F, log.light.borders=c(1,10)) {
	# the function takes the current (2015)
	# .lux format and converts it to .csv format that
	# TAGS service accepts
	warning("\n\nwe have figured out that TAGS service is rounding data, so log=T will produce wrong rounding on upload\n\nsetting log to FALSE\n")
	log=F
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
	
	if (!is.null(start.date)) lig.raw<-lig.raw[as.numeric(lig.raw$datetime) > as.numeric(as.POSIXct(start.date , tz="UTC")),]
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
#Filtered_tw$id<-0
#Data$d$type<-0
#Data$d<-rbind(Data$d, data.frame(id=Filtered_tw$id, gmt= Filtered_tw$datetime, light=Filtered_tw$light, type=Filtered_tw$type))

#All.p<-Data$d[order(Data$d$gmt),]
#rownames(All.p)<-1:nrow(All.p)

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

