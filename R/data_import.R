#' 
#' read TAGS formatted data
#' 
#' Reads the data frame with detected twilight events into the FLightR
#' 
#' @param filename the name of the file which the data are to be read from. File is supposed to be comma separated file of TAGS format. If it does not contain an absolute path, the file name is relative to the current working directory, getwd(). Tilde-expansion is performed where supported. This can be a compressed file (see \code{\link[base]{file}}). Alternatively, file can be a readable text-mode connection (which will be opened for reading if necessary, and if so closed (and hence destroyed) at the end of the function call). File can also be a complete URL. For the supported URL schemes, see help for \code{\link[base]{url}}.
#' @param start.date date of beginning of relevant data collection in \code{\link[base]{POSIXct}} format.
#' @param end.date date of end of relevant data collection in \code{\link[base]{POSIXct}} format.
#' @param log.light.borders Numeric vector with length of 2 for minimum and maximum log(light) levels to use. Alternatively character value 'auto', that will allow FLightR to assign these values according to detected tag type.
#' @param log.irrad.borders Numeric vector with length of 2 for minimum and maximum log(irradiance) values to use. Alternatively character value 'auto', that will allow FLightR to assign these values according to detected tag type.
#' @param saves character values informing FLightR if min or max values were used by logger.
#' @param measurement.period Value in seconds defining how often tag was measuring light levels. If NULL value will be taken from known values for detected tag type.
#' @param impute.on.boundaries logical, if FLightR should approximate values at boundaries. Set it to TRUE only if you have vary few active points at each twilight, e.g if tag was saving every 10 minutes or so.
#' @return list, which is to be further processed with the FLightR.
#' @details The returned object has many parts, the important are: (1) the recorded light data, (2) the detected twilight events, (3) light level data at the moment of each determined sunrise and sunset and around them (24 fixes before and 24 after), and (4) technical parameters of the tag, i. e. its type, saving and measuring period (the periodicity, in seconds, at which a tag measures and saves data).
#' @examples
#' File<-system.file("extdata", "Godwit_TAGS_format.csv", package = "FLightR")
#' Proc.data<-get.tags.data(File)
#' @export
get.tags.data<-function(filename=NULL, start.date=NULL, end.date=NULL, log.light.borders='auto', log.irrad.borders='auto', saves=c("auto", "max", "mean"), measurement.period=NULL,  impute.on.boundaries=FALSE) {
 # measurement.period should be set to saving period for swiss tags
   TAGS.twilights<-utils::read.csv(filename, stringsAsFactors =FALSE)

   # now we have to find tag type and figure out whether data were logtransformed..

   detected<-get.tag.type(TAGS.twilights)

   if (is.null(detected)) {
      if (log.light.borders[1]=='auto') stop("Unrecognized tag type, please supply log.light.borders\n")
      if (log.irrad.borders[1]=='auto') stop("Unrecognized tag type, please supply log.irrad.borders\n")
      if (length(saves)==3 | saves[1]=='auto') stop("Unrecognized tag type, please  tell FLightR what was tag saving - max or mean over saving period\n")
   } else {
      if (detected$log_transformed) TAGS.twilights$light<-exp(TAGS.twilights$light)
      if(detected$tagtype=="Intigeo_Mode_1") {
	     if (log.light.borders[1]=='auto') log.light.borders<-c(1.5, 9)
	     if (log.irrad.borders[1]=='auto') log.irrad.borders<-c(-3,3)
	     if (saves[1] =='auto') saves<-"max"
	  }
      if(detected$tagtype=="Intigeo_Mode_4") {
	     if (log.light.borders[1]=='auto') log.light.borders<-c(1.5, 7)
	     if (log.irrad.borders[1]=='auto') log.irrad.borders<-c(-3,3)
	     if (saves[1] =='auto') saves<-"max"
      }
      
      if(detected$tagtype=="Intigeo_Mode_6_clipped") {
	     if (log.light.borders[1]=='auto') log.light.borders<-c(1.5, 7)
	     if (log.irrad.borders[1]=='auto') log.irrad.borders<-c(-3,3)
	     if (saves[1] =='auto') saves<-"max"
      }
      if(detected$tagtype=="mk") {
	     if (log.light.borders[1]=='auto') log.light.borders<-log(c(4, 61)) 
	     if (log.irrad.borders[1]=='auto') log.irrad.borders<-c(-6.5,4)
	     if (saves[1] =='auto') saves<-"max"
      }
	  if(detected$tagtype=="GDL2v2_GDLpam") {
	     if (log.light.borders[1]=='auto') log.light.borders<-c(2.5, 8)
	     if (log.irrad.borders[1]=='auto') log.irrad.borders<-c(-6.5,1.5)
	     if (saves[1] =='auto') saves<-"mean"
      }
	  if(detected$tagtype=="GDL1") {
	     if (log.light.borders[1]=='auto') log.light.borders<-c(3, 7)
	     if (log.irrad.borders[1]=='auto') log.irrad.borders<-c(-4,1)
	     if (saves[1] =='auto') saves<-"mean"
      }
	  if(detected$tagtype=="GDL2v1") {
	     if (log.light.borders[1]=='auto') log.light.borders<-c(2.2,7)
	     if (log.irrad.borders[1]=='auto') log.irrad.borders<-c(-2.7,0)
	     if (saves[1] =='auto') saves<-"mean"
      }
    if(detected$tagtype=="Lat_2000") {
	     if (log.light.borders[1]=='auto') log.light.borders<-c(100,360)
	     if (log.irrad.borders[1]=='auto') log.irrad.borders<-c(-8,2)
	     if (saves[1] =='auto') saves<-"max"
         measurement.period<-10
    }
      if(detected$tagtype=="mutag") {
        if (log.light.borders[1]=='auto') log.light.borders<-c()
        if (log.irrad.borders[1]=='auto') log.irrad.borders<-c()
        if (saves[1] =='auto') saves<-"max"
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
   
   message("tag saved data every", saving.period, "seconds, and is assumed to measure data every", measurement.period, "seconds, and write down", saves[1], "\n" )
   if (max(TAGS.twilights$light) ==64 & saving.period>500) {
      impute.on.boundaries=TRUE
      warning("saving period was too long for this type of tag, FLightR will impute data\n")
   }
   Proc.data<-process.twilights(FLightR.data$Data,FLightR.data$twilights, 
                             measurement.period=measurement.period, saving.period=saving.period, impute.on.boundaries=impute.on.boundaries)
   if (!is.null(detected)) {
   Proc.data$tagtype<-detected$tagtype
   Proc.data$log_transformed<-detected$log_transformed
   }
   Proc.data$log.light.borders=log.light.borders
   Proc.data$log.irrad.borders=log.irrad.borders
   Proc.data$FLightR.data<-FLightR.data
   return(Proc.data)			
}


get.tag.type<-function(TAGS.twilights) {

   Max_light<-max(TAGS.twilights$light)
   Min_light<-min(TAGS.twilights$light)
   recognized<-FALSE
   
   if(Max_light == 1146.681 &  Min_light == 0.32 ) {
      tagtype<-"Intigeo_Mode_6_clipped"
      log_transformed<-FALSE
	  recognized<-TRUE
   }
   
      
   if(Max_light == log(1146.681) &  Min_light == log(0.32) ) {
      tagtype<-"Intigeo_Mode_6_clipped"
      log_transformed<-TRUE
	  recognized<-TRUE
   }
   
   if(round(Max_light,2) >= 10.05 &  (round(Max_light,2) <= 11.22 | round(Max_light,2) >= 11.22 & round(Max_light,2) <=11.35)) {
      tagtype<-"Intigeo_Mode_1"
      log_transformed<-TRUE
	  recognized<-TRUE
   }
   
   if(round(Max_light,2) >= 10.05 &  (round(Max_light,2) <= 11.22 | round(Max_light,2) >= 11.22 & round(Max_light,2) <=11.35)) {
      tagtype<-"Intigeo_Mode_1"
      log_transformed<-TRUE
	  recognized<-TRUE
   }
   
  #if (round(Max_light/10) >= 2000 & round(Max_light/10)<= 7000) {
   if(round(Max_light,2) >= exp(10.05) &  (round(Max_light,2) <= exp(11.22) | round(Max_light,2) >= exp(11.22) & round(Max_light,2) <= exp(11.35))) {
     tagtype<-"Intigeo_Mode_1"
	 log_transformed<-FALSE
     recognized<-TRUE
   }

   if(round(Max_light,2)==7.06) {
      tagtype<-"Intigeo_Mode_4"
      log_transformed<-TRUE
	  recognized<-TRUE
   }
   if (round(Max_light/10)==116) {
     tagtype<-"Intigeo_Mode_4"
	 log_transformed<-FALSE
     recognized<-TRUE
   }
   
   if (Max_light == 64) {
      tagtype<-"mk"
	  log_transformed<-FALSE
	  recognized<-TRUE
   }

   if (Max_light == log(64)) {
      tagtype<-"mk"
	  log_transformed<-TRUE
	  recognized<-TRUE
   }
   if(round(Max_light,1)==9.2) {
      tagtype<-"GDL2v2_GDLpam"
      log_transformed<-TRUE
	  recognized<-TRUE
   }
   if (round(Max_light/10)==998) {
     tagtype<-"GDL2v2_GDLpam"
	 log_transformed<-FALSE
     recognized<-TRUE
   }
   if(round(Max_light,2)==11.32) {
      tagtype<-"GDL1"
      log_transformed<-TRUE
	  recognized<-TRUE
   }
   if (round(Max_light)==82863) {
     tagtype<-"GDL1"
	 log_transformed<-FALSE
     recognized<-TRUE
   }
   if(round(Max_light,2)==log(63)) {
      tagtype<-"GDL2v1"
      log_transformed<-TRUE
	  recognized<-TRUE
   }
   if (Max_light==63) {
     tagtype<-"GDL2v1"
	 log_transformed<-FALSE
     recognized<-TRUE
   }
   if(log(Max_light) %in% c(4095, 357)) {
      tagtype<-"Lat_2000"
      log_transformed<-FALSE
	  recognized<-TRUE
   }

   if(Max_light %in% c(4095, 357)) {
      tagtype<-"Lat_2000"
      log_transformed<-TRUE
	  recognized<-TRUE
   }

   if (Max_light == 10240) {
     tagtype <- "mutag"
     log_transformed <- FALSE
     recognized <- TRUE
   }

   if (recognized==FALSE) {
    warning("tag type was not recognised!\nmail me details of your tag and I will add them to the list!\n")
    return(NULL)
   } else {
   message("Detected", tagtype, "tag\n")
   if (log_transformed) message("Data found to be logtransformed\n")
   Res<-list(tagtype=tagtype, log_transformed=log_transformed)
   return(Res)
  }
}


convert.lux.to.tags<-function(file, log=FALSE, log.light.borders=c(1,10)) {
	# the function takes the current (2015)
	# .lux format and converts it to .csv format that
	# TAGS service accepts
	warning("\n\nwe have figured out that TAGS service is rounding data, so log=TRUE will produce wrong rounding on upload\n\nsetting log to FALSE\n")
	log=FALSE
	Dat<-utils::read.csv(file, skip=20, sep="\t", stringsAsFactors =FALSE)
	names(Dat)<-c("datetime", "light")
	Dat$datetime<-as.POSIXct(Dat$datetime, tz="GMT", format="%d/%m/%Y %H:%M:%S")

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
	utils::write.csv(Dat_new, file=paste(unlist(strsplit(file, ".lux")), "csv", sep="."), quote =FALSE, row.names=FALSE) 
	message("Success!\n")
	message("file", paste(unlist(strsplit(file, ".lux")), "csv", sep="."), "\nwas saved to", getwd(), "\n")
	return(NULL)
}


read.tags.light.twilight<-function(lig.raw, start.date=NULL, end.date=NULL) {
	lig.raw$datetime<-as.POSIXct(lig.raw$datetime, tz="GMT", format="%Y-%m-%dT%T")

	###
	### I also want to exclude the last days...
	
	if (!is.null(start.date)) lig.raw<-lig.raw[as.numeric(lig.raw$datetime) > as.numeric(as.POSIXct(start.date , tz="GMT")),]
	if (!is.null(end.date)) lig.raw<-lig.raw[as.numeric(lig.raw$datetime) < as.numeric(as.POSIXct(end.date , tz="GMT")),]

	## and also I want to exclude the interpolated and excluded afterwards points
		
	lig<-lig.raw[-which(lig.raw$interp==TRUE),]

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
	utils::write.table(lig.new, file=tmpname, sep="," , row.names = FALSE, col.names =FALSE,  quote=FALSE)
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

Filtered_tw$id<-0
Data$d$type<-0
Data$d<-rbind(Data$d, data.frame(id=Filtered_tw$id, gmt= Filtered_tw$datetime, light=Filtered_tw$light, type=Filtered_tw$type))

Filtered_tw$excluded=0

All.p<-Data$d[order(Data$d$gmt),]

#All.p<-All.p[!duplicated(All.p[,2:3], fromLast=TRUE),]
rownames(All.p)<-1:nrow(All.p)
Res<-list(Data=All.p, twilights=Filtered_tw)
return(Res)
}


process.geolight.output<-function(datetime, light, gl_twl) {

#filter <- loessFilter(gl_twl[,1],gl_twl[,2],gl_twl[,3],k=10)
Filtered_tw <- data.frame(datetime=as.POSIXct(c(gl_twl$tFirst,gl_twl$tSecond),"GMT"),type=c(gl_twl$type,ifelse(gl_twl$type==1,2,1)))

Filtered_tw <- Filtered_tw[!duplicated(Filtered_tw$datetime),]
Filtered_tw <- Filtered_tw[order(Filtered_tw[,1]),]

# now I want to pair data and twilights..		  
Filtered_tw$light<-stats::approx(x=datetime, y=light, xout=Filtered_tw$datetime)$y
Filtered_tw$id<-0
Data<-data.frame(id=1:length(datetime), gmt=datetime, light=light, type=0)

Data<-rbind(Data, data.frame(id=Filtered_tw$id, gmt= Filtered_tw$datetime, light=Filtered_tw$light, type=Filtered_tw$type))

All.p<-Data[order(Data$gmt),]
#All.p<-All.p[!duplicated(All.p[,2:3], fromLast=TRUE),]
rownames(All.p)<-1:nrow(All.p)
Filtered_tw$excluded=0
FLightR.data=list(Data=All.p, twilights=Filtered_tw)
return(FLightR.data)

}


geologger.read.data<-function( file) {
	track <- utils::read.csv(file,header=FALSE, stringsAsFactors=FALSE) 
	names(track)<-c('date','time','light')
	track$datetime <- paste(track$date,track$time,sep=' ') #makes a new column called datetime with date and time concatenated together with a space between
	track$gmt <- as.POSIXct(strptime(track$datetime,'%d/%m/%Y %H:%M:%S'),tz='GMT') #makes a new column called gmt with date and time data in a formate that R can work with
	d <- data.frame(id =1:nrow(track), gmt = track$gmt, light = track$light) 
	Data<-list(d=d)
	return(Data)
}
