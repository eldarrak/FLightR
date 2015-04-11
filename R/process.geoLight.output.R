#######
# here I'll write a function that will take GeoLight output..

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
