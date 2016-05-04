

map.FLightR.ggmap<-function(Result, dates=NULL, plot.cloud=TRUE, map.options=NULL, plot.options=NULL, save.options=NULL, zoom="auto", return.ggobj=FALSE) {
if (!is.null(plot.options)) warning("plot options are not in use yet. Let me know what you would like to have here.")
# dates should be a data.frame with first point - starting dates and last column end dates for periods

# ggsave.options is a list that will be will be directly passed to ggsave

library(fields)
library(ggmap)
	 # select twilights to plot
	 if (is.null(dates)) {
	     twilights.index<-1:length(Result$Results$Points.rle)
	 } else {
	      twilights.index<-c()
	        for (segment in 1:nrow(dates)) {
			twilights.index<-c(twilights.index, which(Result$Results$Quantiles$time>=dates[segment,1] & Result$Results$Quantiles$time<=dates[segment,2]))
	        }
			if (length(twilights.index)==0) stop("dates do not overlap with the track time span!")
	 }
	 
    if (plot.cloud) {
         Points_rle<-Result$Results$Points.rle[twilights.index]
         All.Points<-rep(0, nrow(Result$Spatial$Grid))
         for (twilight in 1:length(twilights.index)) {
             All.Points[Points_rle[[twilight]]$values]<-All.Points[Points_rle[[twilight]]$values] + Points_rle[[twilight]]$lengths
         }
		 Points<-Result$Spatial$Grid[All.Points>0,][sample.int(length(All.Points[All.Points>0]), size = 10000, replace = TRUE, prob = All.Points[All.Points>0]), 1:2]
	}
	
	# background map
	
	# check whether Grid was over dateline:
	overdateline<-ifelse(attr(Result$Spatial$Grid, 'left')>	attr(Result$Spatial$Grid, 'right'), TRUE, FALSE)
	
	if (overdateline) {
	location<-cbind(
	            max(Result$Results$Quantiles$Medianlon[twilights.index]),
                min(Result$Results$Quantiles$Medianlat[twilights.index]),
				min(Result$Results$Quantiles$Medianlon[twilights.index]),
				max(Result$Results$Quantiles$Medianlat[twilights.index]))
	} else {			
	location<-cbind(
	            min(Result$Results$Quantiles$Medianlon[twilights.index]),
                min(Result$Results$Quantiles$Medianlat[twilights.index]),
				max(Result$Results$Quantiles$Medianlon[twilights.index]),
				max(Result$Results$Quantiles$Medianlat[twilights.index]))
	}		
	# there could be problems here if the dateline is passed. I will add check for it later.
	center_lon<-mean(c(location[1], location[3]))
	center_lat<-mean(c(location[2], location[4]))
	
	if (center_lon>180 ) center_lon<-center_lon-360
	if ((location[1])>180 ) location[1]<-location[1]-360
	if ((location[3])>180 ) location[3]<-location[3]-360
	
	# combine location with map options
	if (is.null(map.options)) map.options<-list()
	if (is.null(map.options$location)) map.options$location<-c(center_lon, center_lat)
	if (is.null(map.options$zoom) & is.numeric(zoom)) map.options$zoom=zoom
	if (is.null(map.options$col)) map.options$col="bw"
	
    if (zoom=="auto") {

	    for (zoom_cur in (2:10)) {
		   map.options$zoom=zoom_cur
		   background <-do.call(ggmap::get_map, map.options)
		   bb<-attr(background, 'bb')
	      if (overdateline) {
		   bb[2]<-ifelse(bb[2]< (-180), bb[2]+360, bb[2])
		   bb[4]<-ifelse(bb[4]< (-180), bb[4]+360, bb[4])

		   if (!(
		     location[1]<bb[2] &
			 location[2]>bb[1] &
			 location[3]>bb[4] &
			 location[4]<bb[3] )) {
	       break
	       }
		   } else {
		    if (!(
		     location[1]>bb[2] &
			 location[2]>bb[1] &
			 location[3]<bb[4] &
			 location[4]<bb[3] )) {
	       break
	       }
		   }
        }
	    map.options$zoom=zoom_cur-1
    }
	
	
	
    background <-do.call(ggmap::get_map, map.options)
	
	p<-ggmap::ggmap(background, maprange=T)

	if (plot.cloud) {
	 #if (overdateline) Points[,1]<-ifelse(Points[,1]>0, Points[,1]-360, Points[,1])
	 p<-p+stat_density2d(data=data.frame(Points), aes(fill = ..level.., alpha = ..level.., x=lon, y=lat), size = 0.01,  geom = 'polygon', n=400) 
	 
	 p<-p+ scale_fill_gradient(low = "green", high = "red") 
	 p<-p+ scale_alpha(range = c(0.00, 0.25), guide = FALSE) 
     }
	 p<-p+coord_map(projection="mercator", 
     xlim=c(attr(background, "bb")$ll.lon, attr(background, "bb")$ur.lon),
     ylim=c(attr(background, "bb")$ll.lat, attr(background, "bb")$ur.lat)) +
     theme(legend.position = "none", axis.title = element_blank(), text = element_text(size = 12))
    
    # here I plot track for selected dates
	
	Points2plot<-data.frame(lat=Result$Results$Quantiles$Medianlat,lon=Result$Results$Quantiles$Medianlon)
	if (overdateline | max(Points2plot[,2]>attr(background, "bb")$ur.lon)) Points2plot[,2]=Points2plot[,2]-360
	
    if (is.null(dates)) {
        p<-p+ geom_path(data=Points2plot,aes(x=lon,y=lat),  colour=grey(0.3))
		p<-p+  geom_point(data=Points2plot, shape="+",  colour=grey(0.3))
	} else {
	   for (segment in 1:nrow(dates)) {
	      cur_twilights<-which(Result$Results$Quantiles$time>=dates[segment,1] & Result$Results$Quantiles$time<=dates[segment,2])
          p<-p+ geom_path(data=Points2plot[cur_twilights,],aes(x=lon,y=lat),  colour=grey(0.3))+
          geom_point(data=Points2plot[cur_twilights,], shape="+",  colour=grey(0.3))
	   }
	}
	plot(p)
	
    if (is.null(save.options)) save.options<-list()
	if (is.null(save.options$filename)) save.options$filename<-"FLightR.map.pdf"
	do.call(ggsave, save.options)
	if (return.ggobj) return(p)
	}
	

plot.lon.lat<-function(Result, scheme=c("vertical", "horizontal")) {
   Quantiles<-Result$Results$Quantiles
   if (scheme[1]=="horizontal") {
      par(mfrow=c(1,2)) 
   } else {
      par(mfrow=c(2,1)) 
   }
 
   par(mar=c(2,4,3,1),cex=1)
   Sys.setlocale("LC_ALL", "English")  

   #------here I have create the equinox dates..
   Years<-unique(format(Quantiles$time, format="%Y"))
   eq<-c(as.POSIXct(paste(Years, "-09-22 00:00:00 GMT", sep="")), as.POSIXct(paste(Years, "-03-20 12:00:00 GMT", sep="")))
   eq<-eq[eq>min(Quantiles$time) & eq<max(Quantiles$time)]
   #-------
   
   #------- vert grid
   
   Vert_grid<-seq(as.POSIXct("2000-01-01"), as.POSIXct("2030-01-01"), by="month")
   Vert_grid<-Vert_grid[Vert_grid>=min(Quantiles$time) & Vert_grid<=max(Quantiles$time)]
  
   #Longitude
   plot(Quantiles$Medianlon~Quantiles$time, las=1,col=grey(0.1),pch=16,
      ylab="Longitude",xlab="",lwd=2, ylim=range(c( Quantiles$LCI.lon,
      Quantiles$UCI.lon )), type="n", axes=F)
   axis(2)
   axis.POSIXct(1, x=Quantiles$time,  format="1-%b")
   box()

   #################
   # add vertical lines for the first day of every month
   
   abline(v=Vert_grid, col=grey(0.5), lty=2)
   abline(h=seq(-180, 180, by=10), col=grey(0.5), lty=2)
   
   abline(v=eq, col=2, lwd=2, lty=1)
  

   polygon(x=c(Quantiles$time, rev(Quantiles$time)), 
      y=c(Quantiles$LCI.lon, rev(Quantiles$UCI.lon)),
      col=grey(0.9), border=grey(0.5))

   polygon(x=c(Quantiles$time, rev(Quantiles$time)),
   y=c(Quantiles$TrdQu.lon, rev(Quantiles$FstQu.lon)),
   col=grey(0.7), border=grey(0.5))

   lines(Quantiles$Medianlon~Quantiles$time, col=grey(0.1),lwd=2)
   

   #Latitude
   par(mar=c(3,4,1,1))
   
   plot(Quantiles$Medianlat~Quantiles$time, las=1,col=grey(0.1),
     pch=16,ylab="Latitude",xlab="",lwd=2,
     ylim=range(c( Quantiles$UCI.lat, Quantiles$LCI.lat )), type="n")

   abline(v=Vert_grid, col=grey(0.5), lty=2)
   abline(h=seq(-80, 80, by=10), col=grey(0.5), lty=2)

   abline(v=eq, col=2, lwd=2, lty=1)

   polygon(x=c(Quantiles$time, rev(Quantiles$time)), y=c(Quantiles$LCI.lat, rev(Quantiles$UCI.lat)),
           col=grey(0.9), border=grey(0.5))

   polygon(x=c(Quantiles$time, rev(Quantiles$time)), y=c(Quantiles$TrdQu.lat, rev(Quantiles$FstQu.lat)),
           col=grey(0.7), border=grey(0.5))

   lines(Quantiles$Medianlat~Quantiles$time, col=grey(0.1),lwd=2)

}

