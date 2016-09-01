#' plots result over map
#'
#' plots track over map with probability cloud. Can plot only part of the track if dates are specified
#' @param Result FLightR result object obtained from \code{\link{run.particle.filter}}
#' @param dates either NULL if all twilights should be included or data.frame with first colum - start of the period and second end of the period. Each line represents a new period
#' @param plot.cloud Shlould probability cloud be plotted? If TRUE cloud is estimated by \code{\link[ggplot2]{stat_density2d}}
#' @param map.options options passed to \code{\link[ggmap]{get_map}}, note that \code{zoom} option is defined separately
#' @param plot.options plotting options. Not defined yet!
#' @param save.options ptions passed to \code{\link[ggplot2]{ggsave}}. Filename should be defined here.
#' @param zoom Zoom for map. If 'auto' FLightR will try to find optimal zoom level by downloading different size maps and checking whether all the points fit the map.
#' @param return.ggobj Should ggobj be returned for subsequent checks and/or replotting
#' @param seasonal.colors if true points of the track will have seasonal colors
#' @param seasonal.donut.location if NULL - no color wheel placed, otherwise select one of 'bottomleft', 'bottomright', 'topleft'
#' @param seasonal.donut.proportion how much of X axis should color wheel occupy.
#' return either NULL or ggplot2 class object
#' @export map.FLightR.ggmap
map.FLightR.ggmap<-function(Result, dates=NULL, plot.cloud=TRUE, map.options=NULL, plot.options=NULL, save.options=NULL, zoom="auto", return.ggobj=FALSE, seasonal.colors=TRUE, seasonal.donut.location='topleft', seasonal.donut.proportion=0.5) {
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
	            min(Result$Results$Quantiles$Medianlon[twilights.index][Result$Results$Quantiles$Medianlon[twilights.index]>0]),
				min(Result$Results$Quantiles$Medianlat[twilights.index]),
				max(Result$Results$Quantiles$Medianlon[twilights.index][Result$Results$Quantiles$Medianlon[twilights.index]<0]),
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
	
	if (overdateline) center_lon<-mean(c(location[1], location[3]+360))-360
	
	if (center_lon < (-180)) center_lon<-center_lon+360
	
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
	      
		   bb[2]<-ifelse(bb[2]< (-180), bb[2]+360, bb[2])
		   bb[4]<-ifelse(bb[4]< (-180), bb[4]+360, bb[4])
		   bb[2]<-ifelse(bb[2]> (180), bb[2]-360, bb[2])
		   bb[4]<-ifelse(bb[4]> (180), bb[4]-360, bb[4])
		  
		if (bb[4]<bb[2]) {
           lonisinbb<-(location[1] >= bb[2] || location[1] <= bb[4]) & (location[3] >= bb[2] || location[3] <= bb[4]) 
        } else {
	       lonisinbb<-(location[1] >= bb[2] && location[1] <= bb[4]) & (location[3] >= bb[2] && location[3] <= bb[4])
	    }
	  
	    isinbb<-lonisinbb & location[2]>=bb[1] & location[4]<=bb[3]
	    if (!isinbb) break
	

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
	if (overdateline) Points2plot[,2][Points2plot[,2]>attr(background, "bb")$ur.lon] <- Points2plot[,2][Points2plot[,2]>attr(background, "bb")$ur.lon]-360
	if (overdateline)  Points2plot[,2][Points2plot[,2]<attr(background, "bb")$ll.lon] <- Points2plot[,2][Points2plot[,2]<attr(background, "bb")$ll.lon] + 360
	
    if (is.null(dates)) {
	    cur_twilights<-1:length(Result$Results$Quantiles$time)
		Total_segments<-1
	} else {
	   cur_twilights<-which(Result$Results$Quantiles$time>=dates[segment,1] & Result$Results$Quantiles$time<=dates[segment,2])
	   Total_segments<-nrow(dates)
	
	}
	for (segment in 1:Total_segments) {
       p<-p+ geom_path(data=Points2plot[cur_twilights,],aes(x=lon,y=lat),  colour=grey(0.3))
	   if (seasonal.colors) {
		  Seasonal_palette<-colorRampPalette(hsv(1-((1:365)+(365/4))%%365/365, s=0.8, v=0.8), space="Lab")
          Days<-as.integer(format(Result$Results$Quantiles$time[cur_twilights], "%j"))
		  p<-p+geom_point(data=Points2plot[cur_twilights,], shape=21,  colour=grey(0.3), bg=Seasonal_palette(365)[Days])
	   } else {
		  p<-p+geom_point(data=Points2plot[cur_twilights,], shape="+",  colour=grey(0.3))
	   }
	}
	
	if (!is.null(seasonal.donut.location)) {
       Xrange<-c(attr(background, "bb")$ll.lon, attr(background, "bb")$ur.lon)
       Yrange<-c(attr(background, "bb")$ll.lat, attr(background, "bb")$ur.lat)
       d<-seasonal_donut()  
       g = ggplotGrob(d)

    if (seasonal.donut.location=='bottomleft') {
        p<-p + inset(grob = g, xmin = Xrange[1], xmax = Xrange[1]+seasonal.donut.proportion*(Yrange[2]-Yrange[1]), ymin = Yrange[1], ymax = Yrange[1]+seasonal.donut.proportion*(Yrange[2]-Yrange[1]))
    }

    if (seasonal.donut.location=='topleft') {
        p<-p + inset(grob = g, xmin = Xrange[1], xmax = Xrange[1]+seasonal.donut.proportion*(Yrange[2]-Yrange[1]), ymin = Yrange[2]-seasonal.donut.proportion*(Yrange[2]-Yrange[1]), ymax = Yrange[2])
    }
	if (seasonal.donut.location=='bottomright') {
        p<-p + inset(grob = g, xmin = Xrange[2]-seasonal.donut.proportion*(Yrange[2]-Yrange[1]), xmax = Xrange[2], ymin = Yrange[1], ymax = Yrange[1]+seasonal.donut.proportion*(Yrange[2]-Yrange[1]))
    }

    }
	plot(p)
	
    if (is.null(save.options)) save.options<-list()
	if (is.null(save.options$filename)) save.options$filename<-"FLightR.map.pdf"
	do.call(ggsave, save.options)
	if (return.ggobj) return(p)
	}

#' plots result by longitude and latitude
#'
#' This function plots result by latitude and longitude in either vertical or horizontal layout
#' @param Result FLightR result object obtained from \code{\link{run.particle.filter}}
#' @param scheme either 'vertical' or 'horizontal' layouts
#' return NULL
#' @export plot.lon.lat
plot.lon.lat<-function(Result, scheme=c("vertical", "horizontal")) {

   Quantiles<-Result$Results$Quantiles
   
   
   	# check whether Grid was over dateline:
	overdateline<-ifelse(attr(Result$Spatial$Grid, 'left')>	attr(Result$Spatial$Grid, 'right'), TRUE, FALSE)

   
   if (overdateline) {
   
   Quantiles$Medianlon[Quantiles$Medianlon<0]<-Quantiles$Medianlon[Quantiles$Medianlon<0]+360
   Quantiles$LCI.lon[Quantiles$LCI.lon<0]<-Quantiles$LCI.lon[Quantiles$LCI.lon<0]+360
   Quantiles$UCI.lon[Quantiles$UCI.lon<0]<-Quantiles$UCI.lon[Quantiles$UCI.lon<0]+360
   Quantiles$TrdQu.lon[Quantiles$TrdQu.lon<0]<-Quantiles$TrdQu.lon[Quantiles$TrdQu.lon<0]+360
   Quantiles$FstQu.lon[Quantiles$FstQu.lon<0]<-Quantiles$FstQu.lon[Quantiles$FstQu.lon<0]+360
   
    }
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
   axis(2, las=1)
   axis.POSIXct(1, x=Quantiles$time,  format="1-%b")
   box()

   #################
   # add vertical lines for the first day of every month
   
   abline(v=Vert_grid, col=grey(0.5), lty=2)
   abline(h=seq(-180, 360, by=10), col=grey(0.5), lty=2)
   
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
     ylim=range(c( Quantiles$UCI.lat, Quantiles$LCI.lat )), type="n", axes=F)
   axis(2, las=1)
   axis.POSIXct(1, x=Quantiles$time,  format="1-%b")
   box()
   abline(v=Vert_grid, col=grey(0.5), lty=2)
   abline(h=seq(-80, 80, by=10), col=grey(0.5), lty=2)

   abline(v=eq, col=2, lwd=2, lty=1)

   polygon(x=c(Quantiles$time, rev(Quantiles$time)), y=c(Quantiles$LCI.lat, rev(Quantiles$UCI.lat)),
           col=grey(0.9), border=grey(0.5))

   polygon(x=c(Quantiles$time, rev(Quantiles$time)), y=c(Quantiles$TrdQu.lat, rev(Quantiles$FstQu.lat)),
           col=grey(0.7), border=grey(0.5))

   lines(Quantiles$Medianlat~Quantiles$time, col=grey(0.1),lwd=2)

}


get_buffer<-function(coords, r){
  p = SpatialPoints(matrix(coords, ncol=2), proj4string=CRS("+proj=longlat +datum=WGS84"))
  aeqd <- sprintf("+proj=aeqd +lat_0=%s +lon_0=%s +x_0=0 +y_0=0",
                    p@coords[[2]], p@coords[[1]])
  projected <- spTransform(p, CRS(aeqd))
  buffered <- gBuffer(projected, width=r, byid=TRUE)
  buffer_lonlat <- spTransform(buffered, CRS=p@proj4string)
  return(buffer_lonlat)
  }
  
get_gunion_r<-function(Result) {
    Distances=  spDists(Result$Spatial$Grid[1:which.min(c(nrow(Result$Spatial$Grid), 1000)),1:2], longlat=T)
    # ok, distances go up to 51.2.. the next step is 62.. 
    # so if I round them 
    Selected_dist<-unique(sort(round(Distances/10)*10))[2]
    r<-Selected_dist*0.85*1000
	return(r)
  }
  
get_time_spent_buffer<-function(Result, dates=NULL, percentile=0.5, r=NULL) {
# r in meters.. 
# dates could be either NULL - for plotting all available dates,
# or numeric with length of one - to plot a specific twilight by number (row number in Quantiles)
# or data.frame with two columns - one for the first date in period and one for the last date in period.
	 if (is.null(dates)) {
	     twilights.index<-1:length(Result$Results$Points.rle)
	 } else {
	    if (is.numeric(dates) & length(dates)==1) {
		  twilights.index<-dates
		} else {
	      twilights.index<-c()
	      for (segment in 1:nrow(dates)) {
			twilights.index<-c(twilights.index, which(Result$Results$Quantiles$time>=dates[segment,1] & Result$Results$Quantiles$time<=dates[segment,2]))
	        }
			if (length(twilights.index)==0) stop("dates do not overlap with the track time span!")
		}
	  }
	 cat('function will plot', length(twilights.index), 'twilights\n')
         Points_rle<-Result$Results$Points.rle[twilights.index]
         All.Points<-rep(0, nrow(Result$Spatial$Grid))
         for (twilight in 1:length(twilights.index)) {
             All.Points[Points_rle[[twilight]]$values]<-All.Points[Points_rle[[twilight]]$values] + Points_rle[[twilight]]$lengths
         }
		 
  nParticles<- sum(Result$Results$Points.rle[[1]]$lengths)
  
  Order<-order(All.Points, decreasing=TRUE)
  
  Order_rle<-rle(All.Points[Order]) # 
  
  Cum_probs<-cumsum(Order_rle$values/nParticles/length(twilights.index))
  
  #Percentile=percentile
  tmp<-which(Cum_probs<=percentile)
  if (length(tmp)==0) { 
  Selected_rle<-1 # if there is no such bin - return the most likely one
  #return(ceiling(min(Cum_probs)*100)/100)
  } else {
  Selected_rle<-max(which(Cum_probs<=percentile))
  }
  #Points<-which(All.Points/nParticles/length(twilights.index)>percentile)
  Points<-Order[1:sum(Order_rle$lengths[1:Selected_rle])]
 
  #if (length(Points)==0) { 
       #warning("bird dod not spend ", percentile, " of time in any point:\n     maximum time spent in on grid cell is ", floor(max(All.Points/nParticles/length(twilights.index))*100)/100, "\n")

  
  
  # so this is the area where bird spent almost all the time (95%)
  # how should we show it?
  # yeah it could be some kind of a polygon or 
  # if will combind points by distance
  #  
  
  Points_selected<-(1:length(All.Points))[Points]
 if (is.null(r)) {
     r=get_gunion_r(Result)
     }
  #Now I want to create a spatial buffers around points with this radius
  
  #spoints = SpatialPoints(Result$Spatial$Grid[Points_selected,1:2], proj4string= CRS("+proj=longlat +datum=WGS84"))
  
  Buffers<-apply(matrix(Result$Spatial$Grid[Points_selected,1:2], ncol=2), 1, get_buffer, r=r)
   Buff_comb<-Buffers[[1]] 
   if (length(Buffers)>1) {
   for (i in 2:length(Buffers)) {
       Buff_comb<-gUnion(Buff_comb, Buffers[[i]])
  }
  }
  Buff_comb_simpl<-gSimplify(Buff_comb, tol=0.01, topologyPreserve=T)
  return(list(Buffer=Buff_comb_simpl, nPoints=length(Points)))
  }

#' plots resulting track over map with uncertainty shown by spave utilisation distribution
#' 
#' May be use not only for the whole track but for a set of specific dates, e.g. to show spatial uncertainty during migration
#' 
#' @param Result FLightR result object obtained from \code{\link{run.particle.filter}}
#' @param dates Use NULL if all twilights will be used for plotting, one integer if specific twilight should be plotted (line number in Result$Results$Quantiles). Use data.frame with first colum - start of the period and second - end of the period and each line represents a new period to plot specific periods, e.g. wintering or migration.
#' @param map.options options passed to \code{\link[ggmap]{get_map}}, note that \code{zoom} option is defined separately
#' @param percentiles Probability breaks for utilisation distribution
#' @param zoom Zoom for map. If 'auto' FLightR will try to find optimal zoom level by downloading different size maps and checking whether all the points fit the map.
#' @param geom_polygon.options passed to \code{\link[ggplot2]{geom_polygon}}
#' @param save.options ptions passed to \code{\link[ggplot2]{ggsave}}. Filename should be defined here.
#' @param color.palette colors for probability contours. Either NULL or \code{\link[grDevices]{colorRampPalette}} object
#' @param use.palette should the same colors be used for polygon boundaries as for polygon filling?
#' @param background if provided will be used as a background. Must be created by \code{link[ggmap]{get_map}}
#' @param plot should function produce a plot?
#' @param save should function save results with \code{\link[ggplot2]{ggsave}}?
#' @param add.scale.bar will add scalebar with the \code{\link[ggsn]{scalebar}}
#' @param scalebar.options options passed to \code{\link[ggsn]{scalebar}}
#' @return list with two parts 
#'         \item{res_buffers}{spatial buffers for defined probability values}
#'         \item{p}{\code{\link[ggplot2]{ggplot}} object}
#' @author Eldar Rakhimberdiev
#' @export plot.util.distr
plot.util.distr<-function(Result, dates=NULL, map.options=NULL, percentiles=c(0.4, 0.6, 0.8), zoom="auto", geom_polygon.options=NULL, save.options=NULL, color.palette=NULL, use.palette=TRUE, background=NULL, plot=TRUE, save=TRUE, add.scale.bar=FALSE, scalebar.options=NULL) {

if (is.null(color.palette)) color.palette <- colorRampPalette(rev(c("#edf8fb",
                                    "#b3cde3",
									"#8c96c6",
									"#88419d")))
# percentiles select area where bird spend at least xx percent of it's time
percentiles<-sort(percentiles)

r=get_gunion_r(Result)

#--------------
# Now I want to figure out minimum percentile plossible

res_buffers<-c()
nPoints<-c()
for (percentile in percentiles) {
res_cur<- get_time_spent_buffer(Result, dates, percentile, r)
res_buffers<-c(res_buffers, res_cur$Buffer)
nPoints<-c(nPoints, res_cur$nPoints)
}

# Now I want to take only the upper levels from identical nPoints
Index_not_dupl<-(1:length(percentiles))[!duplicated(nPoints)]

#-----------------#
# and now we want to plot it

if (is.null(background)) {

# combine location with map options
	if (is.null(map.options)) map.options<-list()
	#if (is.null(map.options$location)) map.options$location<-location
	if (is.null(map.options$zoom) & is.numeric(zoom)) map.options$zoom=zoom
	if (is.null(map.options$col)) map.options$col="bw"
	Extent<-extent(res_buffers[[length(res_buffers)]])
	location<-c(Extent@xmin, Extent@ymin, Extent@xmax, Extent@ymax)
    background <-do.call(ggmap::get_map, map.options)
    if (is.null(map.options$location)) map.options$location<-location

	if (zoom=="auto") {

	for (zoom_cur in (2:10)) {
	map.options$zoom=zoom_cur
	
	background <-do.call(ggmap::get_map, map.options)
	if (!(
	location[1]>attr(background, 'bb')[2] &
	location[2]>attr(background, 'bb')[1] &
	location[3]<attr(background, 'bb')[4] &
	location[4]<attr(background, 'bb')[3] )) {
	break
	}
    }
	map.options$zoom=zoom_cur-1
    }
	
	background <-do.call(ggmap::get_map, map.options)
    }
    geom_polygon.options.external=geom_polygon.options
	
	p<-ggmap::ggmap(background)
	#p<-ggmap::ggmap(background, extent = "normal", maprange=FALSE)
	#p<-ggmap::ggmap(background, extent = "normal")
	
	Colors<-color.palette(length(res_buffers))
	
	for (i in length(res_buffers):1) {
	
	buff_cur<-res_buffers[[i]]
	
	if (is.null(geom_polygon.options)) geom_polygon.options=list()
	if (is.null(geom_polygon.options$alpha)) geom_polygon.options$alpha=0.5
	    
    geom_path.options<-list()
    geom_path.options$color<- geom_polygon.options$color
	if (is.null(geom_path.options$color)) geom_path.options$color =ifelse(use.palette, Colors[i], '#756bb1')

	if (is.null(geom_polygon.options$fill)) geom_polygon.options$fill =geom_path.options$color

	geom_path.options$size<- geom_polygon.options$size
	if (is.null(geom_path.options$size)) geom_path.options$size=2
		
	geom_path.options$size<- geom_polygon.options$size
	if (is.null(geom_path.options$size)) geom_path.options$size=2

	geom_path.options$alpha=0.8
	
	geom_path.options$linetype<- geom_polygon.options$linetype
	geom_polygon.options$linetype=0
	
	geom_polygon.options$mapping=aes(x=long, y=lat, group=group)
	geom_path.options$mapping=geom_polygon.options$mapping
	
	geom_polygon.options$data=buff_cur
	geom_path.options$data=geom_polygon.options$data
	p<-p+do.call(ggplot2::geom_polygon, geom_polygon.options)
	p<-p+do.call(ggplot2::geom_path, geom_path.options)
	
	geom_polygon.options=geom_polygon.options.external
	
	}
    if (plot) print(p)
	if (add.scale.bar) {
	if (is.null(scalebar.options)) scalebar.options=list()
	BB<-attr(background, 'bb')

	if (is.null(scalebar.options$dist)) scalebar.options$dist=max(((abs(as.numeric(BB[1, 4])-as.numeric(BB[1, 2]))*100*0.1)%/%25)*25, 25)
	
	scalebar.options$dd2km <- TRUE
	scalebar.options$model <- 'WGS84'
	if (is.null(scalebar.options$location)) scalebar.options$location <- 'topright'
	scalebar.options$y.min=as.numeric(BB[1, 1])+(as.numeric(BB[1, 3])-as.numeric(BB[1, 1]))*0.05
	scalebar.options$y.max=as.numeric(BB[1, 3])-(as.numeric(BB[1, 3])-as.numeric(BB[1, 1]))*0.05
	scalebar.options$x.min=as.numeric(BB[1, 2])-(as.numeric(BB[1, 2])-as.numeric(BB[1, 4]))*0.05
	scalebar.options$x.max=as.numeric(BB[1, 4])+(as.numeric(BB[1, 2])-as.numeric(BB[1, 4]))*0.05
	
	p=p+do.call(ggsn::scalebar, scalebar.options)
	}
	
	if (save) {
	   if (is.null(save.options))  save.options=list()
	   if (is.null(save.options$filename)) save.options$filename<-"time_coverage.png"
	   save.options$plot <-p
	   save.options$dpi <-600
	   tmp<-do.call(ggplot2::ggsave, save.options)
	}
	return(list(res_buffers=res_buffers, p=p))
}


 
seasonal_donut<-function() {
   pie.data<-data.frame(group=as.factor(c(1:12)), value=rep(1,12))  
   pie.data$fraction = pie.data$value / sum(pie.data$value)
   pie.data$ymax = cumsum(pie.data$fraction)
   pie.data$ymin = c(0, head(pie.data$ymax, n = -1))
   
   donut<-ggplot(data = pie.data, aes(fill = group, ymax = ymax, ymin = ymin, xmax = 1, xmin = 2))  +
      geom_rect(colour = "grey30", show.legend = FALSE) +
      coord_polar(theta = "y", start=pi) +
      xlim(c(0, 2)) +
      theme_bw() + 
	  theme(plot.background = element_rect(fill = "transparent" ,colour = NA))+
	  theme(panel.background = element_blank())+
	  theme(panel.grid=element_blank()) +
      theme(axis.text=element_blank()) +
      theme(axis.ticks=element_blank()) +
      theme(panel.border=element_blank()) +
	  scale_fill_manual(values=c(Seasonal_palette(13)))+
      geom_text(aes(x = 1.5, y = ((ymin+ymax)/2), label = group), colour=grey(0.99), size=5) +
      xlab("") +
      ylab("") 
	  return(donut)
}


# function below is at the development stage..
plot.likelihood<-function(object, date=NULL, twilight.number=NULL) {
my.golden.colors <- colorRampPalette(c("white","#FF7100"))

image.plot(as.image(object$Spatial$Phys.Mat[,twilight.number], x=object$Spatial$Grid[,1:2],nrow=60, ncol=60),
                   col=my.golden.colors(64), main=paste("twilight number",twilight.number ))          
library(maps)
map('world', add=T)
map('state', add=T)
 
}


