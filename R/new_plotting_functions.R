#' plots result over map
#'
#' plots track over map with probability cloud. Can plot only part of the track if dates are specified. Note that you can use it only after obtaining and registering in you current session Google Api Key. For details on the API key check [here](http://ornithologyexchange.org/forums/topic/38315-mapflightrggmap-error).
#' @param Result FLightR result object obtained from \code{\link{run.particle.filter}}
#' @param dates either NULL if all twilights should be included or data.frame with first column - start of the period and second end of the period. Each line represents a new period
#' @param plot.cloud Should probability cloud be plotted? If TRUE cloud is estimated by \code{\link[ggplot2]{stat_density2d}}
#' @param map.options options passed to \code{\link[ggmap]{get_map}}, note that \code{zoom} option is defined separately
#' @param plot.options plotting options. Not defined yet!
#' @param save.options options passed to \code{\link[ggplot2]{ggsave}}. Filename should be defined here.
#' @param zoom Zoom for map. If 'auto' FLightR will try to find optimal zoom level by downloading different size maps and checking whether all the points fit the map.
#' @param return.ggobj Should ggobj be returned for subsequent checks and/or replotting
#' @param seasonal.colors if true points of the track will have seasonal colors
#' @param seasonal.donut.location if NULL - no color wheel placed, otherwise select one of 'bottomleft', 'bottomright', 'topleft'
#' @param seasonal.donut.proportion how much of X axis should color wheel occupy.
#' return either NULL or ggplot2 class object
#' @param save should function save results with \code{\link[ggplot2]{ggsave}}?
#' @return if 'return.ggobj=TRUE' return ggplot object otherwise returns 'NULL'.
#' @examples
#' File<-system.file("extdata", "Godwit_TAGS_format.csv", package = "FLightR")
#' # to run example fast we will cut the real data file by 2013 Aug 20
#' Proc.data<-get.tags.data(File, end.date=as.POSIXct('2013-06-25', tz='GMT'))
#' Calibration.periods<-data.frame(
#'        calibration.start=as.POSIXct(c(NA, "2014-05-05"), tz='GMT'),
#'        calibration.stop=as.POSIXct(c("2013-08-20", NA), tz='GMT'),
#'        lon=5.43, lat=52.93) 
#'        #use c() also for the geographic coordinates, if you have more than one calibration location
#'        # (e. g.,  lon=c(5.43, 6.00), lat=c(52.93,52.94))
#' 
#' # NB Below likelihood.correction is set to FALSE for fast run! 
#' # Leave it as default TRUE for real examples
#' Calibration<-make.calibration(Proc.data, Calibration.periods, likelihood.correction=FALSE)
#' 
#' Grid<-make.grid(left=0, bottom=50, right=10, top=56,
#'   distance.from.land.allowed.to.use=c(-Inf, Inf),
#'   distance.from.land.allowed.to.stay=c(-Inf, Inf))
#'
#' all.in<-make.prerun.object(Proc.data, Grid, start=c(5.43, 52.93),
#'                              Calibration=Calibration, threads=2)
#' # here we will run only 1e4 partilces for a very short track.
#' # One should use 1e6 particles for the full run
#' Result<-run.particle.filter(all.in, threads=1,
#'            nParticles=1e3, known.last=TRUE,
#'            precision.sd=25, check.outliers=FALSE)
#'
#'\dontrun{
#' map.FLightR.ggmap(Result, seasonal.donut.location=NULL, zoom=6, save=FALSE)
#'} 
#' # for this short track without variance seasonal donut does not work,
#' # but for normall track it will.
#' @author Eldar Rakhimberdiev
#' @export map.FLightR.ggmap
map.FLightR.ggmap<-function(Result, dates=NULL, plot.cloud=TRUE, map.options=NULL, plot.options=NULL, save.options=NULL, zoom="auto", return.ggobj=FALSE, seasonal.colors=TRUE, seasonal.donut.location='topleft', seasonal.donut.proportion=0.5, save=TRUE) {
if (!is.null(plot.options)) warning("plot options are not in use yet. Let me know what you would like to have here.")

if (utils::packageVersion('ggmap')[1]<2.7) { stop('map.FLightR.ggmap function works only with ggmap >= 2.7.x') }
if (!ggmap::has_google_key()) stop('From August 2018 Google allows to use Google maps only for users with the API key, please get one and proceed as described here: http://ornithologyexchange.org/forums/topic/38315-mapflightrggmap-error/')
# dates should be a data.frame with first point - starting dates and last column end dates for periods

# ggsave.options is a list that will be will be directly passed to ggsave

Opt<-options('ggmap')
if(is.null(Opt$ggmap$google$account_type)) Opt$ggmap$google$account_type<-'standard'
if(is.null(Opt$ggmap$google$day_limit)) Opt$ggmap$google$day_limit<-2500
if(is.null(Opt$ggmap$google$second_limit)) Opt$ggmap$google$day_limit<-50
if(is.null(Opt$ggmap$google$client)) Opt$ggmap$google$client<-NA
if(is.null(Opt$ggmap$google$signature)) Opt$ggmap$google$signature<-NA
options('ggmap'= Opt$ggmap)

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
		 Points<-Result$Spatial$Grid[All.Points>0,,drop=FALSE][sample.int(length(All.Points[All.Points>0]), size = 10000, replace = TRUE, prob = All.Points[All.Points>0]), 1:2]
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
	if (!is.finite(location[1])) location[1] <-180
	if (!is.finite(location[3])) location[3] <- -180
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

	p<-ggmap::ggmap(background)#, maprange=TRUE)

	if (plot.cloud) {
	 #if (overdateline) Points[,1]<-ifelse(Points[,1]>0, Points[,1]-360, Points[,1])
	 
	 # for r cmd check:
	 ..level..<-NA
	 lon<-NA
	 lat<-NA
	 #################
	 p<-  p+ ggplot2::stat_density2d(data=data.frame(Points), ggplot2::aes(fill = ..level.., alpha = ..level.., x=lon, y=lat), size = 0.01,  geom = 'polygon', n=400) 
	 
	 p<-p+ ggplot2::scale_fill_gradient(low = "green", high = "red") 
	 p<-p+ ggplot2::scale_alpha(range = c(0.00, 0.25), guide = FALSE) 
     }

	#p<-p+ggplot2::coord_map(projection="mercator", 
    # xlim=c(attr(background, "bb")$ll.lon, attr(background, "bb")$ur.lon),
    # ylim=c(attr(background, "bb")$ll.lat, attr(background, "bb")$ur.lat))
    p<-p+ ggplot2::theme(legend.position = "none", axis.title = ggplot2::element_blank(), text = ggplot2::element_text(size = 12))
    
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
       p<-p+ ggplot2::geom_path(data=Points2plot[cur_twilights,],ggplot2::aes(x=lon,y=lat),  colour=grDevices::grey(0.3))
	   if (seasonal.colors) {
		  Seasonal_palette<-grDevices::colorRampPalette(grDevices::hsv(1-((1:365)+(365/4))%%365/365, s=0.8, v=0.8), space="Lab")
          Days<-as.integer(format(Result$Results$Quantiles$time[cur_twilights], "%j"))
		  p<-p+ggplot2::geom_point(data=Points2plot[cur_twilights,], shape=21,  colour=grDevices::grey(0.3), bg=Seasonal_palette(365)[Days])
	   } else {
		  p<-p+ggplot2::geom_point(data=Points2plot[cur_twilights,], shape="+",  colour=grDevices::grey(0.3))
	   }
	}
	
	if (!is.null(seasonal.donut.location)) {
       Xrange<-c(attr(background, "bb")$ll.lon, attr(background, "bb")$ur.lon)
       Yrange<-c(attr(background, "bb")$ll.lat, attr(background, "bb")$ur.lat)
       d<-seasonal_donut()  
       g = ggplot2::ggplotGrob(d)

    if (seasonal.donut.location=='bottomleft') {
        p<-p + ggmap::inset(grob = g, xmin = Xrange[1], xmax = Xrange[1]+seasonal.donut.proportion*(Yrange[2]-Yrange[1]), ymin = Yrange[1], ymax = Yrange[1]+seasonal.donut.proportion*(Yrange[2]-Yrange[1]))
    }

    if (seasonal.donut.location=='topleft') {
        p<-p + ggmap::inset(grob = g, xmin = Xrange[1], xmax = Xrange[1]+seasonal.donut.proportion*(Yrange[2]-Yrange[1]), ymin = Yrange[2]-seasonal.donut.proportion*(Yrange[2]-Yrange[1]), ymax = Yrange[2])
    }
	if (seasonal.donut.location=='bottomright') {
        p<-p + ggmap::inset(grob = g, xmin = Xrange[2]-seasonal.donut.proportion*(Yrange[2]-Yrange[1]), xmax = Xrange[2], ymin = Yrange[1], ymax = Yrange[1]+seasonal.donut.proportion*(Yrange[2]-Yrange[1]))
    }

    }
	graphics::plot(p)
	if (save) {
    if (is.null(save.options)) save.options<-list()
	if (is.null(save.options$filename)) save.options$filename<-"FLightR.map.pdf"
	do.call(ggplot2::ggsave, save.options)
	}
	if (return.ggobj) {return(p)} else {return(NULL)}
	}

#' plots result by longitude and latitude
#'
#' This function plots result by latitude and longitude in either vertical or horizontal layout.
#' @param Result FLightR result object obtained from \code{\link{run.particle.filter}}
#' @param scheme either 'vertical' or 'horizontal' layouts
#' return 'NULL'
#' @examples
#' File<-system.file("extdata", "Godwit_TAGS_format.csv", package = "FLightR")
#' # to run example fast we will cut the real data file by 2013 Aug 20
#' Proc.data<-get.tags.data(File, end.date=as.POSIXct('2013-07-02', tz='GMT'))
#' Calibration.periods<-data.frame(
#'        calibration.start=as.POSIXct(c(NA, "2014-05-05"), tz='GMT'),
#'        calibration.stop=as.POSIXct(c("2013-08-20", NA), tz='GMT'),
#'        lon=5.43, lat=52.93) 
#'        #use c() also for the geographic coordinates, if you have more than one calibration location
#'        # (e. g.,  lon=c(5.43, 6.00), lat=c(52.93,52.94))
#' 
#' # NB Below likelihood.correction is set to FALSE for fast run! 
#' # Leave it as default TRUE for real examples
#' Calibration<-make.calibration(Proc.data, Calibration.periods, likelihood.correction=FALSE)
#' 
#' Grid<-make.grid(left=0, bottom=50, right=10, top=56,
#'   distance.from.land.allowed.to.use=c(-Inf, Inf),
#'   distance.from.land.allowed.to.stay=c(-Inf, Inf))
#'
#' all.in<-make.prerun.object(Proc.data, Grid, start=c(5.43, 52.93),
#'                              Calibration=Calibration, threads=2)
#' # here we will run only 1e4 partilces for a very short track.
#' # One should use 1e6 particles for the full run
#' Result<-run.particle.filter(all.in, threads=1,
#'            nParticles=1e3, known.last=TRUE,
#'            precision.sd=25, check.outliers=FALSE)
#'
#' plot_lon_lat(Result)
#'
#' @author Eldar Rakhimberdiev
#' @export plot_lon_lat
plot_lon_lat<-function(Result, scheme=c("vertical", "horizontal")) {
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
      graphics::par(mfrow=c(1,2)) 
   } else {
      graphics::par(mfrow=c(2,1)) 
   }
 
   graphics::par(mar=c(2,4,3,1),cex=1)
   suppressWarnings(Sys.setlocale("LC_ALL", "English")) 

   #------here I have create the equinox dates..
   Years<-unique(format(Quantiles$time, format="%Y"))
   eq<-c(as.POSIXct(paste(Years, "-09-22 00:00:00 GMT", sep="")), as.POSIXct(paste(Years, "-03-20 12:00:00 GMT", sep="")))
   eq<-eq[eq>min(Quantiles$time) & eq<max(Quantiles$time)]

   #------- vert grid
   Vert_grid<-seq(as.POSIXct("2000-01-01", tz='GMT'), as.POSIXct("2030-01-01", tz='GMT'), by="month")
   Vert_grid<-Vert_grid[Vert_grid>=min(Quantiles$time) & Vert_grid<=max(Quantiles$time)]
  
   #Longitude
   graphics::plot(Quantiles$Medianlon~Quantiles$time, las=1,col=grDevices::grey(0.1),pch=16,
      ylab="Longitude",xlab="",lwd=2, ylim=range(c( Quantiles$LCI.lon,
      Quantiles$UCI.lon )), type="n", axes=FALSE)
   graphics::axis(2, las=1)
   graphics::axis.POSIXct(1, x=Quantiles$time,  format="1-%b")
   graphics::box()

   # add vertical lines for the first day of every month
   
   graphics::abline(v=Vert_grid, col=grDevices::grey(0.5), lty=2)
   graphics::abline(h=seq(-180, 360, by=10), col=grDevices::grey(0.5), lty=2)
   
   graphics::abline(v=eq, col=2, lwd=2, lty=1)
  

   graphics::polygon(x=c(Quantiles$time, rev(Quantiles$time)), 
      y=c(Quantiles$LCI.lon, rev(Quantiles$UCI.lon)),
      col=grDevices::grey(0.9), border=grDevices::grey(0.5))

   graphics::polygon(x=c(Quantiles$time, rev(Quantiles$time)),
   y=c(Quantiles$TrdQu.lon, rev(Quantiles$FstQu.lon)),
   col=grDevices::grey(0.7), border=grDevices::grey(0.5))

   graphics::lines(Quantiles$Medianlon~Quantiles$time, col=grDevices::grey(0.1),lwd=2)

   #Latitude
   graphics::par(mar=c(3,4,1,1))
   
   graphics::plot(Quantiles$Medianlat~Quantiles$time, las=1,col=grDevices::grey(0.1),
     pch=16,ylab="Latitude",xlab="",lwd=2,
     ylim=range(c( Quantiles$UCI.lat, Quantiles$LCI.lat )), type="n", axes=FALSE)
   graphics::axis(2, las=1)
   graphics::axis.POSIXct(1, x=Quantiles$time,  format="1-%b")
   graphics::box()
   graphics::abline(v=Vert_grid, col=grDevices::grey(0.5), lty=2)
   graphics::abline(h=seq(-80, 80, by=10), col=grDevices::grey(0.5), lty=2)

   graphics::abline(v=eq, col=2, lwd=2, lty=1)

   graphics::polygon(x=c(Quantiles$time, rev(Quantiles$time)), y=c(Quantiles$LCI.lat, rev(Quantiles$UCI.lat)),
           col=grDevices::grey(0.9), border=grDevices::grey(0.5))

   graphics::polygon(x=c(Quantiles$time, rev(Quantiles$time)), y=c(Quantiles$TrdQu.lat, rev(Quantiles$FstQu.lat)),
           col=grDevices::grey(0.7), border=grDevices::grey(0.5))

   graphics::lines(Quantiles$Medianlat~Quantiles$time, col=grDevices::grey(0.1),lwd=2)
   return(NULL)
}


get_buffer<-function(coords, r){
  p = sp::SpatialPoints(matrix(coords, ncol=2), proj4string=sp::CRS("+proj=longlat +datum=WGS84"))
  aeqd <- sprintf("+proj=aeqd +lat_0=%s +lon_0=%s +x_0=0 +y_0=0",
                    p@coords[[2]], p@coords[[1]])
  projected <- sp::spTransform(p, sp::CRS(aeqd))
  buffered <- rgeos::gBuffer(projected, width=r, byid=TRUE)
  buffer_lonlat <- sp::spTransform(buffered, CRS=p@proj4string)
  return(buffer_lonlat)
  }
  
get_gunion_r<-function(Result) {
    Distances=  sp::spDists(Result$Spatial$Grid[1:min(c(nrow(Result$Spatial$Grid), 1000)),1:2], longlat=TRUE)
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

  Points_selected<-get_utilisation_points(Result, twilights.index, percentile)
  
 if (is.null(r)) {
     r=get_gunion_r(Result)
     }
  #Now I want to create a spatial buffers around points with this radius
  
  #spoints = SpatialPoints(Result$Spatial$Grid[Points_selected,1:2], proj4string= CRS("+proj=longlat +datum=WGS84"))
  
  Buffers<-apply(matrix(Result$Spatial$Grid[Points_selected,1:2], ncol=2), 1, get_buffer, r=r)
   Buff_comb<-Buffers[[1]] 
   if (length(Buffers)>1) {
   for (i in 2:length(Buffers)) {
       Buff_comb<-rgeos::gUnion(Buff_comb, Buffers[[i]])
  }
  }
  Buff_comb_simpl<-rgeos::gSimplify(Buff_comb, tol=0.01, topologyPreserve=TRUE)
  return(list(Buffer=Buff_comb_simpl, nPoints=length(Points)))
  }

#' plots resulting track over map with uncertainty shown by space utilisation distribution
#' 
#' May be use not only for the whole track but for a set of specific dates, e.g. to show spatial uncertainty during migration. Note that you can use it only after obtaining and registering in you current session Google Api Key. For details on the API key check [here](http://ornithologyexchange.org/forums/topic/38315-mapflightrggmap-error).
#' 
#' @param Result FLightR result object obtained from \code{\link{run.particle.filter}}
#' @param dates Use NULL if all twilights will be used for plotting, one integer if specific twilight should be plotted (line number in Result$Results$Quantiles). Use data.frame with first column - start of the period and second - end of the period and each line represents a new period to plot specific periods, e.g. wintering or migration.
#' @param map.options options passed to \code{\link[ggmap]{get_map}}, note that \code{zoom} option is defined separately
#' @param percentiles Probability breaks for utilisation distribution
#' @param zoom Zoom for map. If 'auto' FLightR will try to find optimal zoom level by downloading different size maps and checking whether all the points fit the map.
#' @param geom_polygon.options options passed to \code{\link[ggplot2]{geom_polygon}}
#' @param save.options options passed to \code{\link[ggplot2]{ggsave}}. Filename should be defined here.
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
#' @examples
#' File<-system.file("extdata", "Godwit_TAGS_format.csv", package = "FLightR")
#' # to run example fast we will cut the real data file by 2013 Aug 20
#' Proc.data<-get.tags.data(File, end.date=as.POSIXct('2013-06-25', tz='GMT'))
#' Calibration.periods<-data.frame(
#'        calibration.start=as.POSIXct(c(NA, "2014-05-05"), tz='GMT'),
#'        calibration.stop=as.POSIXct(c("2013-08-20", NA), tz='GMT'),
#'        lon=5.43, lat=52.93) 
#'        #use c() also for the geographic coordinates, if you have more than one calibration location
#'        # (e. g.,  lon=c(5.43, 6.00), lat=c(52.93,52.94))
#' 
#' # NB Below likelihood.correction is set to FALSE for fast run! 
#' # Leave it as default TRUE for real examples
#' Calibration<-make.calibration(Proc.data, Calibration.periods, likelihood.correction=FALSE)
#' 
#' Grid<-make.grid(left=0, bottom=50, right=10, top=56,
#'   distance.from.land.allowed.to.use=c(-Inf, Inf),
#'   distance.from.land.allowed.to.stay=c(-Inf, Inf))
#'
#' all.in<-make.prerun.object(Proc.data, Grid, start=c(5.43, 52.93),
#'                              Calibration=Calibration, threads=1)
#' # here we will run only 1e4 partilces for a very short track.
#' # One should use 1e6 particles for the full run
#' Result<-run.particle.filter(all.in, threads=1,
#'            nParticles=1e3, known.last=TRUE,
#'            precision.sd=25, check.outliers=FALSE)
#'
#'\dontrun{
#' plot_util_distr(Result, zoom=6, save=FALSE)
#'}
#'
#' @author Eldar Rakhimberdiev
#' @export plot_util_distr
plot_util_distr<-function(Result, dates=NULL, map.options=NULL, percentiles=c(0.4, 0.6, 0.8), zoom="auto", geom_polygon.options=NULL, save.options=NULL, color.palette=NULL, use.palette=TRUE, background=NULL, plot=TRUE, save=TRUE, add.scale.bar=FALSE, scalebar.options=NULL) {

if (utils::packageVersion('ggmap')[1]<2.7) {stop('plot_util_distr function works only with ggmap >= 2.7.x')}

if (!ggmap::has_google_key()) stop('From August 2018 Google allows to use Google maps only for users with the API key, please get one and proceed as described here: http://ornithologyexchange.org/forums/topic/38315-mapflightrggmap-error/')


Opt<-options('ggmap')
if(is.null(Opt$ggmap$google$account_type)) Opt$ggmap$google$account_type<-'standard'
if(is.null(Opt$ggmap$google$day_limit)) Opt$ggmap$google$day_limit<-2500
if(is.null(Opt$ggmap$google$second_limit)) Opt$ggmap$google$day_limit<-50
if(is.null(Opt$ggmap$google$client)) Opt$ggmap$google$client<-NA
if(is.null(Opt$ggmap$google$signature)) Opt$ggmap$google$signature<-NA
options('ggmap'= Opt$ggmap)


   # for r cmd check
   long<-NA
   lat<-NA
   group<-NA
   ## end ##########
   
if (is.null(color.palette)) color.palette <- grDevices::colorRampPalette(rev(c("#edf8fb",
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
#Index_not_dupl<-(1:length(percentiles))[!duplicated(nPoints)]

#-----------------#
# and now we want to plot it

#if (is.null(background)) {

# combine location with map options
	if (is.null(map.options)) map.options<-list()
	#if (is.null(map.options$location)) map.options$location<-location
	if (is.null(map.options$zoom) & is.numeric(zoom)) map.options$zoom=zoom
	if (is.null(map.options$col)) map.options$col="bw"
	location<-res_buffers[[length(res_buffers)]]@bbox

    if (is.null(map.options$location)) map.options$location<-rowMeans(res_buffers[[length(res_buffers)]]@bbox)

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

    #}
    geom_polygon.options.external=geom_polygon.options
	
	p<-ggmap::ggmap(background, maprange=TRUE)
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
	
	geom_polygon.options$mapping=ggplot2::aes(x=long, y=lat, group=group)
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
	
    my.rbind.SpatialPolygons = function(..., makeUniqueIDs = FALSE) {
       dots = list(...)
       names(dots) <- NULL
       stopifnot(sp::identicalCRS(dots))
       # checkIDSclash(dots)
       pl = do.call(c, lapply(dots, function(x) methods::slot(x, "polygons")))
       if (makeUniqueIDs)
               pl = makeUniqueIDs(pl)
       sp::SpatialPolygons(pl, proj4string = sp::CRS(sp::proj4string(dots[[1]])))
	}
    makeUniqueIDs <- function(lst) {
	   ids = sapply(lst, function(i) methods::slot(i, "ID"))
	   if (any(duplicated(ids))) {
		  ids <- make.unique(as.character(unlist(ids)), sep = "")
		  for (i in seq(along = ids))
			lst[[i]]@ID = ids[i]
	   }
	   lst
    }
	
	b<-do.call(my.rbind.SpatialPolygons,  c(res_buffers, list(makeUniqueIDs=TRUE))) 
    SPDF = sp::SpatialPolygonsDataFrame(b, data.frame(percentile = percentiles, row.names=names(b)))
	return(list(res_buffers=SPDF, p=p, bg=background))
}


 
seasonal_donut<-function() {
   Seasonal_palette<-grDevices::colorRampPalette(grDevices::hsv(1-((1:365)+(365/4))%%365/365, s=0.8, v=0.8), space="Lab")
   pie.data<-data.frame(group=as.factor(c(1:12)), value=rep(1,12))  
   pie.data$fraction = pie.data$value / sum(pie.data$value)
   pie.data$ymax = cumsum(pie.data$fraction)
   pie.data$ymin = c(0, utils::head(pie.data$ymax, n = -1))
   
   donut<-ggplot2::ggplot(data = pie.data, ggplot2::aes(fill = pie.data$group, ymax = pie.data$ymax, ymin = pie.data$ymin, xmax = 1, xmin = 2))  +
      ggplot2::geom_rect(colour = "grey30", show.legend = FALSE) +
      ggplot2::coord_polar(theta = "y", start=pi) +
      ggplot2::xlim(c(0, 2)) +
      ggplot2::theme_bw() + 
	  ggplot2::theme(plot.background = ggplot2::element_rect(fill = "transparent" ,colour = NA))+
	  ggplot2::theme(panel.background = ggplot2::element_blank())+
	  ggplot2::theme(panel.grid=ggplot2::element_blank()) +
      ggplot2::theme(axis.text=ggplot2::element_blank()) +
      ggplot2::theme(axis.ticks=ggplot2::element_blank()) +
      ggplot2::theme(panel.border=ggplot2::element_blank()) +
	  ggplot2::scale_fill_manual(values=c(Seasonal_palette(13)))+
      ggplot2::geom_text(ggplot2::aes(x = 1.5, y = (( pie.data$ymin+ pie.data$ymax)/2), label =  pie.data$group), colour=grDevices::grey(0.99), size=5) +
      ggplot2::xlab("") +
      ggplot2::ylab("") 
	  return(donut)
}

#' plot likelihood surface over map
#'
#' plots specific likelihood surface over map
#' @details function plots likelihoods before particle filter run, so these are pure results of calibrations without any movement model
#' @param object either output from \code{\link{make.prerun.object}} or \code{\link{run.particle.filter}}
#' @param date either NULL or a date (possibly with time) closest to the twilight you wan to be plotted
#' @param twilight.index number of likelihood surface to be plotted 
#' @return 'NULL'
#' @examples
#' File<-system.file("extdata", "Godwit_TAGS_format.csv", package = "FLightR")
#' # to run example fast we will cut the real data file by 2013 Aug 20
#' Proc.data<-get.tags.data(File, end.date=as.POSIXct('2013-07-02', tz='GMT'))
#' Calibration.periods<-data.frame(
#'        calibration.start=as.POSIXct(c(NA, "2014-05-05"), tz='GMT'),
#'        calibration.stop=as.POSIXct(c("2013-08-20", NA), tz='GMT'),
#'        lon=5.43, lat=52.93) 
#'        #use c() also for the geographic coordinates, if you have more than one calibration location
#'        # (e. g.,  lon=c(5.43, 6.00), lat=c(52.93,52.94))
#' 
#' # NB Below likelihood.correction is set to FALSE for fast run! 
#' # Leave it as default TRUE for real examples
#' Calibration<-make.calibration(Proc.data, Calibration.periods, likelihood.correction=FALSE)
#' 
#' Grid<-make.grid(left=0, bottom=50, right=10, top=56,
#'   distance.from.land.allowed.to.use=c(-Inf, Inf),
#'   distance.from.land.allowed.to.stay=c(-Inf, Inf))
#'
#' all.in<-make.prerun.object(Proc.data, Grid, start=c(5.43, 52.93),
#'                              Calibration=Calibration, threads=2)
#' plot_likelihood(all.in, twilight.index=10)
#'
#' @author Eldar Rakhimberdiev
#' @export plot_likelihood
plot_likelihood<-function(object, date=NULL, twilight.index=NULL) {
   my.golden.colors <- grDevices::colorRampPalette(c("white","#FF7100"))
   if (!is.null(date)) {
      date<-as.POSIXct(date, tz='GMT')
      twilight.index<-which.min(abs(object$Indices$Matrix.Index.Table$time-date))
	}

	  
   fields::image.plot(fields::as.image(object$Spatial$Phys.Mat[,twilight.index], x=object$Spatial$Grid[,1:2],nrow=60, ncol=60),
                   col=my.golden.colors(64), main=paste("twilight number",twilight.index, format(object$Indices$Matrix.Index.Table$time[twilight.index], tz='UTC')))     
   wrld_simpl<-NA				   
   load(system.file("data", "wrld_simpl.rda", package = "maptools"))
   sp::plot(wrld_simpl, add=TRUE)
   return(NULL)
}

get_points_distribution<-function(Result, twilights.index) {
    Points_rle<-Result$Results$Points.rle[twilights.index]
	All.Points<-rep(0, nrow(Result$Spatial$Grid))
	for (twilight in 1:length(twilights.index)) {
        All.Points[Points_rle[[twilight]]$values]<-All.Points[Points_rle[[twilight]]$values] + Points_rle[[twilight]]$lengths
    }
	return(All.Points)	 
}


get_utilisation_points<-function(Result, twilights.index, percentile) { 
  
  All.Points<-get_points_distribution(Result, twilights.index)
  Order<-order(All.Points, decreasing=TRUE)
 
  Order_rle<-rle(All.Points[Order]) # 
  
  nParticles<- sum(Result$Results$Points.rle[[1]]$lengths)
  
  Cum_probs<-cumsum(Order_rle$values/nParticles/length(twilights.index))
  
  tmp<-which(Cum_probs<=percentile)
  if (length(tmp)==0) { 
     Selected_rle<-1 # if there is no such bin - return the most likely one
  } else {
     Selected_rle<-max(which(Cum_probs<=percentile))
  }
  Points<-Order[1:sum(Order_rle$lengths[1:Selected_rle])]

  Points_selected<-(1:length(All.Points))[Points]
  return(Points_selected) 
}
 
