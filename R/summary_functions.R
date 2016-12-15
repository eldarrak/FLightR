# this is addition for FLightR 0.4.3

#' find potential stationary periods and estimates their location and movement schedule
#'
#' This function will find any sites where birds stayed longer than \code{min.stay}. Potential movement is detected by the minimum probability of movement \code{prob.cutoff}.
#' @param Result FLightR result object obtained from \code{\link{run.particle.filter}}
#' @param prob.cutoff Minimum probability that defines movement
#' @param min.stay Minimum duration of stationary period (in twilights)
#' @return list with stationary and movement statistics
#' @author Eldar Rakhimberdiev
#' @examples
#' File<-system.file("extdata", "Godwit_TAGS_format.csv", package = "FLightR")
#' # to run example fast we will cut the real data file by 2013 Aug 20
#' Proc.data<-get.tags.data(File, end.date=as.POSIXct('2013-07-02', tz='GMT'))
#' Calibration.periods<-data.frame(
#'        calibration.start=as.POSIXct(c(NA, "2014-05-05")),
#'        calibration.stop=as.POSIXct(c("2013-08-20", NA)),
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
#' # One should use 1e6 particles for the full run.
#' Result<-run.particle.filter(all.in, threads=1,
#'            nParticles=1e3, known.last=TRUE,
#'            precision.sd=25, check.outliers=FALSE)
#'
#' Summary<-stationary.migration.summary(Result, prob.cutoff=1)
#' # Use lower cut offs for real runs!
#' @export
stationary.migration.summary<-function(Result, prob.cutoff=0.1, min.stay=3) {
   Cutoffs<-which(Result$Results$Movement.results$Decision>=prob.cutoff)
   #cat('Cutoffs\n')
   #print(str(Cutoffs))
   #print(length(Cutoffs))
   #print(data.frame(start=unique(c(0, Cutoffs)), end=unique(c(Cutoffs-1, Total_length))))
   #cat('-----\n')
   if (length(Cutoffs)==0) { 
      cat('bird likely dod not move, exiting without result\n')
	  return(NULL)
   } else {
   Total_length<-length(Result$Results$Movement.results$Decision)
   Potential_stat_periods<-data.frame(start=unique(c(0, Cutoffs)), end=unique(c(Cutoffs-1, Total_length)))
   Potential_stat_periods$Duration<-Potential_stat_periods[,2]-Potential_stat_periods[,1]
   Potential_stat_periods<-Potential_stat_periods[-which(Potential_stat_periods[,3]<min.stay),]
   Quantiles<-c()
   Schedules<-c()
   for (period in 1:nrow(Potential_stat_periods)) {
      cur_period<-Potential_stat_periods[period,1]:Potential_stat_periods[period,2]
      cur_period<-cur_period[cur_period>0]
      Quantiles<-rbind(Quantiles, get_stationary_stats(Result,cur_period))
      Schedules<-rbind(Schedules, get_time_boundaries(Result, cur_period, 0.95))
   }
   ZIDist<-get_ZI_distances(Result) 
   # ok, now we look only at the distances between the periods..
   Potential_movement_periods<-data.frame(start=Potential_stat_periods$end[-nrow(Potential_stat_periods)]+1, end=Potential_stat_periods$start[-1])
   if (rev(Potential_stat_periods$end)[1]<Total_length) {
      Potential_movement_periods<-rbind(Potential_movement_periods, cbind(start=rev(Potential_stat_periods$end)[1], end=Total_length))
      Quantiles<-rbind(Quantiles, NA)
	  Quantiles$Meanlon[nrow(Quantiles)]<- rev(Result$Results$Quantiles$Meanlon)[1]
	  Quantiles$Meanlat[nrow(Quantiles)]<- rev(Result$Results$Quantiles$Meanlat)[1]
      Schedules<-rbind(Schedules, NA)
   }
   
   Distance_Moved<-c()
   for (i in 1:nrow(Potential_movement_periods)) {
      Distance_Moved<-c(Distance_Moved, sum(ZIDist[Potential_movement_periods$start[i]:Potential_movement_periods$end[i], 6]))
   }
   Quantiles$Distance.Moved<-c(0, Distance_Moved)
   Quantiles$Cumul.Distance.Moved<-cumsum(Quantiles$Distance.Moved)

   # another approach would be just to take locations and estimate movements from them
   Distance2<-c()
   for (period in 2:nrow(Quantiles)) {
      Longitudes<-c(Quantiles$Meanlon[period-1], Result$Results$Quantiles$Meanlon[Potential_movement_periods$start[period-1]:Potential_movement_periods$end[period-1]], Quantiles$Meanlon[period])
      Latitudes<-c(Quantiles$Meanlat[period-1], Result$Results$Quantiles$Meanlat[Potential_movement_periods$start[period-1]:Potential_movement_periods$end[period-1]], Quantiles$Meanlat[period])
	  Distance2<-c(Distance2, sum(sp::spDists(cbind(Longitudes, Latitudes), longlat=TRUE, segments=TRUE)))
   }
   Distance2cumulative<-c(0, cumsum(Distance2))
   Quantiles$Distance2<-c(0, Distance2)
   Quantiles$Distance2cumulative<-Distance2cumulative
   Res<-list(Stationary.periods=cbind(Quantiles, Schedules), Potential_stat_periods=Potential_stat_periods, Potential_movement_periods=Potential_movement_periods)
   # summary
   Duration<-diff(c(Result$Indices$Matrix.Index.Table$time[1], rev(Result$Indices$Matrix.Index.Table$time)[1]))
   Total_twilights<-length(Result$Indices$Matrix.Index.Table$time)
   Total_sites<-nrow(Res$Potential_stat_periods)
   
   cat('\n\n\nFLightR used',round(Duration), attr(Duration, 'units'), 'from', format(Result$Indices$Matrix.Index.Table$time[1], format='%d-%b-%Y'), 'to', format(rev(Result$Indices$Matrix.Index.Table$time)[1], format='%d-%b-%Y') , 'with', Total_twilights, 'twilights\n')
   cat('During this time tag moved approximately', round(max(Res$Stationary.periods$Distance2cumulative)), 'km with', nrow(Res$Stationary.periods), 'statiuonary periods from which', length(which(Res$Potential_stat_periods$Duration>30)), 'were longer than two weeks\n')
   cat('\n\n Detected statonary periods (without any merging!):\n')
   print(Res$Stationary.periods[,c(4,7, 12,15,23, 24, 27, 32)])
   return(Res)
   }
}


get_time_boundaries<-function(Result, twilights.index, utilisation.distribution.percentile=0.95) {
 
   Points<-get_utilisation_points(Result, twilights.index, utilisation.distribution.percentile)

   # ok now time boundaries...
 
   Schedule<-find.times.distribution(Result,Points)
 
   # Arrival   
   Cur_cross_a<-which.min(abs(Schedule$Q.50-Result$Indices$Matrix.Index.Table$time[(twilights.index)[1]]))
   Arrival<-Schedule[Cur_cross_a,]
   names(Arrival)<-paste('Arrival', names(Arrival), sep='.')
   
   if (twilights.index[1] <=1) Arrival[1,]<-NA
   
   # Departure
   Cur_cross_d<-which.min(abs(Result$Indices$Matrix.Index.Table$time[rev(twilights.index)[1]]-Schedule$Q.50))
   Departure<-Schedule[Cur_cross_d,]
   names(Departure)<-paste('Departure', names(Departure), sep='.')
   
   if (rev(twilights.index)[1] >=length(Result$Indices$Matrix.Index.Table$time)) Departure[1,]<-NA
   Res<-cbind(Arrival, Departure) 
   return(Res)
   }
 

get_stationary_stats<-function(Result, twilights.index) {

   All.Points<-get_points_distribution(Result, twilights.index)

   Mode<-t(as.matrix(Result$Spatial$Grid[which.max(All.Points),1:2]))

   # calculate stats for longitude:
   Points.Index<-inverse.rle(list(values=(1:length(All.Points)), lengths=All.Points))
   Points.Index<-sample(Points.Index, min(length(Points.Index), 1e6))
   All.longs<-Result$Spatial$Grid[Points.Index,1]
   All.lats<-Result$Spatial$Grid[Points.Index,2]

   All_tmp_1<-(Mode[1,1]-All.longs)%%360
   All_longs.shifted<-ifelse(All_tmp_1>180, All_tmp_1-360,  All_tmp_1)
   All_longs_sd<-stats::sd(All_longs.shifted)
   Sum_lons_shifted<-summary(All_longs.shifted)
   Sum_lons<-ifelse(Mode[1,1]-Sum_lons_shifted>180, Mode[1,1]-Sum_lons_shifted-360, Mode[1,1]-Sum_lons_shifted)
   CI_lon_shifted<-stats::quantile(All_longs.shifted, probs = c(0.025, 0.975))
   CI_lons<-ifelse(Mode[1,1]-CI_lon_shifted>180, Mode[1,1]-CI_lon_shifted-360, Mode[1,1]-CI_lon_shifted)
   
   Quantiles<-c(summary(All.lats), SD=stats::sd(All.lats), Mode=Mode[1,2], Sum_lons, SD=All_longs_sd, Mode=Mode[1,1], stats::quantile(All.lats, probs = c(0.025, 0.975)), CI_lons)
   Quantiles<- as.data.frame(t(Quantiles))
   names(Quantiles)[1:7]<-paste(names(Quantiles)[1:7], "lat", sep="")
   names(Quantiles)[9:15]<-paste(names(Quantiles)[9:15], "lon", sep="")

   names(Quantiles)<-gsub("\\s","", names(Quantiles))
   names(Quantiles)<-gsub("1","F", names(Quantiles))
   names(Quantiles)<-gsub("3","T", names(Quantiles))

   names(Quantiles)[17:20]<-c("LCI.lat", "UCI.lat", "LCI.lon", "UCI.lon")
   
   return(Quantiles)
   }
   
#' Estimate distances moved between twilights
#' 
#' This function estimate distances with all zeros from stationary periods. This means many of the resulting movements will have 0 as the distance
#' @param Result An object created by \code{\link{run.particle.filter}}.
#' @return a data frame containing median and quartiles for the distances and also departure and arrival time
#'
#' @author Eldar Rakhimberdiev
#' @export
get_ZI_distances<-function(Result) {   
     Distances<-Result$Results$Transitions.rle
   for (i in 1:length(Result$Results$Transitions.rle)) {
       Distances[[i]]$values<-	sapply(Result$Results$Transitions.rle[[i]]$values, FUN=function(x) dist.fun(x, Result))
    }
 
    DistancesZI<-c()
  
    for (i in 1:length(Result$Results$Transitions.rle)) {
	   Inverse<-inverse.rle(Distances[[i]])
       DistancesZI<-rbind(DistancesZI,   c(stats::quantile(Inverse, c(0.25, 0.5, 0.75)), Mean=mean(Inverse)))
    }
  
    DistancesZI<-cbind(Departure=as.POSIXct(c(NA,Result$Results$Movement.results$time[-length(Result$Results$Movement.results$time)]), origin=c('1970-01-01'), tz='UTC'), Arrival= as.POSIXct(as.numeric(Result$Results$Movement.results$time), tz='UTC', origin='1970-01-01'), as.data.frame(DistancesZI))

	return(DistancesZI)
}