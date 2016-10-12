# this is addition for FLightR 0.4.3

#' find potential stationary periods and estimates their location and movement schedule
#'
#' This function will find any sites where birds stayed longer than \code{min.stay}. Potential movement is detected by the minimum probability of movement \code{prob.cutoff}.
#' @param Result FLightR result object obtained from \code{\link{run.particle.filter}}
#' @param prob.cutoff Minimum probability that defines movement
#' @param min.stay Minimum duration of stationary period (in twilights)
#' @return list with stationary and movement statistics
#' @author Eldar Rakhimberdiev
#' @export
stationary.migration.summary<-function(Result, prob.cutoff=0.1, min.stay=3) {
   Cutoffs<-which(Result$Results$Movement.results$Decision>=prob.cutoff)
   Total_length<-length(Result$Results$Movement.results$Decision)
   Potential_stat_periods<-data.frame(start=unique(c(0, Cutoffs)), end=unique(c(Cutoffs-1, Total_length)))
   Potential_stat_periods$Duration<-Potential_stat_periods[,2]-Potential_stat_periods[,1]
   Potential_stat_periods<-Potential_stat_periods[-which(Potential_stat_periods[,3]<min.stay),]
   Quantiles<-c()
   Schedules<-c()
   for (period in 1:nrow(Potential_stat_periods)) {
      cur_period<-Potential_stat_periods[period,1]:Potential_stat_periods[period,2]
      cur_period<-cur_period[cur_period>0]
      Quantiles<-rbind(Quantiles,get_stationary_stats(Result,cur_period))
      Schedules<-rbind(Schedules, get_time_boundaries(Result, cur_period, 0.95))
   }

   ZIDist<-get_ZI_distances(Result) 

   # ok, now we look only at the distances between the periods..
   Potential_movement_periods<-data.frame(start=Potential_stat_periods$end[-nrow(Potential_stat_periods)]+1, end=Potential_stat_periods$start[-1])
   Distance_Moved<-c()
   for (i in 1:nrow(Potential_movement_periods)) {
      Distance_Moved<-c(Distance_Moved, sum(ZIDist[Potential_movement_periods$start[i]:Potential_stat_periods$end[i], 6]))
   }
   Quantiles$Distance.Moved<-c(0, Distance_Moved)
   Quantiles$Cumul.Distance.Moved<-cumsum(Quantiles$Distance.Moved)

   Res<-list(Stationary.periods=cbind(Quantiles, Schedules), Potential_stat_periods=Potential_stat_periods, Potential_movement_periods=Potential_movement_periods)
   # summary
   Duration<-diff(c(Result$Indices$Matrix.Index.Table$time[1], rev(Result$Indices$Matrix.Index.Table$time)[1]))
   Total_twilights<-length(Result$Indices$Matrix.Index.Table$time)
   Total_sites<-nrow(Res$Potential_stat_periods)
   
   cat('\n\n\nFLightR used',round(Duration), attr(Duration, 'units'), 'from', format(Result$Indices$Matrix.Index.Table$time[1], format='%d-%b-%Y'), 'to', format(rev(Result$Indices$Matrix.Index.Table$time)[1], format='%d-%b-%Y') , 'with', Total_twilights, 'twilights\n')
   cat('During this time tag moved approximately', round(max(Res$Stationary.periods$Cumul.Distance.Moved)), 'km with', nrow(Res$Stationary.periods), 'statiuonary periods from which', length(which(Res$Potential_stat_periods$Duration>30)), 'were longer than two weeks\n')
   cat('\n\n Detected statonary periods (without any merging!):\n')
   print(Res$Stationary.periods[,c(4,7, 12,15,21, 22, 24, 30)])
   return(Res)
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
   All_longs_sd<-sd(All_longs.shifted)
   Sum_lons_shifted<-summary(All_longs.shifted)
   Sum_lons<-ifelse(Mode[1,1]-Sum_lons_shifted>180, Mode[1,1]-Sum_lons_shifted-360, Mode[1,1]-Sum_lons_shifted)
   CI_lon_shifted<-quantile(All_longs.shifted, probs = c(0.025, 0.975))
   CI_lons<-ifelse(Mode[1,1]-CI_lon_shifted>180, Mode[1,1]-CI_lon_shifted-360, Mode[1,1]-CI_lon_shifted)
   
   Quantiles<-c(summary(All.lats), SD=sd(All.lats), Mode=Mode[1,2], Sum_lons, SD=All_longs_sd, Mode=Mode[1,1], quantile(All.lats, probs = c(0.025, 0.975)), CI_lons)
   Quantiles<- as.data.frame(t(Quantiles))
   names(Quantiles)[1:7]<-paste(names(Quantiles)[1:7], "lat", sep="")
   names(Quantiles)[9:15]<-paste(names(Quantiles)[9:15], "lon", sep="")

   names(Quantiles)<-gsub("\\s","", names(Quantiles))
   names(Quantiles)<-gsub("1","F", names(Quantiles))
   names(Quantiles)<-gsub("3","T", names(Quantiles))

   names(Quantiles)[17:20]<-c("LCI.lat", "UCI.lat", "LCI.lon", "UCI.lon")
   
   return(Quantiles)
   }
   
