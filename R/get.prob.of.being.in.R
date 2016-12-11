###############
# These functions extract from results probabilities of being inside a specific set of the grid points and also dates for these probabilities
# first version added 11-05-2015

get.prob.of.being.in<-function(Result, Index ) {
Prob.of.being.in<-c()

nParticles<-length(inverse.rle(Result$Results$Points.rle[[1]]))

	for (i in 1:nrow(Result$Indices$Matrix.Index.Table)) {
	Prob.of.being.in<-c(Prob.of.being.in, length(which(inverse.rle(Result$Results$Points.rle[[i]]) %in% Index))/nParticles)
	}
	return(Prob.of.being.in)
}

find.time<-function(Prob.of.being.in, time, prob.cutoff=0.95, plot=FALSE) {
P<-Prob.of.being.in-prob.cutoff
if (plot) {
graphics::plot(Prob.of.being.in~time, type="b", pch="+")
graphics::abline(h= prob.cutoff, col="red")
}
Transitions= which(sign(P[-length(P)]* P[-1])==-1)
if (length(Transitions) ==0 ) { 
	#warning("bird has never crossed the boundary of the region!")
	    return(NA)
	} else {
		Crossing_time=time[Transitions] + (time[Transitions+1] - time[Transitions])* 
		(prob.cutoff-Prob.of.being.in[Transitions])/(Prob.of.being.in[Transitions+1]-Prob.of.being.in[Transitions])
		if (plot) graphics::abline(v=Crossing_time, col="red")
		return(Crossing_time)
	}
}

#' extracts times of arrival and departure to/from spatial extent
#'
#' Idea of this functions is to extract schedules for known location
#' @param Result FLightR result object obtained from \code{\link{run.particle.filter}}
#' @param Spatial.Index Row numbers for spatial grid (\code{Result$Spatial$Grid}) to estimate schedules for.
#' @return dataframe with columns for 0.025, 0.25, 0.5, 0.75, 0.975 probability of line crossing and rows for every crossing.
#' @author Eldar Rakhimberdiev
#' @export
find.times.distribution<-function(Result, Spatial.Index) {
   Prob.of.being.in<-get.prob.of.being.in(Result, Spatial.Index)
   time<-Result$Indices$Matrix.Index.Table$time
   # let's start from 0.5 
   quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975)

   Q1<-find.time(Prob.of.being.in, time, quantiles[1], plot=FALSE)
   Q2<-find.time(Prob.of.being.in, time, quantiles[2], plot=FALSE)
   Q3<-find.time(Prob.of.being.in, time, quantiles[3], plot=FALSE)
   Q4<-find.time(Prob.of.being.in, time, quantiles[4], plot=FALSE)
   Q5<-find.time(Prob.of.being.in, time, quantiles[5], plot=FALSE)
   Q1.5<-sort(c(Q1, Q5))
   Q2.4<-sort(c(Q2, Q4))
   
   if (is.na(Q3[1])) {
      return(NA)
   } else {
   Q.50<-as.POSIXct(Q3,origin='1970-01-01', tz="UTC")
      
   Res<-data.frame(Q.025=as.POSIXct(NA, origin='1970-01-01', tz="UTC"),
                   Q.25=as.POSIXct(NA, origin='1970-01-01', tz="UTC"),
                   Q.50=Q.50,
				   Q.75=as.POSIXct(NA, origin='1970-01-01', tz="UTC"),
				   Q.975=as.POSIXct(NA, origin='1970-01-01', tz="UTC"))

   for (i in 1:length(Q.50)) {
      Q1_cur<-Q1.5[Q1.5<=Q3[i]]
	  if (i>1) Q1_cur<-Q1_cur[Q1_cur>Q3[i-1]]
      if (length(Q1_cur)>0) {
	     Res$Q.025[i]<-Q1_cur[length(Q1_cur)]
      }
	  
      Q2_cur<-Q2.4[Q2.4<=Q3[i]]
	  if (i>1) Q2_cur<-Q2_cur[Q2_cur>Q3[i-1]]
      if (length(Q2_cur)>0) {
	     Res$Q.25[i]<-Q2_cur[length(Q2_cur)]
      }	  

	  Q4_cur<-Q2.4[Q2.4>=Q3[i]]
	  if (i<length(Q.50)) Q4_cur<-Q4_cur[Q4_cur<Q3[i+1]]
      if (length(Q4_cur)>0) {
	     Res$Q.75[i]<-Q4_cur[1]
      }	
	  
 	  Q5_cur<-Q1.5[Q1.5>=Q3[i]]
	  if (i<length(Q.50)) Q5_cur<-Q5_cur[Q5_cur<Q3[i+1]]
      if (length(Q5_cur)>0) {
	  
	     Res$Q.975[i]<-Q5_cur[1]
      }	 
   } 
   return(Res)
   }
}
