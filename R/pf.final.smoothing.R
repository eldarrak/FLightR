
pf.final.smoothing<-function(in.Data, All.results, precision.sd=25, nParticles=1e6, save.memory=F, last.particles=NA) {
  # this function simply resamples final points proportionally to the distance to known finish.
  Final.point.real<-in.Data$stop.point
  # now we want to get distances.. I'll not index it as we will do this only once..
  Final.points.modeled=last.particles
  Weights<-dnorm(in.Data$distance[Final.points.modeled, Final.point.real], mean=0, sd=precision.sd)
  Rows<- suppressWarnings(sample.int(nParticles, replace = TRUE, prob = Weights/sum(Weights)))
  if (save.memory) {
    return(All.results[Rows])
  } else {		
    return(All.results[Rows,])
  }
}
