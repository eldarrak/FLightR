
generate.points.dirs<-function(x , in.Data, Current.Proposal, a=45, b=500) {
  # this function is needed to generate new points - it works as from input point Index and biological proposal
  ################
  # x has 3 columns
  # Index
  # number
  # nMoving
  
  resample <- function(x, ...) x[sample.int(length(x), ...)]
  
  if (x[[3]]>0) {
    Dists.distr<-in.Data$distance[x[[1]],]
    Dists.probs<-truncnorm:::dtruncnorm(Dists.distr, a=a, b=b, Current.Proposal$M.mean, Current.Proposal$M.sd)
    ###
    #  library(fields)
    # dists
    #image.plot(as.image(Dists.probs, x=in.Data$Points.Land, nrow=50, ncol=50))
    #
    if (Current.Proposal$Kappa>0) {
      Angles.dir<-in.Data$Azimuths[x[[1]],]
      Angles.probs<-as.numeric(suppressWarnings(circular:::dvonmises(Angles.dir/180*pi, mu=Current.Proposal$Direction/180*pi, kappa=Current.Proposal$Kappa)))
      Angles.probs[is.na(Angles.probs)]<-0
    }
    else {
      Angles.probs<-0.1591549
    }
    # this is neede to catch an error 
    if (which.max(c(min(c(Angles.probs, Dists.probs)), 0))==2) {
      cat("PF produced weird probs \n")
      tmp.out<-list(Angles.probs=Angles.probs, Dists.probs=Dists.probs, Current.Proposal=Current.Proposal, Data=x)
      save(tmp.out, file=paste("tmp.out", round(runif(n=1,min=1, max=1000)), "RData", sep="."))
      Biol.proposal<-pmax.int(Dists.probs*Angles.probs, 0)
    } else {
      Biol.proposal<-Dists.probs*Angles.probs
    }
    pos.biol<-suppressWarnings(sample.int(length(Biol.proposal), size=x[[3]], replace=T, prob=Biol.proposal))
    # ok, now we want to return more complicated stuff - final indexes!
    return(resample(as.integer(c(pos.biol, rep(x[[1]],(x[[2]]-x[[3]]))))))
  }
  else {
    return(as.integer(rep(x[[1]],(x[[2]]))))
  }
}

