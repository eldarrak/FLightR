
update.proposal.PF<-function(output.matrix, Trans, in.Data, save.rle=F, fixed.parameters=NA, a=45, b=500, parallel=F, existing.cluster=NULL) {
 mycl=existing.cluster
  # main function for proposal update
  # from the version 1.5 this function will do estimation of distance and it's SD..
  # save.rle will save the rle object that could be used for futher calculations...
  
  # fixed.parameters could be list(Kappa=1, M.sd=250)
  require("truncnorm")
  require("circular")
  require("CircStats")
  #=========================================
  # here is another place to change...
  # ver. 1.8
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #Trans<-vector(mode = "list", length = (dim(output.matrix)[2]-1))
  #cat("   extracting indexes\n")
  #for (i in 1:(dim(output.matrix)[2]-1)) {
  #  Trans[[i]]<-get.transition.rle(output.matrix[,i], output.matrix[,i+1])
  #}
  cat("   estimating distances\n")
  #####   let's try to get distance distribution:
  Distances<-Trans
  dist.fun<-function(x) {
    in.Data$distance[x%/%1e5, x%%1e5]
  }
   for (i in 1:length(Trans)) {
    Distances[[i]]$values<-sapply(Trans[[i]]$values, FUN=function(x) dist.fun(x))
  }
    cat("   estimating directions\n")
  #ok, now we want to get directions
  
  Directions<-Trans
  direct.fun<-function(x) {
    in.Data$Azimuths[x%/%1e5, x%%1e5]
  }
  for (i in 1:length(Trans)) {
    Directions[[i]]$values<-sapply(Trans[[i]]$values, FUN=function(x) direct.fun(x))
  }
  cat("   estimating mean directions\n")
  Mean.Directions<-unlist(lapply(Directions, FUN=function(x) CircStats:::circ.mean(inverse.rle(list(lengths=x$lengths[!is.na(x$values)], values=x$values[!is.na(x$values)]))*pi/180)*180/pi))
  #plot(Mean.Directions)
  cat("   estimating kappas\n")
  Kappas<-unlist(lapply(Directions, FUN=function(x) CircStats:::est.kappa(inverse.rle(list(lengths=x$lengths[!is.na(x$values)], values=x$values[!is.na(x$values)]))*pi/180)))
  
  
  # now we want to get mean ditance
  cat("   estimating mean dists\n")
  #
  # this takes all the RAM.. something should be changed in the parallel version..
  #if (parallel) {
	#Mean.and.Sigma<-parLapply(cl = mycl, Distances, fun=function(x) mu.sigma.truncnorm(inverse.rle(list(lengths=x$lengths[x$values!=0], values=x$values[x$values!=0])), a=a, b=b))
  #} else {
	Mean.and.Sigma<-lapply(Distances, FUN=function(x) mu.sigma.truncnorm(inverse.rle(list(lengths=x$lengths[x$values!=0], values=x$values[x$values!=0])), a=a, b=b))
  #}
  Mean.Dists<-sapply(Mean.and.Sigma, "[[", i=1)
  cat("   estimating dists SD\n")
  Mean.SD<-sapply(Mean.and.Sigma, "[[", i=2)
  #for (i in 1:length(Trans)) {
    #============================================================
  #  Distances[[i]]$values<-sapply(Trans[[i]]$values, FUN=function(x) dist.fun(x))
	# prepare data for new mu.sigma.truncnorm
	# points.current, dists, azimuths, phys.proposal, Points.Land, point_ID, a=45, b=500, Current.mean, Current.sd
  #points.current<-output.matrix[,i+1]
  # now I need to get closest point
  #=====================================================
  # here is the place to change in v 1.8 we need to go for Trans as it will be created outside...
  # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #dists.real<-inverse.rle(Distances[[i]])
  #nonZero.Points<-which(dists.real !=0)
  #point_ID<-which.min(spDistsN1(in.Data$Points.Land, cbind(in.Data$Final.Means[i,"CENTRE.x"], in.Data$Final.Means[i,"CENTRE.y"]),longlat=T))
  #dists<-in.Data$distance[point_ID,] # these are distances from central point
  #azimuths<-in.Data$Azimuths[point_ID,]
  #phys.proposal<-in.Data$Phys.Mat[,i]
  #Points.Land<-in.Data$Points.Land
  #Current.mean<-in.Data$Matrix.Index.Table$M.mean[i]
  #Current.sd<-in.Data$Matrix.Index.Table$M.sd[i]
  #Mean.and.Sigma<-rbind(Mean.and.Sigma, mu.sigma.truncnorm(points.current, dists, dists.real ,azimuths, phys.proposal, Points.Land, point_ID, a=45, b=500, Current.mean, Current.sd, nonZero.Points))
  #Mean.and.Sigma<-rbind(Mean.and.Sigma, mu.sigma.truncnorm(points.current, dists, dists.real ,azimuths, phys.proposal, Points.Land, point_ID, a=45, b=500, Current.mean, Current.sd, nonZero.Points))
  #}
  
  #Mean.Dists<-Mean.and.Sigma[,1]
 
  #cat("   estimating SD dists SD\n")
  #Mean.SD<-Mean.and.Sigma[,2]
  
  #cat("   estimating mean and SD to report dists SD\n")
#This is temporal change!!!
# attempt to go just for simple distance estimation...
  Mean2report<-unlist(lapply(Distances, FUN=function(x) mean(inverse.rle(list(lengths=x$lengths[x$values!=0], values=x$values[x$values!=0])))))
  SD2report<-unlist(lapply(Distances, FUN=function(x) sd(inverse.rle(list(lengths=x$lengths[x$values!=0], values=x$values[x$values!=0])))))
  cat("   estimating probs of migration\n")
  # ok, now we want to get parameters for distances
  #
  Probability.of.migration<-unlist(lapply(Distances, FUN=function(x) sum(x$lengths[x$values!=0])))/(dim(output.matrix)[1])
  
  ## 
  
  #unlist(lapply(Distances, FUN=function(x) mean(inverse.rle(list(lengths=x$lengths[x$values!=0], values=x$values[x$values!=0])))))
  #plot(Mean.Dists)
  cat("   estimating median dists\n")
  Median.Dists<-unlist(lapply(Distances, FUN=function(x) median(inverse.rle(list(lengths=x$lengths[x$values!=0], values=x$values[x$values!=0])))))
   #unlist(lapply(Distances, FUN=function(x) sd(inverse.rle(list(lengths=x$lengths[x$values!=0], values=x$values[x$values!=0])))))
  #plot(Kappas)
  cat("   creating output")
  ################
  # create new proposal from the posteriors
  #New.Matrix.Index.Table<-data.frame(Direction=Mean.Directions, M.mean=Mean.Dists, M.sd=Mean.SD, Decision=Probability.of.migration, Kappa=Kappas)
  ################
  
  all.arrays.object<-in.Data
  
  all.arrays.object$Matrix.Index.Table$Decision<-Probability.of.migration
  all.arrays.object$Matrix.Index.Table$Direction<-Mean.Directions
  all.arrays.object$Matrix.Index.Table$Kappa<-Kappas
  all.arrays.object$Matrix.Index.Table$M.mean<-Mean.Dists
  all.arrays.object$Matrix.Index.Table$M.sd<-Mean.SD
  all.arrays.object$Matrix.Index.Table$M.medians<-Median.Dists
  all.arrays.object$Matrix.Index.Table$Mean2report<-Mean2report
  all.arrays.object$Matrix.Index.Table$SD2report<-SD2report
  
  #all.arrays.object$Matrix.Index.Table<-New.Matrix.Index.Table
  
  all.arrays.object$Main.Index$Biol.Prev=1:(nrow(all.arrays.object$Matrix.Index.Table))
  if(save.rle) all.arrays.object$Transitions.rle<-Trans
  
  
  if (is.list(fixed.parameters)) {
    # this part can be moved upper not to estimate parametrs in case they are fixed.
    if ("M.sd" %in% names(fixed.parameters)) all.arrays.object$Matrix.Index.Table$M.sd<-fixed.parameters$M.sd
    
    if ("M.mean" %in% names(fixed.parameters)) all.arrays.object$Matrix.Index.Table$M.mean<-fixed.parameters$M.mean
    
    if ("Kappa" %in% names(fixed.parameters)) all.arrays.object$Matrix.Index.Table$Kappa<-fixed.parameters$Kappa
    
    if ("Direction" %in% names(fixed.parameters)) all.arrays.object$Matrix.Index.Table$Kappa<-fixed.parameters$Direction
  }
  #all.arrays.object$Main.Index$Biol.Next=2:(nrow(all.arrays.object$Matrix.Index.Table))
  cat(" DONE!\n")
  return(all.arrays.object)
}