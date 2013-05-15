
create.spatial.sun.matrices<-function(in.Data, cpus=7, cluster.type="SOCK", existing.cluster=NULL) {	
	List<-suppressWarnings(split(in.Data$Matrix.Index.Table, 1:nrow(in.Data$Matrix.Index.Table)))
	Points.Land<-in.Data$Points.Land
	if (cpus >1) {
		if (is.null(existing.cluster)) {
		require("parallel")
		mycl <- parallel:::makeCluster(cpus, type=cluster.type)
		parallel:::clusterSetRNGStream(mycl)
		### we don' need to send all parameters to node. so keep it easy..
		} else {mycl<-existing.cluster}
		WD<-getwd()
		parallel:::clusterExport(mycl, "WD", envir=environment())
		parallel:::clusterEvalQ(mycl, setwd(WD)) 
		parallel:::clusterExport(mycl,"Points.Land", envir=environment())
		#parallel:::clusterEvalQ(mycl, library("tripEstimation"))
		parallel:::clusterExport(mycl,c("sun.matrix.internal","node.run", "elevation", "solar"))
		Res<-simplify2array(parLapply(cl = mycl, List, node.run))
		if (is.null(existing.cluster)) stopCluster(mycl)
	} else {
	Res<-simplify2array(lapply(List, FUN=function(x) sun.matrix.internal(x, Points.Land)))
	}
	return(Res)
}
