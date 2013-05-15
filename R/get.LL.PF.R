
get.LL.PF<-function(in.Data, All.results.mat) {
  # needed to estimate log likelihood of the optimization
  L=0
  for (i in 1:(dim(All.results.mat)[2]-1)) {
    L=L+log(mean(in.Data$Phys.Mat[All.results.mat[,i],i]))
  }
  LL=-L
  return(LL)
}
