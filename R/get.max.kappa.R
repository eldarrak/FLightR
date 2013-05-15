## function to get optimization coefficient...
get.max.kappa<-function(max.sd, a=45, b=500, maxIter=10000, verbose=F) {
  require(truncnorm)
  # max(dtruncnorm())/min(dtruncnorm())==max(dvonmises/min(dvonmises))
  try.k<-100
  Run=T
  Iter=1
  while (Iter<maxIter & Run) {
    max.diff.dist<-truncnorm:::dtruncnorm(x=a, mean=a, sd=max.sd, a=a, b=b)/truncnorm:::dtruncnorm(x=b, mean=a, sd=max.sd, a=a, b=b)
    max.diff.dir<-suppressWarnings(dvonmises(x=0, mu=0, kappa=try.k)/dvonmises(x=pi, mu=0, kappa=try.k))
    if (verbose) {
      cat("Iteration:", Iter,  "\n")
      cat("max.diff.dist:", max.diff.dist, "\n")
      cat("max.diff.dir:", max.diff.dir, "\n")
      cat("ratio:", max.diff.dist/max.diff.dir, "\n") }
    Ratio=max.diff.dist/max.diff.dir
    if (abs(log(Ratio))>0.5) {
      if (Ratio>1) { try.k=try.k+0.5		}
      else {try.k=try.k-0.5}
      if (verbose) cat("try.k", try.k, "\n")
    } else {
      cat("OPTIMIZED! in", Iter, "iterations\n")
      Run=F
    }
    Iter=Iter+1
  }
  return(try.k)
}
# get.max.kappa(max.sd=25, a=1, b=500, maxIter=1000)

