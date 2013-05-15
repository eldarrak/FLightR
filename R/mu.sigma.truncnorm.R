# there are 2 functions but only simplest is exported into the package
mu.sigma.truncnorm<-function(x, a=45, b=500) {
  if (length(unique(x))>1) {
    tr.norm<-function(prm) {
      sum(-log(truncnorm:::dtruncnorm(x,a=a,b=b,mean=prm[1],sd=prm[2])))
    }
		Res=try(optim(c(mean(x),sd(x)), tr.norm, method="BFGS"))
		#if (class(Res)=="try-error") Res=try(optim(c(mean(x),sd(x)), tr.norm,method="SANN"))
		if (class(Res)=="try-error") save(x, Res, file="x.RData")
    return(c(Res$par[1], Res$par[2]))
  } else {
    return(c(mean(x), sd(x)))
  }
}