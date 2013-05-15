
pf.par.internal<-function(x, Current.Proposal) {
  # this simple function is needed to save I/0 during parallel run. Being executed at a slave it creates inoput for the next function taking Current proposal from the slave and Parameters from the master
  new.Parameters<-c(x=list(x), get("Parameters"), Current.Proposal=list(Current.Proposal))
  Res<-do.call(generate.points.dirs, new.Parameters)
  return(Res)
}