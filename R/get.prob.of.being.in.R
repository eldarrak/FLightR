###############
# Theese functions extract from results probabilities of being inside a specific set of the grid points and also dates for these probabities
# first version added 11-05-2015

get.prob.of.being.in<-function(Result, Index ) {
Prob.of.being.in<-c()

nParticles<-length(inverse.rle(Result$Points.rle[[1]]))

	for (i in 1:nrow(Result$Matrix.Index.Table)) {
	Prob.of.being.in<-c(Prob.of.being.in, length(which(inverse.rle(Result$Points.rle[[i]]) %in% Index))/nParticles)
	}
	return(Prob.of.being.in)
}

find.time<-function(Prob.of.being.in, time, prob.cutoff=0.95, plot=F) {
P<-Prob.of.being.in-prob.cutoff
if (plot) {
plot(Prob.of.being.in~time, type="b", pch="+")
abline(h= prob.cutoff, col="red")
}
Transitions= which(sign(P[-length(P)]* P[-1])==-1)
if (length(Transitions) ==0 ) { 
	cat("bird has never crossed the boundary of the region!")
	} else {
		Crossing_time=time[Transitions] + (time[Transitions+1] - time[Transitions])* 
		(prob.cutoff-Prob.of.being.in[Transitions])/(Prob.of.being.in[Transitions+1]-Prob.of.being.in[Transitions])
		if (plot) abline(v=Crossing_time, col="red")
		return(Crossing_time)
	}
}


find.times.distribution<-function(Prob.of.being.in, time, plot=F) {
# let's start from 0.5 
quantiles=c(0.025, 0.25, 0.5, 0.75, 0.975)

Q1<-find.time(Prob.of.being.in, time, quantiles[1], plot=F)
Q2<-find.time(Prob.of.being.in, time, quantiles[2], plot=F)
Q3<-find.time(Prob.of.being.in, time, quantiles[3], plot=F)
Q4<-find.time(Prob.of.being.in, time, quantiles[4], plot=F)
Q5<-find.time(Prob.of.being.in, time, quantiles[5], plot=F)
# find pairs for each in Q3
Res<-data.frame(Q.025=Q1[sapply(Q3, FUN=function(x) which.min(abs(x - Q1)))], 
		Q.25=Q2[sapply(Q3, FUN=function(x) which.min(abs(x - Q2)))],
		Q.50=Q3,
		Q.75=Q4[sapply(Q3, FUN=function(x) which.min(abs(x - Q4)))],
		Q.975=Q5[sapply(Q3, FUN=function(x) which.min(abs(x - Q5)))])

	return(Res)
	}

