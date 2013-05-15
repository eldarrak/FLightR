
sun.matrix.internal<-function(Matrix.Index.Table.Row, Points.Land) {
	# this function will preestimate SpSunM for one row of in.data
	# spatial.Points is an Rdata file that contains XXX?? object
	#require(tripEstimation)
	require(LightR)
	require(maptools)
	require(fields)
	
	correct.estimates<-function(Latitudes) {
		# this function estimates correction coef. for twilight matrices
		# it will look for a file "Correction.coefs.RData"
		#Res<-suppressWarnings(try(load("Correction.coefs.RData"), silent=TRUE))
		data(A)
		#if (class(Res)=="try-error") stop("Correction file: \"Correction.coefs.RData\" is needed in the working directory\n")
		if (any(abs(Latitudes)>84)) warning("there is no correction coefs for latitudes over 84 degree - NAs were returned\n")
		return(predict(A, Latitudes))
	}
	
	time.angle.integrate<-function(x) {
	  sum(dnorm(elevation(x[[1]], x[[2]], solar(p.time.gmt)), mean=Angle[1], sd=Angle[2])*p.time.dens*p.dtime)
	}
	# integrating without the last point
	p.num = 100
	p.time.cumprob <- seq(0.001, 0.999, length.out=p.num)                     # in probabilities
	p.time.points  <- qnorm(p.time.cumprob,as.numeric(Matrix.Index.Table.Row$Real.time),Matrix.Index.Table.Row$Loess.se.fit*3600*sqrt(Matrix.Index.Table.Row$Loess.n))  # in time
	p.dtime        <- (p.time.points[-1] - p.time.points[-p.num])           # step in time - dTime 
	p.time.points1 <-p.time.points[-p.num]                                  # deleting last point
	p.time.dens    <- dnorm(p.time.points1,as.numeric(Matrix.Index.Table.Row$Real.time),Matrix.Index.Table.Row$Loess.se.fit*3600*sqrt(Matrix.Index.Table.Row$Loess.n))  # probability density - p(Time)
	p.time.gmt <- as.POSIXct(p.time.points1, tz="GMT", origin="1970-01-01")
	Angle<-as.matrix(Matrix.Index.Table.Row[c("Angle", "Angle.sd")])[1,]
	Res<-apply(Points.Land, 1, FUN=time.angle.integrate)
	Correction.coefs<-correct.estimates(Points.Land[,2])
	Res<-Res/Correction.coefs
	return(Res)
}

