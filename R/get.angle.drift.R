
get.angle.drift<-function(in.Data, plot=F) {
	#require(tripEstimation)
	 # last.point.fun<-function(x) {
	#	in.Data$Points.Land[x%%1e5,]
	  #}
	  Mean.angle<-NULL
	  SD.angle<-NULL
	  Angles.rle<-in.Data$Points.rle
	  for (i in 1:length(in.Data$Matrix.Index.Table$Real.time)) {
		Angles.rle[[i+1]]$values<-apply(matrix(in.Data$Points.Land[in.Data$Points.rle[[i+1]]$values,1:2], ncol=2), 1, FUN=function(x) elevation(x[[1]], x[[2]], solar(in.Data$Matrix.Index.Table$Real.time[i])))
		SD.angle<-c(SD.angle, sd(inverse.rle(Angles.rle[[i+1]])))
		Mean.angle<-c(Mean.angle, weighted.mean(Angles.rle[[i+1]]$values, Angles.rle[[i+1]]$lengths))
		}
	#Angles.rle<-Angles.rle[-1]
	if (plot) {
	plot(Mean.angle[in.Data$Matrix.Index.Table$Dusk]~as.numeric(in.Data$Matrix.Index.Table$Real.time[in.Data$Matrix.Index.Table$Dusk]), col="blue")
	lines(predict(loess(Mean.angle[in.Data$Matrix.Index.Table$Dusk]~as.numeric(in.Data$Matrix.Index.Table$Real.time[in.Data$Matrix.Index.Table$Dusk]), span=150/length(Mean.angle[in.Data$Matrix.Index.Table$Dusk])))~ as.numeric(in.Data$Matrix.Index.Table$Real.time[in.Data$Matrix.Index.Table$Dusk]), col="blue", type="l", lwd=2)
	
	points(Mean.angle[!in.Data$Matrix.Index.Table$Dusk]~as.numeric(in.Data$Matrix.Index.Table$Real.time[!in.Data$Matrix.Index.Table$Dusk]), pch="+", col="red")
	lines(predict(loess(Mean.angle[!in.Data$Matrix.Index.Table$Dusk]~as.numeric(in.Data$Matrix.Index.Table$Real.time[!in.Data$Matrix.Index.Table$Dusk]), span=150/length(Mean.angle[!in.Data$Matrix.Index.Table$Dusk])))~ as.numeric(in.Data$Matrix.Index.Table$Real.time[!in.Data$Matrix.Index.Table$Dusk]), col="red", type="l", lwd=2)
	}
	
	Mean.Dusk.angle.drift<-predict(loess(Mean.angle[in.Data$Matrix.Index.Table$Dusk]~as.numeric(in.Data$Matrix.Index.Table$Real.time[in.Data$Matrix.Index.Table$Dusk]), span=150/length(Mean.angle[in.Data$Matrix.Index.Table$Dusk])))
	SD.Dusk.angle.drift<-sd(Mean.angle[in.Data$Matrix.Index.Table$Dusk] - Mean.Dusk.angle.drift)
	Mean.Dawn.angle.drift<-predict(loess(Mean.angle[!in.Data$Matrix.Index.Table$Dusk]~as.numeric(in.Data$Matrix.Index.Table$Real.time[!in.Data$Matrix.Index.Table$Dusk]), span=150/length(Mean.angle[in.Data$Matrix.Index.Table$Dusk])))
	SD.Dawn.angle.drift<-sd(Mean.angle[!in.Data$Matrix.Index.Table$Dusk] - Mean.Dawn.angle.drift)
	# now we want to comibine them
	Mean.angle.drift<-rep(NA, length(Mean.angle))
	SD.angle.drift<-rep(NA, length(Mean.angle))
	Mean.angle.drift[in.Data$Matrix.Index.Table$Dusk]<-Mean.Dusk.angle.drift
	Mean.angle.drift[!in.Data$Matrix.Index.Table$Dusk]<-Mean.Dawn.angle.drift
	SD.angle.drift[in.Data$Matrix.Index.Table$Dusk]<-SD.Dusk.angle.drift
	SD.angle.drift[!in.Data$Matrix.Index.Table$Dusk]<-SD.Dawn.angle.drift
	return(data.frame(Mean=Mean.angle.drift, SD=SD.angle.drift))
	}
