get.initial.calibration<-function(in.Data, days.stable=list(Dusk=c(1:20), Dawn=c(1:20))) {
	par(mfrow=c(2,2))
	Dusk.angle<-my.calibrate(in.Data$Matrix.Index.Table$Real.time[in.Data$Matrix.Index.Table$Dusk==T][days.stable$Dusk],  in.Data$Points.Land[in.Data$start.point,1],  in.Data$Points.Land[in.Data$start.point,2])
	Dawn.angle<-my.calibrate(in.Data$Matrix.Index.Table$Real.time[in.Data$Matrix.Index.Table$Dusk==F][days.stable$Dawn]	,  in.Data$Points.Land[in.Data$start.point,1],  in.Data$Points.Land[in.Data$start.point,2])
	in.Data$Matrix.Index.Table$Angle[in.Data$Matrix.Index.Table$Dusk==T]<-Dusk.angle[1]
	in.Data$Matrix.Index.Table$Angle.sd[in.Data$Matrix.Index.Table$Dusk==T]<-Dusk.angle[2]
	in.Data$Matrix.Index.Table$Angle[in.Data$Matrix.Index.Table$Dusk==F]<-Dawn.angle[1]
	in.Data$Matrix.Index.Table$Angle.sd[in.Data$Matrix.Index.Table$Dusk==F]<-Dawn.angle[2]
	par(mfrow=c(1,1))
	return(in.Data$Matrix.Index.Table)
}

