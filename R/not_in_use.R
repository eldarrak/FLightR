# not_in_use.R
# this functions are currently not used but may be used at some point....
# I will not release them on CRAN and will archive locally.

get.time.shift<-function(start, Twilight.time.mat.Calib.dawn, Twilight.log.light.mat.Calib.dawn, Twilight.time.mat.Calib.dusk, Twilight.log.light.mat.Calib.dusk,  log.light.borders=NA, diap=c(-600, 600), plot=T,  log.irrad.borders=c(-15, 50), verbose=F) {

	# THIS FUNCTION IS NOT CURRENTLY U
	# this version will try to make slope difference insignificant..
	# will try to make this criteria as the main.

	# for now - very simple function:
	# 3 loops - minutes - 10 seconds - seconds..
	# 	# we may want to do it with optim..
	x<-list(start[1], start[2])
	get.time.shift.internal<-function(start, Twilight.time.mat.Calib.dawn, Twilight.log.light.mat.Calib.dawn, Twilight.time.mat.Calib.dusk, Twilight.log.light.mat.Calib.dusk,  log.light.borders, diap, Step=1, plot,  log.irrad.borders=c(-8, 1.5)) {
		Final<-c()
		for (i in seq(diap[1], diap[2], Step)) {
#	cat("i", i, "Step", Step, "\n")	
			Calibration<-try(logger.template.calibration(Twilight.time.mat.Calib.dawn+i, Twilight.log.light.mat.Calib.dawn, Twilight.time.mat.Calib.dusk+i, Twilight.log.light.mat.Calib.dusk, positions=start, log.light.borders=log.light.borders,  log.irrad.borders= log.irrad.borders))
			
			# and now I want to estimate LL that the bird was in the known location
			LL<-0
			All<- Inf
			
			# now combining all data together
			Twilight.time.mat.Calib<-cbind(Twilight.time.mat.Calib.dawn, Twilight.time.mat.Calib.dusk)
			Dusk<-c(rep(FALSE, dim(Twilight.time.mat.Calib.dawn)[2]), rep(TRUE, dim(Twilight.time.mat.Calib.dusk)[2]))
			Twilight.log.light.mat.Calib<-cbind(Twilight.log.light.mat.Calib.dawn, Twilight.log.light.mat.Calib.dusk)
			
			for (Twilight.ID in 1:(dim(Twilight.time.mat.Calib)[2])) {
			#print(Twilight.ID )
			Twilight.solar.vector<-solar(as.POSIXct(Twilight.time.mat.Calib[c(1:24, 26:49), Twilight.ID]+i, tz="gmt", origin="1970-01-01"))
			Twilight.log.light.vector<-Twilight.log.light.mat.Calib[c(1:24, 26:49), Twilight.ID]
			Res<-get.current.slope.prob(x, calibration=Calibration,  Twilight.solar.vector=Twilight.solar.vector, Twilight.log.light.vector=Twilight.log.light.vector, plot=F, verbose=verbose,  log.light.borders=log.light.borders,  log.irrad.borders= log.irrad.borders, dusk=Dusk[Twilight.ID])
			#cat("Twilight.ID", Twilight.ID, "Res", Res, "\n")
			if (is.na(Res) | Res==0 ) {
			Res<-get.current.slope.prob(x, calibration=Calibration,  Twilight.solar.vector=Twilight.solar.vector, Twilight.log.light.vector=Twilight.log.light.vector, plot=F, verbose=T,  log.light.borders=log.light.borders,  log.irrad.borders= log.irrad.borders, dusk=Dusk[Twilight.ID])
			
			break

			} 
			LL<-LL+log(Res)
			All<-c(All, Res)
		}
#print(length(All))
#print(All)
#print((dim(Twilight.time.mat.Calib)[2]))		
		if (length(All)<=(dim(Twilight.time.mat.Calib)[2])) {
		All=-Inf;
		LL=-Inf
		}
		Final=rbind(Final, cbind(i, LL, min(All),  Calibration$Significance.of.dusk.dawn.diff))
		if (TRUE) print(Final)
		}
		return(Final)
		}
		
	# minutes
	Final.min<-get.time.shift.internal(start, Twilight.time.mat.Calib.dawn, Twilight.log.light.mat.Calib.dawn, Twilight.time.mat.Calib.dusk, Twilight.log.light.mat.Calib.dusk,  log.light.borders, diap=diap, Step=60, plot=F, log.irrad.borders=log.irrad.borders)
	cat("\n")
	if (max(Final.min[,2])==-Inf) stop("something went wrong\n more likely that the logger was not stable during the calibration time\n exclude some of the twilights and try again\n")
	#Diap.10.sec<-c(Final.min[which.max(Final.min[,4]),1]-60, Final.min[which.max(Final.min[,4]),1]+60)
	Diap.10.sec<-c(Final.min[which.max(Final.min[,2]),1]-60, Final.min[which.max(Final.min[,2]),1]+60)
	Final.10secs.min<-get.time.shift.internal(start, Twilight.time.mat.Calib.dawn, Twilight.log.light.mat.Calib.dawn, Twilight.time.mat.Calib.dusk, Twilight.log.light.mat.Calib.dusk,  log.light.borders, diap=Diap.10.sec, Step=10,  log.irrad.borders= log.irrad.borders)
	
	Diap.1.sec<-c(Final.10secs.min[which.max(Final.10secs.min[,2]),1]-10, Final.10secs.min[which.max(Final.10secs.min[,2]),1]+10)
	Final.1sec.min<-get.time.shift.internal(start, Twilight.time.mat.Calib.dawn, Twilight.log.light.mat.Calib.dawn, Twilight.time.mat.Calib.dusk, Twilight.log.light.mat.Calib.dusk,  log.light.borders, diap=Diap.1.sec, Step=1,  log.irrad.borders= log.irrad.borders)
	
	return(Final.1sec.min[which.max(Final.1sec.min[,2]),1])
	}
	
