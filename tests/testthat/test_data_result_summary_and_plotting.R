
test_that('plot_lon_lat_works',  {
   File<-system.file("extdata", "Godwit_TAGS_format.csv", package = "FLightR")
   Proc.data<-get.tags.data(File, end.date=as.POSIXct('2013-08-21', tz='GMT'))
   Calibration.periods<-data.frame(
        calibration.start=NA,
        calibration.stop=as.POSIXct("2013-08-20"),
		lon=5.43, lat=52.93) 
   Calibration<-make.calibration(Proc.data, Calibration.periods, likelihood.correction=FALSE)
   Grid<-make.grid(left=0, bottom=50, right=10, top=56,
     distance.from.land.allowed.to.use=c(-Inf, Inf),
     distance.from.land.allowed.to.stay=c(-Inf, Inf))
   all.in<-make.prerun.object(Proc.data, Grid, start=c(5.43, 52.93), Calibration=Calibration, threads=1)
   Result<-run.particle.filter(all.in, threads=1,
           nParticles=1e4, known.last=TRUE, check.outliers=FALSE)
   # not too clever but it looks like that other ways of graphics test are not stable..
   expect_silent(plot_lon_lat(Result))
   }
)

test_that('map_flightr_ggmap_works',  {
   File<-system.file("extdata", "Godwit_TAGS_format.csv", package = "FLightR")
   Proc.data<-get.tags.data(File, end.date=as.POSIXct('2013-08-21', tz='GMT'))
   Calibration.periods<-data.frame(
        calibration.start=NA,
        calibration.stop=as.POSIXct("2013-08-20"),
		lon=5.43, lat=52.93) 
   Calibration<-make.calibration(Proc.data, Calibration.periods, likelihood.correction=FALSE)
   Grid<-make.grid(left=0, bottom=50, right=10, top=56,
     distance.from.land.allowed.to.use=c(-Inf, Inf),
     distance.from.land.allowed.to.stay=c(-Inf, Inf))
   all.in<-make.prerun.object(Proc.data, Grid, start=c(5.43, 52.93), Calibration=Calibration, threads=1)
   Result<-run.particle.filter(all.in, threads=1,
           nParticles=1e4, known.last=TRUE, check.outliers=FALSE)
   # not too clever but it looks like that other ways of graphics test are not stable..
   p<-map.FLightR.ggmap(Result, seasonal.donut.location=NULL,return.ggobj=TRUE)    
   expect_equal(length(p), 9)
   }
)
