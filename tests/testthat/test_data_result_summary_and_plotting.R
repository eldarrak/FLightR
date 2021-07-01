
test_that('plot_lon_lat_works',  {
   #skip_on_os('solaris')
   File<-system.file("extdata", "Godwit_TAGS_format.csv", package = "FLightR")
   Proc.data<-get.tags.data(File, end.date=as.POSIXct('2013-10-02', tz='GMT'))
   Calibration.periods<-data.frame(
        calibration.start=NA,
        calibration.stop=as.POSIXct("2013-08-20", tz='GMT'),
		lon=5.43, lat=52.93) 
   Calibration<-make.calibration(Proc.data, Calibration.periods, likelihood.correction=FALSE)
   Grid<-make.grid(left=0, bottom=50, right=10, top=56,
     distance.from.land.allowed.to.use=c(-Inf, Inf),
     distance.from.land.allowed.to.stay=c(-Inf, Inf))
   all.in<-make.prerun.object(Proc.data, Grid, start=c(5.43, 52.93), Calibration=Calibration, threads=1)
   Result<-run.particle.filter(all.in, threads=1,
           nParticles=1e3, known.last=TRUE, check.outliers=FALSE)
   # not too clever but it looks like that other ways of graphics test are not stable..
   expect_silent(plot_lon_lat(Result))
   }
)

test_that('map_flightr_ggmap_works',  {
   File<-system.file("extdata", "Godwit_TAGS_format.csv", package = "FLightR")
   Proc.data<-get.tags.data(File, end.date=as.POSIXct('2013-10-02', tz='GMT'))
   Calibration.periods<-data.frame(
        calibration.start=NA,
        calibration.stop=as.POSIXct("2013-08-20", tz='GMT'),
        lon=5.43, lat=52.93) 
   Calibration<-make.calibration(Proc.data, Calibration.periods, likelihood.correction=FALSE)
   Grid<-make.grid(left=0, bottom=50, right=10, top=56,
      distance.from.land.allowed.to.use=c(-Inf, Inf),
      distance.from.land.allowed.to.stay=c(-Inf, Inf))
   all.in<-make.prerun.object(Proc.data, Grid, start=c(5.43, 52.93), Calibration=Calibration, threads=1)
   Result<-run.particle.filter(all.in, threads=1,
           nParticles=1e3, known.last=TRUE, check.outliers=FALSE)
   if (utils::packageVersion('ggmap')[1]>=2.7) {
   if (ggmap::has_google_key()) {
         # not too clever but it looks like that other ways of graphics test are not stable..
         p<-map.FLightR.ggmap(Result, seasonal.donut.location=NULL,return.ggobj=TRUE, zoom=5, save=FALSE)    
         expect_equal(length(p), 9)
   } else {
      expect_error({
         p<-map.FLightR.ggmap(Result, seasonal.donut.location=NULL,return.ggobj=TRUE, zoom=5, save=FALSE)    
      })
      }
      } else {
      expect_error({
         p<-map.FLightR.ggmap(Result, seasonal.donut.location=NULL,return.ggobj=TRUE, zoom=5, save=FALSE)    
      })
   } 
   }
 )


 test_that('plot_util_distr_works',  {
   File<-system.file("extdata", "Godwit_TAGS_format.csv", package = "FLightR")
   Proc.data<-get.tags.data(File, end.date=as.POSIXct('2013-10-02', tz='GMT'))
   Calibration.periods<-data.frame(
        calibration.start=NA,
        calibration.stop=as.POSIXct("2013-08-20", tz='GMT'),
		lon=5.43, lat=52.93) 
   Calibration<-make.calibration(Proc.data, Calibration.periods, likelihood.correction=FALSE)
   Grid<-make.grid(left=0, bottom=50, right=10, top=56,
     distance.from.land.allowed.to.use=c(-Inf, Inf),
     distance.from.land.allowed.to.stay=c(-Inf, Inf))
   all.in<-make.prerun.object(Proc.data, Grid, start=c(5.43, 52.93), Calibration=Calibration, threads=1)
   Result<-run.particle.filter(all.in, threads=1,
           nParticles=1e3, known.last=TRUE, check.outliers=FALSE)
   if (utils::packageVersion('ggmap')[1]>=2.7) {
   if (ggmap::has_google_key()) {
   # not too clever but it looks like that other ways of graphics test are not stable..
   expect_output(plot_util_distr(Result, zoom=5, save=FALSE))
   } else {
      expect_error(plot_util_distr(Result, zoom=5, save=FALSE))
   }
   } else {
      expect_error({
         p<-map.FLightR.ggmap(Result, seasonal.donut.location=NULL,return.ggobj=TRUE, zoom=5, save=FALSE)    
      })
   } 
}   
)


 test_that('summary produces null foor stationary tag',  {
    File<-system.file("extdata", "Godwit_TAGS_format.csv", package = "FLightR")
    Proc.data<-get.tags.data(File, end.date=as.POSIXct('2013-07-02', tz='GMT'))
    Calibration.periods<-data.frame(
         calibration.start=NA,
         calibration.stop=as.POSIXct("2013-08-20", tz='GMT'),
		 lon=5.43, lat=52.93) 
    Calibration<-make.calibration(Proc.data, Calibration.periods, likelihood.correction=FALSE)
    Grid<-make.grid(left=0, bottom=50, right=10, top=56,
      distance.from.land.allowed.to.use=c(-Inf, Inf),
      distance.from.land.allowed.to.stay=c(-Inf, Inf))
    all.in<-make.prerun.object(Proc.data, Grid, start=c(5.43, 52.93), Calibration=Calibration, threads=1)
    Result<-run.particle.filter(all.in, threads=1,
            nParticles=1e3, known.last=TRUE, check.outliers=FALSE)
    expect_null(stationary.migration.summary(Result, prob.cutoff=1))
    }
 )
