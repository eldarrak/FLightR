test_that('plot_slopes_by_location works',  {
   File<-system.file("extdata", "Godwit_TAGS_format.csv", package = "FLightR")
   Proc.data<-get.tags.data(File)
   expect_null(plot_slopes_by_location(Proc.data=Proc.data, location=c(5.43, 52.93)))
   expect_silent(abline(v=as.POSIXct("2013-08-20"))) # end of first calibration period
   expect_silent(abline(v=as.POSIXct("2014-05-05"))) # start of the second calibration period
   }
)

test_that('get.tags.data works',  {
   File<-system.file("extdata", "Godwit_TAGS_format.csv", package = "FLightR")
   Proc.data<-get.tags.data(File)
   Calibration.periods<-data.frame(
        calibration.start=as.POSIXct(c(NA, "2014-05-05")),
        calibration.stop=as.POSIXct(c("2013-08-20", NA)),
        lon=5.43, lat=52.93) 
   Calibration<-make.calibration(Proc.data, Calibration.periods, likelihood.correction=FALSE)
   expect_equal(0.15, round(Calibration$Parameters$LogSlope[1],2))
   expect_equal(0.1, round(Calibration$Parameters$LogSlope[2],1))
   expect_equal('parametric.slope', Calibration$Parameters$calibration.type)
   expect_equal(c(1.5,9), Calibration$Parameters$log.light.borders)
   expect_equal(c(-3,3), Calibration$Parameters$log.irrad.borders)
   expect_false(Calibration$Parameters$impute.on.boundaries)
   expect_null(Calibration$Parameters$c_fun)
   }
)

test_that('find.stationary.location works',  {
   File<-system.file("extdata", "Godwit_TAGS_format.csv", package = "FLightR")
   Proc.data<-get.tags.data(File)
   Location<-find.stationary.location(Proc.data, '2013-07-20', '2013-08-20', initial.coords=c(10, 50), reltol=0.03, plot=FALSE,  print.optimization = FALSE)
   expect_equal(c(5,55), Location)
   }
)

test_that('make.grid works',  {
   Grid<-make.grid(left=0, bottom=50, right=10, top=56,
     distance.from.land.allowed.to.use=c(-100, Inf),
     distance.from.land.allowed.to.stay=c(0, 100))
   expect_equal(c(145, 3), dim(Grid))
   expect_equal(111.5, sum(Grid[,3]))
   }
)

test_that('parallel setup works works',  {
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
   all.in<-make.prerun.object(Proc.data, Grid, start=c(5.43, 52.93), Calibration=Calibration, threads=4)
   expect_equal(dim(all.in$Spatial$Phys.Mat), c(180, 126))
   expect_true(max(all.in$Spatial$Phys.Mat)>1)
}
)
