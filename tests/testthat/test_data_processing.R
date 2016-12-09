
test_that('plot_slopes_by_location works',  {
   File<-system.file("extdata", "Godwit_TAGS_format.csv", package = "FLightR")
   Proc.data<-get.tags.data(File)
   png('plot_slopes_by_location.png')
   plot_slopes_by_location(Proc.data=Proc.data, location=c(5.43, 52.93))
   abline(v=as.POSIXct("2013-08-20")) # end of first calibration period
   abline(v=as.POSIXct("2014-05-05")) # start of the second calibration period
   dev.off()
   expect_equal(89, round(file.size('plot_slopes_by_location.png')/100))
   }
)
