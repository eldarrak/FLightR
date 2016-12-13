test_that('BAStag2TAGS can return',  {
   raw_data<-data.frame(twilight=as.POSIXct(seq(as.numeric(as.POSIXct('2016-12-09 10:55'), tz='UTC'), by=300, length.out=1000), origin='1970-01-01', tz='UTC'), light=runif(1000, 0,9))
   twl<-data.frame(Twilight=as.POSIXct(seq(as.numeric(as.POSIXct('2016-12-09 10:55'), tz='UTC'), by=43200, length.out=8), origin='1970-01-01', tz='UTC'))
   twl$Rise<-rep(c(1,2),4)
   twl$Deleted=0
   twl$Deleted[1]<-1
   threshold=1.5
   tmp<-BAStag2TAGS(raw_data, twl, threshold) 
   expect_equal(dim(tmp), c(1008,5))
  }
)

test_that('BAStag2TAGS can write',  {
   raw_data<-data.frame(twilight=as.POSIXct(seq(as.numeric(as.POSIXct('2016-12-09 10:55'), tz='UTC'), by=300, length.out=1000), origin='1970-01-01', tz='UTC'), light=runif(1000, 0,9))
   twl<-data.frame(Twilight=as.POSIXct(seq(as.numeric(as.POSIXct('2016-12-09 10:55'), tz='UTC'), by=43200, length.out=8), origin='1970-01-01', tz='UTC'))
   twl$Rise<-rep(c(1,2),4)
   twl$Deleted=0
   twl$Deleted[1]<-1
   threshold=1.5
   BAStag2TAGS(raw_data, twl, threshold, 'BAStag2TAGS.csv') 
   tmp<-read.csv('BAStag2TAGS.csv')
   expect_equal(dim(tmp), c(1008,5))
  }
)

test_that('twGeos2TAGS can return',  {
   raw_data<-data.frame(twilight=as.POSIXct(seq(as.numeric(as.POSIXct('2016-12-09 10:55'), tz='UTC'), by=300, length.out=1000), origin='1970-01-01', tz='UTC'), light=runif(1000, 0,9))
   twl<-data.frame(Twilight=as.POSIXct(seq(as.numeric(as.POSIXct('2016-12-09 10:55'), tz='UTC'), by=43200, length.out=8), origin='1970-01-01', tz='UTC'))
   twl$Rise<-rep(c(1,2),4)
   twl$Deleted=0
   twl$Deleted[1]<-1
   threshold=1.5
   tmp<-twGeos2TAGS(raw_data, twl, threshold) 
   expect_equal(dim(tmp), c(1008,5))
  }
)

test_that('twGeos2TAGS can write',  {
   raw_data<-data.frame(twilight=as.POSIXct(seq(as.numeric(as.POSIXct('2016-12-09 10:55'), tz='UTC'), by=300, length.out=1000), origin='1970-01-01', tz='UTC'), light=runif(1000, 0,9))
   twl<-data.frame(Twilight=as.POSIXct(seq(as.numeric(as.POSIXct('2016-12-09 10:55'), tz='UTC'), by=43200, length.out=8), origin='1970-01-01', tz='UTC'))
   twl$Rise<-rep(c(1,2),4)
   twl$Deleted=0
   twl$Deleted[1]<-1
   threshold=1.5
   twGeos2TAGS(raw_data, twl, threshold, 'twGeos2TAGS.csv') 
   tmp<-read.csv('twGeos2TAGS.csv')
   expect_equal(dim(tmp), c(1008,5))
  }
)


test_that('GeoLight2TAGS can return',  {
   raw_data<-data.frame(twilight=as.POSIXct(seq(as.numeric(as.POSIXct('2016-12-09 10:55'), tz='UTC'), by=300, length.out=1000), origin='1970-01-01', tz='UTC'), light=runif(1000, 0,9))
   twl<-data.frame(Twilight=as.POSIXct(seq(as.numeric(as.POSIXct('2016-12-09 10:55'), tz='UTC'), by=43200, length.out=8), origin='1970-01-01', tz='UTC'))
   gl_twl<-data.frame(tFirst=twl$Twilight[-nrow(twl)], tSecond=twl$Twilight[-1], type=c(rep(c(0,1),3),0))
   threshold=1.5
   tmp<-GeoLight2TAGS(raw_data, gl_twl, threshold) 
   expect_equal(dim(tmp), c(1008,5))
  }
)

test_that('GeoLight2TAGS can write',  {
   raw_data<-data.frame(twilight=as.POSIXct(seq(as.numeric(as.POSIXct('2016-12-09 10:55'), tz='UTC'), by=300, length.out=1000), origin='1970-01-01', tz='UTC'), light=runif(1000, 0,9))
   twl<-data.frame(Twilight=as.POSIXct(seq(as.numeric(as.POSIXct('2016-12-09 10:55'), tz='UTC'), by=43200, length.out=8), origin='1970-01-01', tz='UTC'))
   gl_twl<-data.frame(tFirst=twl$Twilight[-nrow(twl)], tSecond=twl$Twilight[-1], type=c(rep(c(0,1),3),0))
   threshold=1.5
   GeoLight2TAGS(raw_data, gl_twl, threshold, 'GeoLight2TAGS.csv') 
   tmp<-read.csv('GeoLight2TAGS.csv')
   expect_equal(dim(tmp), c(1008,5))
  }
)

test_that('get.tags.data works',  {
   File<-system.file("extdata", "Godwit_TAGS_format.csv", package = "FLightR")
   Proc.data<-get.tags.data(File)
   expect_equal(nrow(Proc.data$FLightR.data$Data), 24211)
   expect_equal(nrow(Proc.data$FLightR.data$twilights), 162)
  }
)

