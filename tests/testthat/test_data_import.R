test_that('BAStag2TAGS can return',  {
   raw_data<-data.frame(twilight=as.POSIXct(seq(as.numeric(as.POSIXct('2016-12-09 10:55', tz='GMT')), by=300, length.out=1000), origin='1970-01-01', tz='GMT'), light=runif(1000, 0,9))
   twl<-data.frame(Twilight=as.POSIXct(seq(as.numeric(as.POSIXct('2016-12-09 10:55', tz='GMT')), by=43200, length.out=8), origin='1970-01-01', tz='GMT'))
   twl$Rise<-rep(c(1,2),4)
   twl$Deleted=0
   twl$Deleted[1]<-1
   threshold=1.5
   tmp<-BAStag2TAGS(raw_data, twl, threshold) 
   expect_equal(dim(tmp), c(1008,5))
  }
)

test_that('BAStag2TAGS can write',  {
   raw_data<-data.frame(twilight=as.POSIXct(seq(as.numeric(as.POSIXct('2016-12-09 10:55', tz='GMT')), by=300, length.out=1000), origin='1970-01-01', tz='GMT'), light=runif(1000, 0,9))
   twl<-data.frame(Twilight=as.POSIXct(seq(as.numeric(as.POSIXct('2016-12-09 10:55', tz='GMT')), by=43200, length.out=8), origin='1970-01-01', tz='GMT'))
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
   raw_data<-data.frame(twilight=as.POSIXct(seq(as.numeric(as.POSIXct('2016-12-09 10:55', tz='GMT')), by=300, length.out=1000), origin='1970-01-01', tz='GMT'), light=runif(1000, 0,9))
   twl<-data.frame(Twilight=as.POSIXct(seq(as.numeric(as.POSIXct('2016-12-09 10:55', tz='GMT')), by=43200, length.out=8), origin='1970-01-01', tz='GMT'))
   twl$Rise<-rep(c(1,2),4)
   twl$Deleted=0
   twl$Deleted[1]<-1
   threshold=1.5
   tmp<-twGeos2TAGS(raw_data, twl, threshold) 
   expect_equal(dim(tmp), c(1008,5))
  }
)

test_that('twGeos2TAGS can write',  {
   raw_data<-data.frame(twilight=as.POSIXct(seq(as.numeric(as.POSIXct('2016-12-09 10:55', tz='GMT')), by=300, length.out=1000), origin='1970-01-01', tz='GMT'), light=runif(1000, 0,9))
   twl<-data.frame(Twilight=as.POSIXct(seq(as.numeric(as.POSIXct('2016-12-09 10:55', tz='GMT')), by=43200, length.out=8), origin='1970-01-01', tz='GMT'))
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
   raw_data<-data.frame(twilight=as.POSIXct(seq(as.numeric(as.POSIXct('2016-12-09 10:55', tz='GMT')), by=300, length.out=1000), origin='1970-01-01', tz='GMT'), light=runif(1000, 0,9))
   twl<-data.frame(Twilight=as.POSIXct(seq(as.numeric(as.POSIXct('2016-12-09 10:55', tz='GMT')), by=43200, length.out=8), origin='1970-01-01', tz='GMT'))
   gl_twl<-data.frame(tFirst=twl$Twilight[-nrow(twl)], tSecond=twl$Twilight[-1], type=c(rep(c(0,1),3),0))
   threshold=1.5
   tmp<-GeoLight2TAGS(raw_data, gl_twl, threshold) 
   expect_equal(dim(tmp), c(1008,5))
  }
)

test_that('GeoLight2TAGS can write',  {
   raw_data<-data.frame(twilight=as.POSIXct(seq(as.numeric(as.POSIXct('2016-12-09 10:55', tz='GMT')), by=300, length.out=1000), origin='1970-01-01', tz='GMT'), light=runif(1000, 0,9))
   twl<-data.frame(Twilight=as.POSIXct(seq(as.numeric(as.POSIXct('2016-12-09 10:55', tz='GMT')), by=43200, length.out=8), origin='1970-01-01', tz='GMT'))
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

write_fake_tags_file<-function(light, file=tempfile(fileext=".csv")) {
   datetime<-as.POSIXct("2020-01-01 00:00:00", tz="GMT") + seq_along(light) * 300
   twilight<-rep(0, length(light))
   twilight[c(200, 400, 600, 800)]<-c(1, 2, 1, 2)
   interp<-rep(FALSE, length(light))
   interp[1]<-TRUE
   out<-data.frame(
      datetime=format(datetime, "%Y-%m-%dT%H:%M:%S.000Z", tz="GMT"),
      light=light,
      twilight=twilight,
      interp=interp,
      excluded=FALSE
   )
   utils::write.csv(out, file=file, row.names=FALSE, quote=FALSE)
   file
}

test_that('explicit raw logger mode bypasses fragile max-light detection', {
   File<-write_fake_tags_file(rep(seq(4, 500, length.out=1000)))
   Proc.data<-get.tags.data(File, tag.type="mk", light.scale="raw")
   expect_equal(Proc.data$tagtype, "mk")
   expect_false(Proc.data$log_transformed)
   expect_equal(Proc.data$log.light.borders, log(c(4, 61)))
   expect_equal(max(Proc.data$FLightR.data$Data$light), 500)
})

test_that('explicit log logger mode is converted internally and not double logged', {
   File<-write_fake_tags_file(log(rep(seq(4, 64, length.out=1000))))
   Proc.data<-get.tags.data(File, tag.type="mk", light.scale="log")
   expect_equal(Proc.data$tagtype, "mk")
   expect_true(Proc.data$log_transformed)
   expect_equal(max(Proc.data$FLightR.data$Data$light), 64, tolerance=1e-8)
   expect_true(max(Proc.data$Twilight.log.light.mat.dusk) <= log(64) + 1e-8)
})

test_that('logger mode aliases select tag type and light scale', {
   File<-write_fake_tags_file(rep(seq(4, 64, length.out=1000)))
   Proc.data<-get.tags.data(File, logger.mode="mk-real")
   expect_equal(Proc.data$tagtype, "mk")
   expect_false(Proc.data$log_transformed)

   File.log<-write_fake_tags_file(log(rep(seq(4, 64, length.out=1000))))
   Proc.data.log<-get.tags.data(File.log, logger.mode="mk-log")
   expect_equal(Proc.data.log$tagtype, "mk")
   expect_true(Proc.data.log$log_transformed)
})

test_that('failed automatic import reports clear error by default', {
   File<-write_fake_tags_file(rep(seq(4, 500, length.out=1000)))
   expect_error(
      get.tags.data(File),
      "Automatic tag/light-scale detection was failed.*Observed light range"
   )
})

test_that('ambiguous automatic import stops by default', {
   File<-write_fake_tags_file(rep(seq(4, 11.32, length.out=1000)))
   expect_error(
      get.tags.data(File),
      "Automatic tag/light-scale detection was ambiguous.*Candidate modes"
   )
})

test_that('explicit generic TAGS mode requires logger-specific defaults from the user', {
   File<-write_fake_tags_file(rep(seq(4, 500, length.out=1000)))
   expect_error(
      get.tags.data(File, tag.type="tags", light.scale="raw"),
      "no default log.light.borders"
   )
   Proc.data<-get.tags.data(
      File,
      tag.type="tags",
      light.scale="raw",
      log.light.borders=log(c(4, 500)),
      log.irrad.borders=c(-6, 4),
      saves="max"
   )
   expect_equal(Proc.data$tagtype, "tags")
   expect_false(Proc.data$log_transformed)
})

test_that('custom tracker raw bounds are stored as metadata only', {
   File<-write_fake_tags_file(rep(seq(4, 500, length.out=1000)))
   Proc.data<-get.tags.data(
      File,
      tag.type="tags",
      light.scale="raw",
      raw.light.bounds=c(0, 65535),
      log.light.borders=log(c(5, 500)),
      log.irrad.borders=c(-6, 4),
      saves="max"
   )
   expect_equal(Proc.data$tagtype, "tags")
   expect_equal(Proc.data$log.light.borders, log(c(5, 500)))
   expect_equal(Proc.data$tag.settings$raw.light.bounds, c(0, 65535))
   expect_equal(Proc.data$tag.settings$log.light.borders, log(c(5, 500)))
   expect_false(Proc.data$tag.settings$inferred$raw.light.bounds)
   expect_false(Proc.data$tag.settings$inferred$log.light.borders)
   expect_true(max(Proc.data$Twilight.log.light.mat.dusk) <= log(500) + 1e-8)
   expect_equal(Proc.data$Metadata$Import$requested.tag.type, "tags")
   expect_equal(Proc.data$Metadata$Import$requested.light.scale, "raw")
   expect_equal(Proc.data$Metadata$Import$resolved.tag.type, "tags")
   expect_equal(Proc.data$Metadata$Import$resolved.light.scale, "raw")
   expect_equal(Proc.data$Metadata$Import$observed.max.light, 500)
   expect_equal(Proc.data$Metadata$Import$detection.status, "confident")
   expect_false(Proc.data$Metadata$Import$inferred$raw.light.bounds)
})

test_that('removed future import arguments are not part of the public API', {
   formals.now<-names(formals(get.tags.data))
   expect_false("zero.offset" %in% formals.now)
   expect_false("saturation.values" %in% formals.now)
   expect_false("auto.fail" %in% formals.now)
   expect_false("validate.light.range" %in% formals.now)
   expect_false("log.light.bounds" %in% formals.now)
})

test_that('structured tag detection reports confident, failed, and ambiguous cases', {
   mk.raw<-data.frame(light=c(0, 64))
   detected.mk<-detect.tag.type(mk.raw)
   expect_equal(detected.mk$status, "confident")
   expect_equal(detected.mk$tag.type, "mk")
   expect_equal(detected.mk$light.scale, "raw")
   expect_equal(detected.mk$candidate.modes, "mk-raw")

   unknown<-detect.tag.type(data.frame(light=c(4, 500)))
   expect_equal(unknown$status, "failed")
   expect_equal(unknown$candidate.modes, character(0))
   expect_match(format_detection_failure(unknown), "Observed light range")

   ambiguous<-detect.tag.type(data.frame(light=c(4, 11.32)))
   expect_equal(ambiguous$status, "ambiguous")
   expect_true(length(ambiguous$candidate.modes) > 1)
   expect_match(format_detection_failure(ambiguous), "Candidate modes")
})

test_that('tag detection can be resolved with partial explicit settings', {
   ambiguous<-detect.tag.type(data.frame(light=c(4, 11.32)))

   filtered.by.type<-filter.detection.candidates(ambiguous, tagtype="GDL1")
   expect_equal(filtered.by.type$status, "confident")
   expect_equal(filtered.by.type$tag.type, "GDL1")

   filtered.failed<-filter.detection.candidates(ambiguous, tagtype="mk")
   expect_equal(filtered.failed$status, "failed")

   expect_message(
      resolved.type<-resolve.tag.type(
         data.frame(light=c(4, 11.32)),
         tag.type="gdl1",
         light.scale="auto"
      ),
      "user-specified GDL1"
   )
   expect_equal(resolved.type$tagtype, "GDL1")
   expect_equal(resolved.type$light.scale, "log")

   expect_message(
      resolved.scale<-resolve.tag.type(
         data.frame(light=c(4, 64)),
         tag.type="auto",
         light.scale="raw"
      ),
      "user-specified light.scale=raw"
   )
   expect_equal(resolved.scale$tagtype, "mk")
   expect_false(resolved.scale$log_transformed)
})

test_that('explicit tag settings bypass automatic inference and warn on suspicious ranges', {
   expect_warning(
      explicit.log<-resolve.tag.type(data.frame(light=c(0, 64)),
         tag.type="mk", light.scale="log"),
      "observed max light"
   )
   expect_equal(explicit.log$tagtype, "mk")
   expect_true(explicit.log$log_transformed)
   expect_equal(explicit.log$detection$status, "confident")

   expect_warning(
      explicit.raw<-resolve.tag.type(data.frame(light=c(0, log(10))),
         tag.type="mk", light.scale="raw"),
      "already contain log-light"
   )
   expect_equal(explicit.raw$light.scale, "raw")

   expect_error(
      resolve.tag.type(data.frame(light=c(0, 64)), tag.type="unknown", light.scale="raw"),
      "Unknown tag.type"
   )
   expect_error(
      resolve.tag.type(data.frame(light=c(0, 64)), logger.mode="mk"),
      "logger.mode should combine"
   )
})

test_that('tag defaults cover known logger families', {
   expect_equal(get.tag.type.defaults("mk")$saves, "max")
   expect_equal(get.tag.type.defaults("GDL1")$saves, "mean")
   expect_equal(get.tag.type.defaults("GDL2v1")$log.irrad.borders, c(-2.7, 0))
   expect_equal(get.tag.type.defaults("Lat_2000")$log.light.borders, c(100, 360))
   expect_error(get.tag.type.defaults("custom"), "No default import settings")
})

test_that('legacy get.tag.type wrapper keeps warning-on-failure behaviour', {
   expect_message(tag<-get.tag.type(data.frame(light=c(0, 64))), "Detectedmk")
   expect_equal(tag$tagtype, "mk")

   expect_warning(
      failed<-get.tag.type(data.frame(light=c(4, 500))),
      "Automatic tag/light-scale detection was failed"
   )
   expect_null(failed)
})

test_that('import metadata records inferred and user supplied settings', {
   detection<-resolve.tag.type(data.frame(light=c(0, 64)),
      tag.type="mk", light.scale="raw")
   metadata<-make.import.metadata(
      requested.tag.type="mk",
      requested.light.scale="raw",
      requested.logger.mode=NULL,
      detected=detection,
      raw.light.bounds=c(0, 65535),
      log.light.borders=log(c(4, 64)),
      observed.light.range=c(0, 64),
      user.raw.light.bounds=TRUE
   )

   expect_equal(metadata$resolved.tag.type, "mk")
   expect_equal(metadata$resolved.light.scale, "raw")
   expect_equal(metadata$observed.max.light, 64)
   expect_false(metadata$inferred$tag.type)
   expect_false(metadata$inferred$raw.light.bounds)
   expect_equal(metadata$transformation.applied, "process.twilights uses log(light)")
})

