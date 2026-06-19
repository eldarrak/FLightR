#' 
#' read TAGS formatted data
#' 
#' Reads the data frame with detected twilight events into the FLightR
#' 
#' @param filename the name of the file which the data are to be read from. File is supposed to be comma separated file of TAGS format. If it does not contain an absolute path, the file name is relative to the current working directory, getwd(). Tilde-expansion is performed where supported. This can be a compressed file. Alternatively, file can be a readable text-mode connection (which will be opened for reading if necessary, and if so closed (and hence destroyed) at the end of the function call). File can also be a complete URL.
#' @param start.date date of beginning of relevant data collection in \code{POSIXct} format.
#' @param end.date date of end of relevant data collection in \code{POSIXct} format.
#' @param log.light.borders Numeric vector with length of 2 for minimum and maximum log(light) levels to use. Alternatively character value 'auto', that will allow FLightR to assign these values according to detected tag type.
#' @param log.irrad.borders Numeric vector with length of 2 for minimum and maximum log(irradiance) values to use. Alternatively character value 'auto', that will allow FLightR to assign these values according to detected tag type.
#' @param saves character values informing FLightR if min or max values were used by logger.
#' @param measurement.period Value in seconds defining how often tag was measuring light levels. If NULL value will be taken from known values for detected tag type.
#' @param impute.on.boundaries logical, if FLightR should approximate values at boundaries. Set it to TRUE only if you have vary few active points at each twilight, e.g if tag was saving every 10 minutes or so.
#' @param tag.type Logger type. The default, \code{"auto"}, keeps the historical automatic detection from light minimum and maximum values. Explicit values such as \code{"mk"}, \code{"intigeo"}, \code{"intigeo_mode_1"}, \code{"gdl1"}, \code{"gdl2v1"}, \code{"gdl2v2_gdlpam"}, and \code{"lat2000"} bypass automatic tag-type detection.
#' @param light.scale Scale of the light values in the TAGS file. Use \code{"auto"} for historical automatic detection, \code{"raw"} for untransformed positive light values, or \code{"log"} when the file already contains log-light values. If \code{"log"} is supplied, FLightR converts values back to raw light internally so downstream processing does not log-transform them twice.
#' @param logger.mode Optional convenience alias combining \code{tag.type} and \code{light.scale}, e.g. \code{"mk-real"}, \code{"mk-raw"}, \code{"mk-log"}, \code{"intigeo-real"}, or \code{"intigeo-log"}. If supplied, it overrides \code{tag.type} and \code{light.scale}.
#' @param raw.light.bounds Optional raw light-value bounds for a custom logger. These are stored in \code{Proc.data$tag.settings} and \code{Proc.data$Metadata$Import} for diagnostics.
#' @param tag.settings Optional list with any of \code{tag.type}, \code{light.scale}, \code{logger.mode}, or \code{raw.light.bounds}. Values in this list override the corresponding arguments.
#' @return list, which is to be further processed with the FLightR.
#' @details The returned object has many parts, the important are: (1) the recorded light data, (2) the detected twilight events, (3) light level data at the moment of each determined sunrise and sunset and around them (24 fixes before and 24 after), and (4) technical parameters of the tag, i. e. its type, saving and measuring period (the periodicity, in seconds, at which a tag measures and saves data). Import settings are stored both in the historical top-level fields, such as \code{tagtype}, \code{log_transformed}, \code{log.light.borders}, and \code{log.irrad.borders}, and in \code{Proc.data$Metadata$Import} for structured provenance.
#' @examples
#' File<-system.file("extdata", "Godwit_TAGS_format.csv", package = "FLightR")
#' Proc.data<-get.tags.data(File)
#' Proc.data<-get.tags.data(File, tag.type = "intigeo", light.scale = "log")
#' @export
get.tags.data<-function(filename=NULL, start.date=NULL, end.date=NULL, log.light.borders='auto', log.irrad.borders='auto', saves=c("auto", "max", "mean"), measurement.period=NULL,  impute.on.boundaries=FALSE, tag.type="auto", light.scale="auto", logger.mode=NULL, raw.light.bounds=NULL, tag.settings=NULL) {
 # measurement.period should be set to saving period for swiss tags
   if (!is.null(tag.settings)) {
      if (!is.null(tag.settings$tag.type)) tag.type<-tag.settings$tag.type
      if (!is.null(tag.settings$light.scale)) light.scale<-tag.settings$light.scale
      if (!is.null(tag.settings$logger.mode)) logger.mode<-tag.settings$logger.mode
      if (!is.null(tag.settings$raw.light.bounds)) raw.light.bounds<-tag.settings$raw.light.bounds
   }
   user.raw.light.bounds<-!is.null(raw.light.bounds)
   if (!is.null(raw.light.bounds) & (length(raw.light.bounds)!=2 | !is.numeric(raw.light.bounds))) stop("raw.light.bounds should be NULL or a numeric vector of length 2\n")
   requested.tag.type<-tag.type
   requested.light.scale<-light.scale
   requested.logger.mode<-logger.mode

   TAGS.twilights<-utils::read.csv(filename, stringsAsFactors =FALSE)
   observed.light.range<-range(TAGS.twilights$light, na.rm=TRUE)

   # now we have to find tag type and figure out whether data were logtransformed..

   detected<-resolve.tag.type(TAGS.twilights, tag.type=tag.type, light.scale=light.scale, logger.mode=logger.mode)

   if (is.null(detected)) {
      if (log.light.borders[1]=='auto') stop("Unrecognized tag type, please supply log.light.borders or set tag.type/light.scale explicitly\n")
      if (log.irrad.borders[1]=='auto') stop("Unrecognized tag type, please supply log.irrad.borders or set tag.type/light.scale explicitly\n")
      if (length(saves)==3 | saves[1]=='auto') stop("Unrecognized tag type, please tell FLightR what the tag was saving - max or mean over saving period\n")
   } else {
      if (detected$log_transformed) TAGS.twilights$light<-exp(TAGS.twilights$light)
      if (detected$tagtype=="tags") {
         if (log.light.borders[1]=='auto') stop("Generic tag.type='tags' has no default log.light.borders; please supply log.light.borders\n")
         if (log.irrad.borders[1]=='auto') stop("Generic tag.type='tags' has no default log.irrad.borders; please supply log.irrad.borders\n")
         if (length(saves)==3 | saves[1]=='auto') stop("Generic tag.type='tags' has no default saves setting; please supply saves='max' or saves='mean'\n")
      } else {
         defaults<-get.tag.type.defaults(detected$tagtype)
         if (log.light.borders[1]=='auto') log.light.borders<-defaults$log.light.borders
         if (log.irrad.borders[1]=='auto') log.irrad.borders<-defaults$log.irrad.borders
         if (saves[1] =='auto') saves<-defaults$saves
      }
      if(detected$tagtype=="Lat_2000" & is.null(measurement.period)) measurement.period<-10
   }

   FLightR.data<-read.tags.light.twilight(TAGS.twilights,
                      start.date=start.date, end.date=end.date)

   # ok, now we want to specify measurement period
   saving.period<-round(diff(as.numeric(FLightR.data$Data$gmt[1:2]))/10)*10
   if (is.null(measurement.period)) {
      if(saves=="mean") {
	     measurement.period<-saving.period	
      } else {
         measurement.period<-60
      }
   }
   
   message("tag saved data every", saving.period, "seconds, and is assumed to measure data every", measurement.period, "seconds, and write down", saves[1], "\n" )
   if (max(TAGS.twilights$light) ==64 & saving.period>500) {
      impute.on.boundaries=TRUE
      warning("saving period was too long for this type of tag, FLightR will impute data\n")
   }
   Proc.data<-process.twilights(FLightR.data$Data,FLightR.data$twilights, 
                             measurement.period=measurement.period, saving.period=saving.period, impute.on.boundaries=impute.on.boundaries)
   if (!is.null(detected)) {
   Proc.data$tagtype<-detected$tagtype
   Proc.data$log_transformed<-detected$log_transformed
   }
   Proc.data$log.light.borders=log.light.borders
   Proc.data$log.irrad.borders=log.irrad.borders
   import.metadata<-make.import.metadata(
      requested.tag.type=requested.tag.type,
      requested.light.scale=requested.light.scale,
      requested.logger.mode=requested.logger.mode,
      detected=detected,
      raw.light.bounds=raw.light.bounds,
      log.light.borders=log.light.borders,
      observed.light.range=observed.light.range,
      user.raw.light.bounds=user.raw.light.bounds
   )
   Proc.data$tag.settings<-list(
      requested.tag.type=import.metadata$requested.tag.type,
      requested.light.scale=import.metadata$requested.light.scale,
      requested.logger.mode=import.metadata$requested.logger.mode,
      tag.type=import.metadata$resolved.tag.type,
      light.scale=import.metadata$resolved.light.scale,
      raw.light.bounds=import.metadata$raw.light.bounds,
      log.light.borders=import.metadata$log.light.borders,
      observed.light.range=import.metadata$observed.light.range,
      detection.status=import.metadata$detection.status,
      detection.reason=import.metadata$detection.message,
      detection.message=import.metadata$detection.message,
      detection.candidate.modes=import.metadata$detection.candidate.modes,
      detection=if (!is.null(detected$detection)) detected$detection else NULL,
      transformation.applied=import.metadata$transformation.applied,
      inferred=import.metadata$inferred
   )
   Proc.data$Metadata<-list(Import=import.metadata)
   attr(Proc.data, "tag.settings")<-Proc.data$tag.settings
   Proc.data$FLightR.data<-FLightR.data
   return(Proc.data)			
}


make.import.metadata<-function(requested.tag.type=NULL, requested.light.scale=NULL, requested.logger.mode=NULL, detected=NULL, raw.light.bounds=NULL, log.light.borders=NULL, observed.light.range=c(NA_real_, NA_real_), user.raw.light.bounds=FALSE) {
   resolved.tag.type<-if (!is.null(detected)) detected$tagtype else requested.tag.type
   resolved.light.scale<-if (!is.null(detected) && isTRUE(detected$log_transformed)) "log" else if (!is.null(detected)) "raw" else requested.light.scale
   detection<-if (!is.null(detected$detection)) detected$detection else NULL
   list(
      requested.tag.type=requested.tag.type,
      requested.light.scale=requested.light.scale,
      requested.logger.mode=requested.logger.mode,
      resolved.tag.type=resolved.tag.type,
      resolved.light.scale=resolved.light.scale,
      raw.light.bounds=raw.light.bounds,
      log.light.borders=log.light.borders,
      observed.light.range=observed.light.range,
      observed.min.light=observed.light.range[1],
      observed.max.light=observed.light.range[2],
      transformation.applied=if (!is.null(detected) && isTRUE(detected$log_transformed)) "exp(light) before processing; process.twilights uses log(light)" else "process.twilights uses log(light)",
      detection.status=if (!is.null(detection)) detection$status else NA_character_,
      detection.message=if (!is.null(detection)) detection$message else NA_character_,
      detection.candidate.modes=if (!is.null(detection)) detection$candidate.modes else character(0),
      inferred=list(
         tag.type=is.null(requested.logger.mode) && identical(requested.tag.type, "auto"),
         light.scale=is.null(requested.logger.mode) && identical(requested.light.scale, "auto"),
         raw.light.bounds=!user.raw.light.bounds,
         log.light.borders=identical(log.light.borders, "auto")
      )
   )
}


resolve.tag.type<-function(TAGS.twilights, tag.type="auto", light.scale="auto", logger.mode=NULL) {

   explicit.mode<-!is.null(logger.mode)
   if (explicit.mode) {
      mode<-tolower(gsub("_", "-", logger.mode))
      parts<-strsplit(mode, "-", fixed=TRUE)[[1]]
      if (length(parts)<2) stop("logger.mode should combine tag type and light scale, for example 'mk-real', 'mk-log', 'intigeo-real', or 'intigeo-log'\n")
      scale.alias<-parts[length(parts)]
      type.alias<-paste(parts[-length(parts)], collapse="-")
      tag.type<-type.alias
      light.scale<-scale.alias
   }

   tag.type<-tolower(gsub("-", "_", tag.type))
   light.scale<-tolower(light.scale)
   if (light.scale=="real") light.scale<-"raw"
   if (!light.scale %in% c("auto", "raw", "log")) stop("light.scale should be one of 'auto', 'raw', or 'log'\n")

   tag.map<-c(
      auto="auto",
      mk="mk",
      intigeo="Intigeo_Mode_1",
      intigeo_mode_1="Intigeo_Mode_1",
      intigeo_mode_4="Intigeo_Mode_4",
      intigeo_mode_6="Intigeo_Mode_6_clipped",
      intigeo_mode_6_clipped="Intigeo_Mode_6_clipped",
      gdl2v2_gdlpam="GDL2v2_GDLpam",
      gdlpam="GDL2v2_GDLpam",
      gdl1="GDL1",
      gdl2v1="GDL2v1",
      lat2000="Lat_2000",
      lat_2000="Lat_2000",
      tags="tags"
   )
   if (!tag.type %in% names(tag.map)) stop("Unknown tag.type '", tag.type, "'. Use 'auto', 'mk', 'intigeo', 'intigeo_mode_1', 'intigeo_mode_4', 'gdl1', 'gdl2v1', 'gdl2v2_gdlpam', 'lat2000', or 'tags'.\n")

   detection<-detect.tag.type(TAGS.twilights)
   if (tag.type=="auto" & light.scale=="auto" & !explicit.mode) return(resolve.auto.detection(detection))
   if (tag.type=="auto" & light.scale!="auto") {
      detection<-filter.detection.candidates(detection, light.scale=light.scale)
      detected<-resolve.auto.detection(detection)
      detected$log_transformed<-light.scale=="log"
      detected$light.scale<-light.scale
      message("Using detected ", detected$tagtype, " tag with user-specified light.scale=", light.scale, "\n")
      return(detected)
   }
   if (tag.type!="auto" & light.scale=="auto") {
      detection<-filter.detection.candidates(detection, tagtype=unname(tag.map[tag.type]))
      detected<-resolve.auto.detection(detection)
      detected$tagtype<-unname(tag.map[tag.type])
      message("Using user-specified ", detected$tagtype, " tag with detected light.scale=", detected$light.scale, "\n")
      return(detected)
   }

   tagtype<-unname(tag.map[tag.type])
   Max_light<-max(TAGS.twilights$light, na.rm=TRUE)
   Min_light<-min(TAGS.twilights$light, na.rm=TRUE)
   consistency<-check.explicit.light.consistency(Max_light, Min_light, light.scale)
   if (!is.null(consistency)) warning(consistency, call.=FALSE)
   message("Using user-specified ", tagtype, " tag\n")
   message("Using user-specified light.scale=", light.scale, "; observed light range: ", signif(Min_light, 6), " to ", signif(Max_light, 6), "\n")
   if (light.scale=="log") message("Input light values are treated as log-transformed and converted internally before processing\n")
   detection<-make.detection.result(
      status="confident",
      reason="tag.type and light.scale were supplied by the user; automatic detection was not used",
      observed.min=Min_light,
      observed.max=Max_light,
      candidates=data.frame(tagtype=tagtype, light.scale=light.scale, log_transformed=light.scale=="log", score=Inf, reason="user-specified", stringsAsFactors=FALSE)
   )
   return(list(tagtype=tagtype, light.scale=light.scale, log_transformed=light.scale=="log", detection=detection))
}


resolve.auto.detection<-function(detection) {
   if (detection$status=="confident") {
      cand<-detection$candidates[1,]
      message("Detected", cand$tagtype, "tag\n")
      if (cand$log_transformed) message("Data found to be logtransformed\n")
      return(list(tagtype=cand$tagtype, light.scale=cand$light.scale, log_transformed=cand$log_transformed, detection=detection))
   }
   stop(format_detection_failure(detection), call.=FALSE)
}


filter.detection.candidates<-function(detection, tagtype=NULL, light.scale=NULL) {
   candidates<-detection$candidates
   if (!is.null(tagtype) && nrow(candidates)>0) candidates<-candidates[candidates$tagtype==tagtype,,drop=FALSE]
   if (!is.null(light.scale) && nrow(candidates)>0) candidates<-candidates[candidates$light.scale==light.scale,,drop=FALSE]
   if (nrow(candidates)==0) {
      return(make.detection.result(
         status="failed",
         reason="automatic detection found no candidates matching the user-specified partial settings",
         observed.min=detection$observed.min,
         observed.max=detection$observed.max,
         candidates=candidates
      ))
   }
   make.detection.result(
      status=if (length(unique(candidates$score))>0 && sum(candidates$score==max(candidates$score))==1) "confident" else "ambiguous",
      reason=if (sum(candidates$score==max(candidates$score))==1) "unique best candidate after applying user-specified partial settings" else "multiple candidate modes match the observed light range after applying user-specified partial settings",
      observed.min=detection$observed.min,
      observed.max=detection$observed.max,
      candidates=candidates
   )
}


check.explicit.light.consistency<-function(Max_light, Min_light, light.scale) {
   if (light.scale=="log" && Max_light>20) return(paste0("Explicit light.scale='log' was supplied, but observed max light is ", signif(Max_light, 6), "; check whether the file may contain raw light values."))
   if (light.scale=="raw" && Max_light<20 && Min_light>=0) return(paste0("Explicit light.scale='raw' was supplied, but observed light range is ", signif(Min_light, 6), " to ", signif(Max_light, 6), "; check whether the file may already contain log-light values."))
   NULL
}


get.tag.type.defaults<-function(tagtype) {
   if(tagtype=="Intigeo_Mode_1") return(list(log.light.borders=c(1.5, 9), log.irrad.borders=c(-3,3), saves="max"))
   if(tagtype=="Intigeo_Mode_4") return(list(log.light.borders=c(1.5, 7), log.irrad.borders=c(-3,3), saves="max"))
   if(tagtype=="Intigeo_Mode_6_clipped") return(list(log.light.borders=c(1.5, 7), log.irrad.borders=c(-3,3), saves="max"))
   if(tagtype=="mk") return(list(log.light.borders=log(c(4, 61)), log.irrad.borders=c(-6.5,4), saves="max"))
   if(tagtype=="GDL2v2_GDLpam") return(list(log.light.borders=c(2.5, 8), log.irrad.borders=c(-6.5,1.5), saves="mean"))
   if(tagtype=="GDL1") return(list(log.light.borders=c(3, 7), log.irrad.borders=c(-4,1), saves="mean"))
   if(tagtype=="GDL2v1") return(list(log.light.borders=c(2.2,7), log.irrad.borders=c(-2.7,0), saves="mean"))
   if(tagtype=="Lat_2000") return(list(log.light.borders=c(100,360), log.irrad.borders=c(-8,2), saves="max"))
   stop("No default import settings are known for tag type '", tagtype, "'. Please supply log.light.borders, log.irrad.borders, and saves explicitly.\n")
}


get.tag.type<-function(TAGS.twilights) {
   detection<-detect.tag.type(TAGS.twilights)
   if (detection$status!="confident") {
      warning(format_detection_failure(detection), call.=FALSE)
      return(NULL)
   }
   cand<-detection$candidates[1,]
   message("Detected", cand$tagtype, "tag\n")
   if (cand$log_transformed) message("Data found to be logtransformed\n")
   list(tagtype=cand$tagtype, light.scale=cand$light.scale, log_transformed=cand$log_transformed, detection=detection)
}


detect.tag.type<-function(TAGS.twilights) {
   Max_light<-max(TAGS.twilights$light, na.rm=TRUE)
   Min_light<-min(TAGS.twilights$light, na.rm=TRUE)
   candidates<-data.frame(
      tagtype=character(0),
      light.scale=character(0),
      log_transformed=logical(0),
      score=numeric(0),
      reason=character(0),
      stringsAsFactors=FALSE
   )
   add_candidate<-function(tagtype, light.scale, score, reason) {
      candidates<<-rbind(candidates, data.frame(
         tagtype=tagtype,
         light.scale=light.scale,
         log_transformed=light.scale=="log",
         score=score,
         reason=reason,
         stringsAsFactors=FALSE
      ))
   }

   if(Max_light == 1146.681 &  Min_light == 0.32) add_candidate("Intigeo_Mode_6_clipped", "raw", 1, "exact raw Intigeo Mode 6 clipped min/max")
   if(Max_light == log(1146.681) &  Min_light == log(0.32)) add_candidate("Intigeo_Mode_6_clipped", "log", 1, "exact log Intigeo Mode 6 clipped min/max")
   if(round(Max_light,2) >= 10.05 & round(Max_light,2) <= 12.50) add_candidate("Intigeo_Mode_1", "log", 0.8, "max light within historical log Intigeo Mode 1 range")
   if(round(Max_light,2) >= exp(10.05) & round(Max_light,2) <= exp(12.50)) add_candidate("Intigeo_Mode_1", "raw", 0.8, "max light within historical raw Intigeo Mode 1 range")
   if(round(Max_light,2)==7.06) add_candidate("Intigeo_Mode_4", "log", 1, "exact log Intigeo Mode 4 maximum")
   if(round(Max_light/10)==116) add_candidate("Intigeo_Mode_4", "raw", 0.9, "rounded raw Intigeo Mode 4 maximum")
   if(Max_light == 64) add_candidate("mk", "raw", 1, "exact mk raw maximum")
   if(Max_light == log(64)) add_candidate("mk", "log", 1, "exact mk log maximum")
   if(round(Max_light,1)==9.2) add_candidate("GDL2v2_GDLpam", "log", 1, "rounded GDL2v2/GDLpam log maximum")
   if(round(Max_light/10)==998) add_candidate("GDL2v2_GDLpam", "raw", 0.9, "rounded GDL2v2/GDLpam raw maximum")
   if(round(Max_light,2)==11.32) add_candidate("GDL1", "log", 1, "rounded GDL1 log maximum")
   if(round(Max_light)==82863) add_candidate("GDL1", "raw", 0.9, "rounded GDL1 raw maximum")
   if(round(Max_light,2)==log(63)) add_candidate("GDL2v1", "log", 1, "exact GDL2v1 log maximum")
   if(Max_light==63) add_candidate("GDL2v1", "raw", 1, "exact GDL2v1 raw maximum")
   if(log(Max_light) %in% c(4095, 357)) add_candidate("Lat_2000", "raw", 1, "historical Lat_2000 raw check")
   if(Max_light %in% c(4095, 357)) add_candidate("Lat_2000", "log", 1, "historical Lat_2000 log check")

   if (nrow(candidates)==0) {
      return(make.detection.result(
         status="failed",
         reason="no known tag/light-scale rule matched the observed light range",
         observed.min=Min_light,
         observed.max=Max_light,
         candidates=candidates
      ))
   }
   make.detection.result(
      status=if (nrow(candidates)==1) "confident" else "ambiguous",
      reason=if (nrow(candidates)==1) candidates$reason[1] else "multiple candidate modes match the observed light range",
      observed.min=Min_light,
      observed.max=Max_light,
      candidates=candidates
   )
}


make.detection.result<-function(status, reason, observed.min, observed.max, candidates) {
   candidate.modes<-if (nrow(candidates)==0) character(0) else paste(candidates$tagtype, candidates$light.scale, sep="-")
   list(
      tag.type=if (status=="confident" && nrow(candidates)>0) candidates$tagtype[1] else NA_character_,
      light.scale=if (status=="confident" && nrow(candidates)>0) candidates$light.scale[1] else NA_character_,
      status=status,
      reason=reason,
      message=reason,
      observed.min=observed.min,
      observed.max=observed.max,
      candidate.modes=candidate.modes,
      candidates=candidates
   )
}


format_detection_failure<-function(detection) {
   candidates<-if (length(detection$candidate.modes)==0) "none" else paste(detection$candidate.modes, collapse=", ")
   paste0(
      "Automatic tag/light-scale detection was ", detection$status, ": ", detection$reason, ".\n",
      "Observed light range: ", signif(detection$observed.min, 6), " to ", signif(detection$observed.max, 6), ".\n",
      "Candidate modes: ", candidates, ".\n",
      "Please specify settings explicitly, for example:\n",
      "  get.tags.data(file, tag.type='mk', light.scale='raw')\n",
      "  get.tags.data(file, tag.type='mk', light.scale='log')\n",
      "  get.tags.data(file, tag.type='tags', light.scale='raw', log.light.borders=log(c(5, 500)), log.irrad.borders=c(-6, 4), saves='max')\n"
   )
}


convert.lux.to.tags<-function(file, log=FALSE, log.light.borders=c(1,10)) {
	# the function takes the current (2015)
	# .lux format and converts it to .csv format that
	# TAGS service accepts
	warning("\n\nwe have figured out that TAGS service is rounding data, so log=TRUE will produce wrong rounding on upload\n\nsetting log to FALSE\n")
	log=FALSE
	Dat<-utils::read.csv(file, skip=20, sep="\t", stringsAsFactors =FALSE)
	names(Dat)<-c("datetime", "light")
	Dat$datetime<-as.POSIXct(Dat$datetime, tz="GMT", format="%d/%m/%Y %H:%M:%S")

	Dat_new<-Dat
	Dat_new$datetime<-format(Dat_new$datetime, format="%Y-%m-%d %H:%M:%S")
	
	if (log) {
	Dat_new$light<-log(Dat_new$light)
	
	if (!is.null(log.light.borders)) {
	Dat_new$light[Dat_new$light<log.light.borders[1]]<-log.light.borders[1]
	Dat_new$light[Dat_new$light>log.light.borders[2]]<-log.light.borders[2]
	}
	}
	# now I need to save as csv...
	utils::write.csv(Dat_new, file=paste(unlist(strsplit(file, ".lux")), "csv", sep="."), quote =FALSE, row.names=FALSE) 
	message("Success!\n")
	message("file", paste(unlist(strsplit(file, ".lux")), "csv", sep="."), "\nwas saved to", getwd(), "\n")
	return(NULL)
}


read.tags.light.twilight<-function(lig.raw, start.date=NULL, end.date=NULL) {
	lig.raw$datetime<-as.POSIXct(lig.raw$datetime, tz="GMT", format="%Y-%m-%dT%T")

	###
	### I also want to exclude the last days...
	
	if (!is.null(start.date)) lig.raw<-lig.raw[as.numeric(lig.raw$datetime) > as.numeric(as.POSIXct(start.date , tz="GMT")),]
	if (!is.null(end.date)) lig.raw<-lig.raw[as.numeric(lig.raw$datetime) < as.numeric(as.POSIXct(end.date , tz="GMT")),]

	## and also I want to exclude the interpolated and excluded afterwards points
		
	lig<-lig.raw[-which(lig.raw$interp==TRUE),]

##########################################################################
##Convert the datetime/light fields into the format that FLightR works on##
##########################################################################

	lig.new<-data.frame(
	format(lig$datetime, "%d/%m/%Y"),
	format(lig$datetime, "%H:%M:%S"),
 	lig$light)

########################################################
##Save that date/light file so it can be read in later##
########################################################
	tmpname<-tempfile(fileext = ".csv")
	utils::write.table(lig.new, file=tmpname, sep="," , row.names = FALSE, col.names =FALSE,  quote=FALSE)
	Data<-geologger.read.data(file=tmpname)
	unlink(x=tmpname)
#########################################################
##Use TAGS output to define twilight periods################
#########################################################
# 

Filtered_tw <- lig.raw[(lig.raw$twilight)>0 & !lig.raw$excluded, c("datetime", "twilight", "light")]
names(Filtered_tw)[2]="type"
Filtered_tw <- Filtered_tw[!duplicated(Filtered_tw$datetime),]
Filtered_tw <- Filtered_tw[order(Filtered_tw[,1]),]

Filtered_tw$id<-0
Data$d$type<-0
Data$d<-rbind(Data$d, data.frame(id=Filtered_tw$id, gmt= Filtered_tw$datetime, light=Filtered_tw$light, type=Filtered_tw$type))

Filtered_tw$excluded=0

All.p<-Data$d[order(Data$d$gmt),]

#All.p<-All.p[!duplicated(All.p[,2:3], fromLast=TRUE),]
rownames(All.p)<-1:nrow(All.p)
Res<-list(Data=All.p, twilights=Filtered_tw)
return(Res)
}


process.geolight.output<-function(datetime, light, gl_twl) {

#filter <- loessFilter(gl_twl[,1],gl_twl[,2],gl_twl[,3],k=10)
Filtered_tw <- data.frame(datetime=as.POSIXct(c(gl_twl$tFirst,gl_twl$tSecond),"GMT"),type=c(gl_twl$type,ifelse(gl_twl$type==1,2,1)))

Filtered_tw <- Filtered_tw[!duplicated(Filtered_tw$datetime),]
Filtered_tw <- Filtered_tw[order(Filtered_tw[,1]),]

# now I want to pair data and twilights..		  
Filtered_tw$light<-stats::approx(x=datetime, y=light, xout=Filtered_tw$datetime)$y
Filtered_tw$id<-0
Data<-data.frame(id=1:length(datetime), gmt=datetime, light=light, type=0)

Data<-rbind(Data, data.frame(id=Filtered_tw$id, gmt= Filtered_tw$datetime, light=Filtered_tw$light, type=Filtered_tw$type))

All.p<-Data[order(Data$gmt),]
#All.p<-All.p[!duplicated(All.p[,2:3], fromLast=TRUE),]
rownames(All.p)<-1:nrow(All.p)
Filtered_tw$excluded=0
FLightR.data=list(Data=All.p, twilights=Filtered_tw)
return(FLightR.data)

}


geologger.read.data<-function( file) {
	track <- utils::read.csv(file,header=FALSE, stringsAsFactors=FALSE) 
	names(track)<-c('date','time','light')
	track$datetime <- paste(track$date,track$time,sep=' ') #makes a new column called datetime with date and time concatenated together with a space between
	track$gmt <- as.POSIXct(strptime(track$datetime,'%d/%m/%Y %H:%M:%S'),tz='GMT') #makes a new column called gmt with date and time data in a formate that R can work with
	d <- data.frame(id =1:nrow(track), gmt = track$gmt, light = track$light) 
	Data<-list(d=d)
	return(Data)
}
