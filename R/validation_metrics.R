# Local validation helpers for private double-tagged data.

validation_yaml_value <- function(x) {
  x <- trimws(sub("#.*$", "", x))
  if (!nzchar(x) || identical(tolower(x), "null")) {
    return(NULL)
  }
  x <- sub("^['\"]", "", x)
  sub("['\"]$", "", x)
}

validation_read_simple_yaml <- function(path) {
  lines <- readLines(path, warn = FALSE)
  out <- list()
  section <- NULL

  for (line in lines) {
    if (!nzchar(trimws(line)) || grepl("^\\s*#", line)) {
      next
    }

    if (grepl("^[^[:space:]][^:]*:\\s*$", line)) {
      section <- trimws(sub(":\\s*$", "", line))
      out[[section]] <- list()
      next
    }

    key_value <- regmatches(line, regexec("^\\s*([^:]+):\\s*(.*)$", line))[[1]]
    if (length(key_value) != 3) {
      next
    }

    key <- trimws(key_value[2])
    value <- validation_yaml_value(key_value[3])
    if (is.null(section)) {
      out[[key]] <- value
    } else {
      out[[section]][[key]] <- value
    }
  }

  out
}

#' Read local validation configuration
#'
#' Reads a YAML configuration file for local, private double-tag validation.
#' The parser uses the yaml package when available and otherwise supports the
#' simple key/value structure used by the example config.
#'
#' @param path Path to a local validation config file.
#' @return A named list with validation configuration.
read_validation_config <- function(path) {
  if (!file.exists(path)) {
    stop("Validation config file does not exist: ", path, call. = FALSE)
  }

  if (requireNamespace("yaml", quietly = TRUE)) {
    return(yaml::read_yaml(path))
  }

  validation_read_simple_yaml(path)
}

#' Read a GPS track for validation
#'
#' Reads a CSV file containing GPS fixes. Required columns are time, lon, and
#' lat. The time column is converted to POSIXct in UTC when needed.
#'
#' @param path Path to a private local GPS CSV file.
#' @return A data.frame with POSIXct time and numeric lon/lat columns.
read_gps_track <- function(path) {
  if (!file.exists(path)) {
    stop("GPS track file does not exist: ", path, call. = FALSE)
  }

  gps <- utils::read.csv(path, stringsAsFactors = FALSE)
  required <- c("time", "lon", "lat")
  missing <- setdiff(required, names(gps))
  if (length(missing) > 0) {
    stop("GPS track is missing required columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  if (!inherits(gps$time, "POSIXct")) {
    gps$time <- as.POSIXct(gps$time, tz = "UTC")
  }
  gps$lon <- as.numeric(gps$lon)
  gps$lat <- as.numeric(gps$lat)

  gps
}

validation_twilight_times <- function(Result) {
  candidates <- list(
    Result$twilight,
    Result$twilights,
    Result$Twilight,
    Result$Data$twilight,
    Result$Data$twilights,
    Result$Results$twilight,
    Result$Results$twilights
  )

  for (candidate in candidates) {
    if (is.null(candidate)) {
      next
    }
    if (inherits(candidate, "POSIXt")) {
      return(as.POSIXct(candidate, tz = "UTC"))
    }
    if (is.data.frame(candidate)) {
      for (name in c("time", "Time", "twilight", "Twilight", "datetime", "DateTime")) {
        if (name %in% names(candidate)) {
          return(as.POSIXct(candidate[[name]], tz = "UTC"))
        }
      }
    }
  }

  stop("Could not find twilight times in Result.", call. = FALSE)
}

#' Match GPS fixes to FLightR twilight times
#'
#' @param Result A FLightR Result-like object containing twilight times.
#' @param gps A data.frame with time, lon, and lat columns.
#' @param max_time_diff_hours Maximum allowed absolute time difference.
#' @return A data.frame of matched twilight and GPS records.
match_gps_to_twilights <- function(Result, gps, max_time_diff_hours = 6) {
  twilight_times <- validation_twilight_times(Result)
  if (!inherits(gps$time, "POSIXct")) {
    gps$time <- as.POSIXct(gps$time, tz = "UTC")
  }

  matched <- lapply(seq_along(twilight_times), function(i) {
    diffs <- abs(as.numeric(difftime(gps$time, twilight_times[i], units = "hours")))
    best <- which.min(diffs)
    if (!length(best) || is.na(diffs[best]) || diffs[best] > max_time_diff_hours) {
      return(NULL)
    }
    data.frame(
      twilight_index = i,
      twilight_time = twilight_times[i],
      gps_time = gps$time[best],
      time_diff_hours = diffs[best],
      gps_lon = gps$lon[best],
      gps_lat = gps$lat[best]
    )
  })

  matched <- do.call(rbind, matched)
  if (is.null(matched)) {
    return(data.frame(
      twilight_index = integer(),
      twilight_time = as.POSIXct(character()),
      gps_time = as.POSIXct(character()),
      time_diff_hours = numeric(),
      gps_lon = numeric(),
      gps_lat = numeric()
    ))
  }

  rownames(matched) <- NULL
  matched
}

validation_points_from_matrix <- function(x, twilight_index) {
  if (length(dim(x)) == 3) {
    slice <- x[, , twilight_index]
    colnames(slice) <- colnames(x)
    x <- slice
  }
  x <- as.data.frame(x)

  lon_name <- intersect(c("lon", "Lon", "longitude", "Longitude", "x"), names(x))[1]
  lat_name <- intersect(c("lat", "Lat", "latitude", "Latitude", "y"), names(x))[1]
  weight_name <- intersect(c("weight", "weights", "prob", "probability", "p"), names(x))[1]

  if (is.na(lon_name) || is.na(lat_name)) {
    return(NULL)
  }

  points <- data.frame(
    lon = as.numeric(x[[lon_name]]),
    lat = as.numeric(x[[lat_name]])
  )
  points$weight <- if (is.na(weight_name)) 1 else as.numeric(x[[weight_name]])
  points
}

validation_result_candidates <- function(Result) {
  list(
    Result$posterior,
    Result$posterior_points,
    Result$particle_locations,
    Result$spatial,
    Result$Results$posterior,
    Result$Results$posterior_points,
    Result$Results$particle_locations
  )
}

#' Extract posterior point distribution for one twilight
#'
#' @param Result A FLightR Result-like object.
#' @param twilight_index Integer twilight index.
#' @return A data.frame with lon, lat, and weight columns.
posterior_point_distribution <- function(Result, twilight_index) {
  for (candidate in validation_result_candidates(Result)) {
    if (is.null(candidate)) {
      next
    }

    if (is.list(candidate) && !is.data.frame(candidate) && length(candidate) >= twilight_index) {
      points <- validation_points_from_matrix(candidate[[twilight_index]], twilight_index)
    } else {
      points <- validation_points_from_matrix(candidate, twilight_index)
    }

    if (!is.null(points)) {
      points <- points[stats::complete.cases(points[, c("lon", "lat", "weight")]), ]
      points <- points[points$weight > 0, ]
      if (nrow(points) > 0) {
        points$weight <- points$weight / sum(points$weight)
        rownames(points) <- NULL
        return(points)
      }
    }
  }

  stop("Could not find posterior point distribution for twilight index ",
       twilight_index, ".", call. = FALSE)
}

validation_distance_km <- function(lon1, lat1, lon2, lat2) {
  if (requireNamespace("geosphere", quietly = TRUE)) {
    return(geosphere::distHaversine(cbind(lon1, lat1), cbind(lon2, lat2)) / 1000)
  }

  radius_km <- 6371.0088
  lon1 <- lon1 * pi / 180
  lat1 <- lat1 * pi / 180
  lon2 <- lon2 * pi / 180
  lat2 <- lat2 * pi / 180
  dlon <- lon2 - lon1
  dlat <- lat2 - lat1
  a <- sin(dlat / 2)^2 + cos(lat1) * cos(lat2) * sin(dlon / 2)^2
  radius_km * 2 * atan2(sqrt(a), sqrt(1 - a))
}

validation_weighted_median <- function(x, w) {
  ord <- order(x)
  x <- x[ord]
  w <- w[ord] / sum(w)
  x[which(cumsum(w) >= 0.5)[1]]
}

validation_mode_point <- function(points) {
  points[which.max(points$weight), c("lon", "lat")]
}

validation_ci95_contains <- function(points, gps_lon, gps_lat) {
  point_dist <- validation_distance_km(points$lon, points$lat, gps_lon, gps_lat)
  ordered <- order(point_dist)
  cumulative <- cumsum(points$weight[ordered])
  cutoff <- point_dist[ordered][which(cumulative >= 0.95)[1]]
  mode_point <- validation_mode_point(points)
  gps_to_mode <- validation_distance_km(gps_lon, gps_lat, mode_point$lon, mode_point$lat)
  gps_to_mode <= cutoff
}

#' Summarize validation against matched GPS fixes
#'
#' @param Result A FLightR Result-like object with posterior point distributions.
#' @param gps_matched Output of match_gps_to_twilights().
#' @return A one-row data.frame with validation metrics.
validation_summary <- function(Result, gps_matched) {
  if (nrow(gps_matched) == 0) {
    return(data.frame(
      n_matched = 0L,
      median_error_km = NA_real_,
      mean_error_km = NA_real_,
      mode_error_km = NA_real_,
      ci95_coverage = NA_real_,
      notes = "No GPS fixes matched to twilight times."
    ))
  }

  errors <- lapply(seq_len(nrow(gps_matched)), function(i) {
    points <- posterior_point_distribution(Result, gps_matched$twilight_index[i])
    mode_point <- validation_mode_point(points)
    point_dist <- validation_distance_km(
      points$lon, points$lat, gps_matched$gps_lon[i], gps_matched$gps_lat[i])
    data.frame(
      median_error_km = validation_weighted_median(point_dist, points$weight),
      mean_error_km = stats::weighted.mean(point_dist, points$weight),
      mode_error_km = validation_distance_km(
        mode_point$lon, mode_point$lat, gps_matched$gps_lon[i], gps_matched$gps_lat[i]),
      ci95_contains = validation_ci95_contains(points, gps_matched$gps_lon[i], gps_matched$gps_lat[i])
    )
  })

  errors <- do.call(rbind, errors)
  notes <- character()
  if (nrow(errors) < 20) {
    notes <- c(notes, "CI95 coverage is approximate because fewer than 20 matches are available.")
  }
  if (any(!stats::complete.cases(errors))) {
    notes <- c(notes, "Some posterior metrics could not be calculated.")
  }

  data.frame(
    n_matched = nrow(gps_matched),
    median_error_km = stats::median(errors$median_error_km, na.rm = TRUE),
    mean_error_km = mean(errors$mean_error_km, na.rm = TRUE),
    mode_error_km = stats::median(errors$mode_error_km, na.rm = TRUE),
    ci95_coverage = mean(errors$ci95_contains, na.rm = TRUE),
    notes = if (length(notes)) paste(unique(notes), collapse = " ") else ""
  )
}
