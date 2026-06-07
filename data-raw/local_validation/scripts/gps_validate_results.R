output_dir <- "D:/GitHub/FLightR/data-raw/local_validation/outputs"
gps_file <- "D:/GitHub/FLightR/data-raw/local_validation/private_data/A1_GPS_positions.csv"
max_time_diff_hours <- 12

if (!dir.exists(output_dir)) {
  stop("Output directory not found: ", output_dir, call. = FALSE)
}
if (!file.exists(gps_file)) {
  stop("GPS file not found: ", basename(gps_file), call. = FALSE)
}
if (!requireNamespace("geosphere", quietly = TRUE)) {
  stop("Package geosphere is required.", call. = FALSE)
}

guess_col <- function(names, candidates) {
  lower <- tolower(names)
  for (candidate in candidates) {
    idx <- which(lower == tolower(candidate))
    if (length(idx) > 0) {
      return(names[idx[1]])
    }
  }
  idx <- grep(paste(candidates, collapse = "|"), lower)
  if (length(idx) > 0) {
    return(names[idx[1]])
  }
  NA_character_
}

gps_raw <- utils::read.csv(gps_file, stringsAsFactors = FALSE)
time_col <- guess_col(names(gps_raw), c("time", "datetime", "date", "timestamp", "gmt"))
lon_col <- guess_col(names(gps_raw), c("lon", "longitude", "x"))
lat_col <- guess_col(names(gps_raw), c("lat", "latitude", "y"))
if (any(is.na(c(time_col, lon_col, lat_col)))) {
  stop("Could not identify GPS time/lon/lat columns safely.", call. = FALSE)
}

gps <- data.frame(
  time = as.POSIXct(gps_raw[[time_col]], tz = "GMT"),
  lon = as.numeric(gps_raw[[lon_col]]),
  lat = as.numeric(gps_raw[[lat_col]])
)
gps <- gps[stats::complete.cases(gps), ]
gps <- gps[order(gps$time), ]
if (nrow(gps) < 2) {
  stop("GPS data contain fewer than two usable fixes.", call. = FALSE)
}

benchmark_combined_path <- file.path(output_dir, "benchmark_combined.csv")
if (!file.exists(benchmark_combined_path)) {
  source("D:/GitHub/FLightR/data-raw/local_validation/scripts/summarize_existing_benchmarks.R")
}
bench <- utils::read.csv(benchmark_combined_path, stringsAsFactors = FALSE)

infer_version <- function(git_commit, git_branch) {
  if (!is.na(git_commit) && git_commit == "86fcb29") return("old_master")
  if ((!is.na(git_commit) && git_commit == "fe959c1") ||
      (!is.na(git_branch) && git_branch == "local-validation-framework")) return("new_branch")
  "unknown"
}

result_files <- list.files(output_dir, pattern = "^Result_.*\\.rds$", full.names = TRUE)
if (length(result_files) == 0) {
  stop("No Result_*.rds files found in outputs.", call. = FALSE)
}

metadata_for_result <- function(path) {
  base <- basename(path)
  row <- bench[basename(bench$result_path) == base, , drop = FALSE]
  if (nrow(row) == 0) {
    row <- bench[grepl(gsub("\\.rds$", "", base), bench$result_path, fixed = TRUE), , drop = FALSE]
  }
  if (nrow(row) > 0) {
    row <- row[1, , drop = FALSE]
    return(row)
  }
  data.frame(
    run_label = sub("^Result_|\\.rds$", "", base),
    git_branch = NA_character_,
    git_commit = NA_character_,
    nParticles = NA_integer_,
    seed = NA_integer_,
    subset = NA_character_,
    nTwilights = NA_integer_,
    LL = NA_real_,
    particle_filter_seconds = NA_real_,
    total_seconds = NA_real_,
    result_path = path,
    version = "unknown",
    stringsAsFactors = FALSE
  )
}

get_col <- function(x, candidates) {
  lower <- tolower(names(x))
  for (candidate in candidates) {
    idx <- which(lower == tolower(candidate))
    if (length(idx) > 0) return(names(x)[idx[1]])
  }
  NA_character_
}

result_times <- function(Result, n) {
  q <- Result$Results$Quantiles
  time_col <- get_col(q, c("time", "Time", "datetime", "DateTime"))
  if (!is.na(time_col)) {
    return(as.POSIXct(q[[time_col]], tz = "GMT"))
  }
  mt <- Result$Indices$Matrix.Index.Table
  for (candidate in c("Real.time", "time")) {
    if (!is.null(mt[[candidate]])) {
      return(as.POSIXct(mt[[candidate]][seq_len(min(n, length(mt[[candidate]])))], tz = "GMT"))
    }
  }
  as.POSIXct(rep(NA_real_, n), origin = "1970-01-01", tz = "GMT")
}

result_points <- function(Result) {
  q <- Result$Results$Quantiles
  lon_col <- get_col(q, c("Medianlon", "MedianLon", "median_lon", "Meanlon", "MeanLon", "lon"))
  lat_col <- get_col(q, c("Medianlat", "MedianLat", "median_lat", "Meanlat", "MeanLat", "lat"))
  if (is.na(lon_col) || is.na(lat_col)) {
    stop("Could not find FLightR longitude/latitude columns in Result quantiles.",
         call. = FALSE)
  }
  n <- nrow(q)
  data.frame(
    twilight_index = seq_len(n),
    twilight_time = result_times(Result, n),
    lon = as.numeric(q[[lon_col]]),
    lat = as.numeric(q[[lat_col]])
  )
}

gps_seconds <- as.numeric(gps$time)

validate_one <- function(path) {
  meta <- metadata_for_result(path)
  Result <- readRDS(path)
  pts <- result_points(Result)
  q_seconds <- as.numeric(pts$twilight_time)
  ok_time <- is.finite(q_seconds)
  if (!any(ok_time)) {
    return(NULL)
  }
  nearest_diff <- vapply(q_seconds, function(x) min(abs(gps_seconds - x), na.rm = TRUE), numeric(1))
  gps_lon <- stats::approx(gps_seconds, gps$lon, xout = q_seconds, rule = 1)$y
  gps_lat <- stats::approx(gps_seconds, gps$lat, xout = q_seconds, rule = 1)$y
  keep <- ok_time & is.finite(gps_lon) & is.finite(gps_lat) &
    nearest_diff <= max_time_diff_hours * 3600
  if (!any(keep)) {
    return(NULL)
  }
  error_km <- geosphere::distGeo(
    cbind(pts$lon[keep], pts$lat[keep]),
    cbind(gps_lon[keep], gps_lat[keep])
  ) / 1000
  data.frame(
    run_label = meta$run_label[1],
    version = infer_version(meta$git_commit[1], meta$git_branch[1]),
    git_commit = meta$git_commit[1],
    git_branch = meta$git_branch[1],
    nParticles = meta$nParticles[1],
    seed = if ("seed" %in% names(meta)) meta$seed[1] else NA_integer_,
    subset = if ("subset" %in% names(meta)) meta$subset[1] else NA_character_,
    n_twilights = if ("nTwilights" %in% names(meta)) meta$nTwilights[1] else nrow(pts),
    result_file = basename(path),
    twilight_index = pts$twilight_index[keep],
    twilight_time = pts$twilight_time[keep],
    nearest_gps_time_diff_hours = nearest_diff[keep] / 3600,
    error_km = as.numeric(error_km)
  )
}

by_result <- do.call(rbind, lapply(result_files, validate_one))
if (is.null(by_result)) {
  by_result <- data.frame()
}

summary <- if (nrow(by_result) > 0) {
  do.call(rbind, lapply(split(by_result, by_result$result_file), function(x) {
    data.frame(
      run_label = x$run_label[1],
      version = x$version[1],
      git_commit = x$git_commit[1],
      git_branch = x$git_branch[1],
      nParticles = x$nParticles[1],
      seed = x$seed[1],
      subset = x$subset[1],
      n_twilights = x$n_twilights[1],
      n_matched = nrow(x),
      median_error_km = stats::median(x$error_km, na.rm = TRUE),
      mean_error_km = mean(x$error_km, na.rm = TRUE),
      q75_error_km = unname(stats::quantile(x$error_km, 0.75, na.rm = TRUE)),
      q95_error_km = unname(stats::quantile(x$error_km, 0.95, na.rm = TRUE)),
      max_error_km = max(x$error_km, na.rm = TRUE)
    )
  }))
} else {
  data.frame()
}

by_result_path <- file.path(output_dir, "gps_error_by_result.csv")
summary_path <- file.path(output_dir, "gps_error_summary_by_result.csv")
utils::write.csv(by_result, by_result_path, row.names = FALSE)
utils::write.csv(summary, summary_path, row.names = FALSE)

message("Validated result files: ", length(unique(by_result$result_file)))
print(summary)
message("Wrote: ", by_result_path)
message("Wrote: ", summary_path)

