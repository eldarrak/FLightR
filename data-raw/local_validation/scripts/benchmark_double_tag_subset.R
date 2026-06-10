# Local-only benchmark for private double-tagged FLightR data.
# Run this from a FLightR checkout with setwd("D:/GitHub/FLightR_old_master")
# or setwd("D:/GitHub/FLightR_new"), then source this script.
#
# Private data and outputs are intentionally outside Git tracking.

message("Starting local FLightR benchmark...")

# ----------------------------
# User-editable settings
# ----------------------------

private_root <- "D:/GitHub/FLightR/data-raw/local_validation/private_data"
output_root  <- "D:/GitHub/FLightR/data-raw/local_validation/outputs"

tags_file <- file.path(private_root, "A3_TAGS_format.csv")

# Optional GPS file. Leave as NA if not ready.
gps_file <- NA_character_
# Example:
# gps_file <- file.path(private_root, "A1_GPS_positions.csv")

start <- c(5.43, 52.93)
stop  <- c(5.43, 52.93)

# Important: for first subset tests, keep start.date = NA and shorten only end.date.
# This preserves early calibration data.
proc_start_date <- Sys.getenv("FLIGHTR_PROC_START_DATE", unset = NA)
proc_end_date   <- Sys.getenv("FLIGHTR_PROC_END_DATE", unset = "2014-01-15")
proc_start_date <- if (identical(tolower(proc_start_date), "na") || !nzchar(proc_start_date)) NA else proc_start_date
proc_end_date <- if (identical(tolower(proc_end_date), "na") || !nzchar(proc_end_date)) NA else proc_end_date

model.ageing <- FALSE
likelihood.correction <- FALSE    # use FALSE for smoke tests; later try "auto" or TRUE

grid_left <- as.numeric(Sys.getenv("FLIGHTR_GRID_LEFT", unset = "-14"))
grid_right <- as.numeric(Sys.getenv("FLIGHTR_GRID_RIGHT", unset = "13"))
grid_bottom <- as.numeric(Sys.getenv("FLIGHTR_GRID_BOTTOM", unset = "30"))
grid_top <- as.numeric(Sys.getenv("FLIGHTR_GRID_TOP", unset = "57"))

nParticles <- as.integer(Sys.getenv("FLIGHTR_NPARTICLES", unset = "100"))
legacy_threads <- Sys.getenv("FLIGHTR_THREADS", unset = NA)
prerun_threads <- as.integer(Sys.getenv("FLIGHTR_PRERUN_THREADS", unset = ifelse(is.na(legacy_threads), "4", legacy_threads)))
pf_threads <- as.integer(Sys.getenv("FLIGHTR_PF_THREADS", unset = ifelse(is.na(legacy_threads), "1", legacy_threads)))
# Prerun/preparation parallelism can be beneficial; particle-filter PSOCK parallelism should be benchmarked before use.

known.last <- tolower(Sys.getenv("FLIGHTR_KNOWN_LAST", unset = "false")) %in% c("1", "true", "yes")
# For a truncated subset, use known.last = FALSE unless the subset ends at a real known location.
# For full track ending at recapture, set TRUE.
precision.sd <- 25
check.outliers <- FALSE

random_seed <- as.integer(Sys.getenv("FLIGHTR_SEED", unset = "123"))
profile.phases <- tolower(Sys.getenv("FLIGHTR_PROFILE_PHASES", unset = "false")) %in% c("1", "true", "yes")
profile.top.level <- tolower(Sys.getenv("FLIGHTR_PROFILE_TOP_LEVEL", unset = "false")) %in% c("1", "true", "yes")
propagation.backend <- Sys.getenv("FLIGHTR_PROPAGATION_BACKEND", unset = "auto")

run_label <- Sys.getenv("FLIGHTR_RUN_LABEL", unset = NA)
if (is.na(run_label) || !nzchar(run_label)) {
  run_label <- paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"))
}

# ----------------------------
# Helpers
# ----------------------------

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
source("D:/GitHub/FLightR/data-raw/local_validation/scripts/prerun_cache_helpers.R")

repo_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

git_value <- function(args) {
  out <- tryCatch(
    system2("git", args, stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_
  )
  if (length(out) == 0) NA_character_ else out[[1]]
}

git_commit <- git_value(c("rev-parse", "--short", "HEAD"))
git_branch <- git_value(c("branch", "--show-current"))
git_status_porcelain <- git_value(c("status", "--porcelain"))
git_dirty <- !is.na(git_status_porcelain) && nzchar(git_status_porcelain)

safe_label <- gsub("[^A-Za-z0-9_\\-]+", "_", run_label)
prefix <- paste0(
  safe_label,
  "_", ifelse(is.na(git_branch), "unknown_branch", git_branch),
  "_", ifelse(is.na(git_commit), "unknown_commit", git_commit),
  "_Np", nParticles
)

time_it <- function(expr) {
  gc()
  start_time <- Sys.time()
  start_proc <- proc.time()
  value <- eval.parent(substitute(expr))
  end_proc <- proc.time()
  end_time <- Sys.time()
  list(
    value = value,
    elapsed_seconds = unname((end_proc - start_proc)[["elapsed"]]),
    wall_seconds = as.numeric(difftime(end_time, start_time, units = "secs"))
  )
}

as_posix_or_null <- function(x) {
  if (length(x) == 0 || is.na(x) || !nzchar(x)) return(NULL)
  as.POSIXct(x, tz = "GMT")
}

load_local_flightR <- function(repo_root) {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(repo_root, quiet = TRUE, export_all = FALSE)
    return("pkgload::load_all")
  }
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(repo_root, quiet = TRUE, export_all = FALSE)
    return("devtools::load_all")
  }
  library(FLightR)
  return("library(FLightR)")
}

distance_km <- function(lon1, lat1, lon2, lat2) {
  if (!requireNamespace("geosphere", quietly = TRUE)) {
    stop("Package 'geosphere' is needed for distance calculations.")
  }
  geosphere::distGeo(cbind(lon1, lat1), cbind(lon2, lat2)) / 1000
}

guess_column <- function(x, candidates) {
  hit <- candidates[candidates %in% names(x)]
  if (length(hit) == 0) NA_character_ else hit[[1]]
}

read_gps_optional <- function(path) {
  if (length(path) == 0 || is.na(path) || !nzchar(path) || !file.exists(path)) {
    return(NULL)
  }
  
  gps <- utils::read.csv(path, stringsAsFactors = FALSE)
  
  time_col <- guess_column(gps, c(
    "datetime", "timestamp", "time", "Time", "date", "Date",
    "gmt", "GMT", "DateTime", "date_time"
  ))
  lon_col <- guess_column(gps, c(
    "lon", "Lon", "longitude", "Longitude", "LONGITUDE", "x"
  ))
  lat_col <- guess_column(gps, c(
    "lat", "Lat", "latitude", "Latitude", "LATITUDE", "y"
  ))
  
  if (any(is.na(c(time_col, lon_col, lat_col)))) {
    warning("GPS file was found, but columns could not be guessed. ",
            "Need time/lon/lat columns. Skipping GPS validation.")
    return(NULL)
  }
  
  gps_out <- data.frame(
    time = as.POSIXct(gps[[time_col]], tz = "GMT"),
    lon = as.numeric(gps[[lon_col]]),
    lat = as.numeric(gps[[lat_col]])
  )
  gps_out <- gps_out[is.finite(gps_out$lon) & is.finite(gps_out$lat) & !is.na(gps_out$time), ]
  gps_out <- gps_out[order(gps_out$time), ]
  gps_out
}

match_gps_to_result <- function(Result, gps, max_time_diff_hours = 6) {
  if (is.null(gps) || nrow(gps) == 0) return(NULL)
  
  q <- Result$Results$Quantiles
  if (!("time" %in% names(q))) {
    if ("time" %in% names(Result$Indices$Matrix.Index.Table)) {
      q$time <- Result$Indices$Matrix.Index.Table$time[seq_len(nrow(q))]
    } else {
      warning("No time column available in Result quantiles. Skipping GPS validation.")
      return(NULL)
    }
  }
  
  out <- vector("list", nrow(q))
  for (i in seq_len(nrow(q))) {
    dt <- abs(as.numeric(difftime(gps$time, q$time[i], units = "hours")))
    j <- which.min(dt)
    if (length(j) == 1 && is.finite(dt[j]) && dt[j] <= max_time_diff_hours) {
      lon_est <- if ("Medianlon" %in% names(q)) q$Medianlon[i] else q$Meanlon[i]
      lat_est <- if ("Medianlat" %in% names(q)) q$Medianlat[i] else q$Meanlat[i]
      out[[i]] <- data.frame(
        twilight_index = i,
        flightr_time = q$time[i],
        gps_time = gps$time[j],
        time_diff_hours = dt[j],
        est_lon = lon_est,
        est_lat = lat_est,
        gps_lon = gps$lon[j],
        gps_lat = gps$lat[j],
        error_km = distance_km(lon_est, lat_est, gps$lon[j], gps$lat[j])
      )
    }
  }
  
  out <- do.call(rbind, out)
  if (is.null(out) || nrow(out) == 0) return(NULL)
  out
}

summarise_gps_validation <- function(matched) {
  if (is.null(matched) || nrow(matched) == 0) {
    return(data.frame(
      n_matched = 0,
      median_error_km = NA_real_,
      mean_error_km = NA_real_,
      q95_error_km = NA_real_
    ))
  }
  
  data.frame(
    n_matched = nrow(matched),
    median_error_km = stats::median(matched$error_km, na.rm = TRUE),
    mean_error_km = mean(matched$error_km, na.rm = TRUE),
    q95_error_km = unname(stats::quantile(matched$error_km, 0.95, na.rm = TRUE))
  )
}

# ----------------------------
# Run benchmark
# ----------------------------

set.seed(random_seed)

load_method <- load_local_flightR(repo_root)
message("Loaded FLightR using: ", load_method)
message("Repo: ", repo_root)
message("Branch: ", git_branch, " commit: ", git_commit)
message("Run label: ", run_label)
message("Particles: ", nParticles, " prerun threads: ", prerun_threads, " PF threads: ", pf_threads)
message("Propagation backend: ", propagation.backend)

if (!file.exists(tags_file)) {
  stop("TAGS file does not exist: ", tags_file)
}

start_date <- as_posix_or_null(proc_start_date)
end_date <- as_posix_or_null(proc_end_date)

proc_step <- time_it({
  get.tags.data(tags_file, start.date = start_date, end.date = end_date)
})
Proc.data <- proc_step$value

Calibration.periods <- data.frame(
  calibration.start = as.POSIXct(c("2000-01-01", "2014-05-05"), tz = "GMT"),
  calibration.stop  = as.POSIXct(c("2013-08-20", "2020-01-01"), tz = "GMT"),
  lon = start[1],
  lat = start[2]
)

cal_step <- time_it({
  make.calibration(
    Proc.data,
    Calibration.periods,
    model.ageing = model.ageing,
    plot.each = FALSE,
    plot.final = FALSE,
    likelihood.correction = likelihood.correction
  )
})
Calibration <- cal_step$value

grid_step <- time_it({
  make.grid(
    left = grid_left,
    bottom = grid_bottom,
    right = grid_right,
    top = grid_top,
    distance.from.land.allowed.to.use = c(-Inf, Inf),
    distance.from.land.allowed.to.stay = c(-Inf, Inf),
    plot = FALSE
  )
})
Grid <- grid_step$value

prerun_key_components <- list(
  cache_schema = "local_validation_prerun_v1",
  tags_file = as.list(file_fingerprint(tags_file)),
  proc_start_date = ifelse(is.null(start_date), NA, format(start_date, "%Y-%m-%d %H:%M:%S %Z")),
  proc_end_date = ifelse(is.null(end_date), NA, format(end_date, "%Y-%m-%d %H:%M:%S %Z")),
  calibration_periods = paste(utils::capture.output(print(Calibration.periods)), collapse = " | "),
  model_ageing = model.ageing,
  likelihood_correction = as.character(likelihood.correction),
  grid_left = grid_left,
  grid_right = grid_right,
  grid_bottom = grid_bottom,
  grid_top = grid_top,
  grid_size = nrow(Grid),
  distance_from_land_allowed_to_use = "-Inf;Inf",
  distance_from_land_allowed_to_stay = "-Inf;Inf",
  start = paste(start, collapse = ";"),
  stop = paste(stop, collapse = ";"),
  make_prerun_threads = prerun_threads,
  propagation_backend_in_prerun = "none_currently_backend_independent",
  git_commit = git_commit,
  git_dirty = git_dirty,
  FLightR_version = as.character(utils::packageVersion("FLightR"))
)
prerun_key <- make_prerun_cache_key(prerun_key_components)
prerun_cached <- load_or_build_prerun(
  key = prerun_key,
  expected = list(grid_size = nrow(Grid)),
  output_root = output_root,
  build_fun = function() {
    make.prerun.object(
      Proc.data,
      Grid,
      start = start,
      end = stop,
      Calibration = Calibration,
      threads = prerun_threads
    )
  }
)
all.in <- prerun_cached$value
prerun_cache <- prerun_cached$cache
prerun_step <- list(
  value = all.in,
  elapsed_seconds = prerun_cache$effective_seconds,
  wall_seconds = prerun_cache$effective_seconds
)
message("Prerun cache enabled: ", prerun_cache$enabled)
message("Prerun cache hit: ", prerun_cache$hit)
message("Prerun cache key: ", prerun_cache$key)
message("Prerun cache file: ", prerun_cache$file)
if (!is.na(prerun_cache$warning)) warning(prerun_cache$warning, call. = FALSE)

set.seed(random_seed)  # Keep PF RNG comparable across prerun cache rebuild/hit paths.
pf_step <- time_it({
  run.particle.filter(
    all.in,
    threads = pf_threads,
    nParticles = nParticles,
    known.last = known.last,
    precision.sd = precision.sd,
    check.outliers = check.outliers,
    profile.phases = profile.phases,
    profile.top.level = profile.top.level,
    propagation.backend = propagation.backend,
    plot = FALSE
  )
})
Result <- pf_step$value

result_path <- file.path(output_root, paste0("Result_", prefix, ".rds"))
saveRDS(Result, result_path)

phase_profile_path <- NA_character_
if (isTRUE(profile.phases) && !is.null(Result$Results$phase_profile)) {
  phase_profile_path <- file.path(output_root, paste0("pf_phase_profile_", prefix, ".csv"))
  utils::write.csv(Result$Results$phase_profile, phase_profile_path, row.names = FALSE)
}

top_level_profile <- data.frame(
  phase = c(
    "get.tags.data",
    "make.calibration",
    "make.grid",
    "make.prerun.object",
    "run.particle.filter",
    if (!is.null(Result$Results$top_level_profile)) {
      paste0("run.particle.filter::", Result$Results$top_level_profile$phase)
    } else {
      character()
    },
    "total_benchmark_runtime"
  ),
  elapsed_seconds = c(
    proc_step$elapsed_seconds,
    cal_step$elapsed_seconds,
    grid_step$elapsed_seconds,
    prerun_step$elapsed_seconds,
    pf_step$elapsed_seconds,
    if (!is.null(Result$Results$top_level_profile)) {
      Result$Results$top_level_profile$elapsed_seconds
    } else {
      numeric()
    },
    proc_step$elapsed_seconds + cal_step$elapsed_seconds +
      grid_step$elapsed_seconds + prerun_step$elapsed_seconds + pf_step$elapsed_seconds
  ),
  stringsAsFactors = FALSE
)

top_level_profile_path <- file.path(output_root, paste0("pf_top_level_timing_", prefix, ".csv"))
utils::write.csv(top_level_profile, top_level_profile_path, row.names = FALSE)

benchmark <- data.frame(
  run_label = run_label,
  git_branch = git_branch,
  git_commit = git_commit,
  repo_root = repo_root,
  tags_file_basename = basename(tags_file),
  proc_start_date = ifelse(is.null(start_date), NA, format(start_date, "%Y-%m-%d %H:%M:%S")),
  proc_end_date = ifelse(is.null(end_date), NA, format(end_date, "%Y-%m-%d %H:%M:%S")),
  nParticles = nParticles,
  threads = pf_threads,
  prerun_threads = prerun_threads,
  pf_threads = pf_threads,
  legacy_threads_env = ifelse(is.na(legacy_threads), NA, legacy_threads),
  known_last = known.last,
  precision_sd = precision.sd,
  check_outliers = check.outliers,
  profile_phases = profile.phases,
  profile_top_level = profile.top.level,
  propagation_backend = propagation.backend,
  model_ageing = model.ageing,
  likelihood_correction = as.character(likelihood.correction),
  grid_left = grid_left,
  grid_right = grid_right,
  grid_bottom = grid_bottom,
  grid_top = grid_top,
  grid_size = nrow(Grid),
  n_twilights = nrow(Result$Indices$Matrix.Index.Table),
  LL = if (!is.null(Result$Results$LL)) Result$Results$LL else NA_real_,
  get_tags_seconds = proc_step$elapsed_seconds,
  calibration_seconds = cal_step$elapsed_seconds,
  grid_seconds = grid_step$elapsed_seconds,
  prerun_seconds = prerun_step$elapsed_seconds,
  prerun_cache_enabled = prerun_cache$enabled,
  prerun_cache_hit = prerun_cache$hit,
  prerun_cache_file = prerun_cache$file,
  prerun_cache_key = prerun_cache$key,
  prerun_build_seconds = prerun_cache$build_seconds,
  prerun_load_seconds = prerun_cache$load_seconds,
  prerun_save_seconds = prerun_cache$save_seconds,
  make_prerun_seconds_effective = prerun_cache$effective_seconds,
  particle_filter_seconds = pf_step$elapsed_seconds,
  total_seconds = proc_step$elapsed_seconds + cal_step$elapsed_seconds +
    grid_step$elapsed_seconds + prerun_step$elapsed_seconds + pf_step$elapsed_seconds,
  phase_profile_path = phase_profile_path,
  top_level_profile_path = top_level_profile_path,
  result_path = result_path
)

benchmark_path <- file.path(output_root, paste0("benchmark_", prefix, ".csv"))
utils::write.csv(benchmark, benchmark_path, row.names = FALSE)

gps <- read_gps_optional(gps_file)
matched <- match_gps_to_result(Result, gps, max_time_diff_hours = 6)
gps_summary <- summarise_gps_validation(matched)

gps_summary_path <- file.path(output_root, paste0("gps_summary_", prefix, ".csv"))
utils::write.csv(gps_summary, gps_summary_path, row.names = FALSE)

if (!is.null(matched)) {
  gps_matched_path <- file.path(output_root, paste0("gps_matched_", prefix, ".csv"))
  utils::write.csv(matched, gps_matched_path, row.names = FALSE)
} else {
  gps_matched_path <- NA_character_
}

message("DONE")
message("Result saved to: ", result_path)
message("Benchmark saved to: ", benchmark_path)
if (!is.na(phase_profile_path)) message("Phase profile saved to: ", phase_profile_path)
message("Top-level timing saved to: ", top_level_profile_path)
message("GPS summary saved to: ", gps_summary_path)
if (!is.na(gps_matched_path)) message("GPS matched errors saved to: ", gps_matched_path)

invisible(list(
  Result = Result,
  benchmark = benchmark,
  gps_summary = gps_summary,
  result_path = result_path,
  benchmark_path = benchmark_path,
  phase_profile_path = phase_profile_path,
  top_level_profile_path = top_level_profile_path,
  gps_summary_path = gps_summary_path,
  gps_matched_path = gps_matched_path
))
