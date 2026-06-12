#!/usr/bin/env Rscript

root_dir <- normalizePath(file.path(getwd()), winslash = "/", mustWork = FALSE)
config_path <- file.path(root_dir, "data-raw/local_validation/validation_config.local.yml")

if (!file.exists(config_path)) {
  stop(
    "Local validation config not found: ", config_path, "\n",
    "Copy data-raw/local_validation/validation_config.example.yml to ",
    "data-raw/local_validation/validation_config.local.yml and fill in private local paths.",
    call. = FALSE
  )
}

source(file.path(root_dir, "R/validation_metrics.R"))

config <- read_validation_config(config_path)
output_dir <- config$outputs$output_dir
if (is.null(output_dir)) {
  output_dir <- "data-raw/local_validation/outputs"
}

output_dir <- normalizePath(file.path(root_dir, output_dir), winslash = "/", mustWork = FALSE)
allowed_output_dir <- normalizePath(
  file.path(root_dir, "data-raw/local_validation/outputs"),
  winslash = "/",
  mustWork = FALSE
)

if (!startsWith(output_dir, allowed_output_dir)) {
  stop("Validation outputs must be saved under data-raw/local_validation/outputs/.", call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

gps_track_path <- config$inputs$gps_track_path
if (is.null(gps_track_path)) {
  stop("Config must define inputs$gps_track_path.", call. = FALSE)
}

gps_track_path <- normalizePath(file.path(root_dir, gps_track_path), winslash = "/", mustWork = TRUE)

result_object_name <- config$inputs$result_object_name
result_path <- config$inputs$result_path

if (!is.null(result_object_name) && exists(result_object_name, envir = .GlobalEnv)) {
  Result <- get(result_object_name, envir = .GlobalEnv)
} else if (!is.null(result_path)) {
  result_path <- normalizePath(file.path(root_dir, result_path), winslash = "/", mustWork = TRUE)
  Result <- readRDS(result_path)
} else {
  stop(
    "Config must define inputs$result_path, or inputs$result_object_name must name an existing Result object.",
    call. = FALSE
  )
}

gps <- read_gps_track(gps_track_path)

max_time_diff_hours <- config$matching$max_time_diff_hours
if (is.null(max_time_diff_hours)) {
  max_time_diff_hours <- 6
}

gps_matched <- match_gps_to_twilights(
  Result,
  gps,
  max_time_diff_hours = as.numeric(max_time_diff_hours)
)
summary <- validation_summary(Result, gps_matched)

summary_name <- config$outputs$summary_csv
if (is.null(summary_name)) {
  summary_name <- "double_tag_validation_summary.csv"
}

summary_path <- file.path(output_dir, basename(summary_name))
matched_path <- file.path(output_dir, "double_tag_validation_matches.csv")

utils::write.csv(summary, summary_path, row.names = FALSE)
utils::write.csv(gps_matched, matched_path, row.names = FALSE)

message("Validation summary written to: ", summary_path)
message("Matched GPS fixes written to: ", matched_path)
