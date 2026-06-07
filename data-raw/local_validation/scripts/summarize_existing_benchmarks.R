output_dir <- "D:/GitHub/FLightR/data-raw/local_validation/outputs"

if (!dir.exists(output_dir)) {
  stop("Output directory not found: ", output_dir, call. = FALSE)
}

benchmark_files <- list.files(
  output_dir,
  pattern = "^benchmark_.*\\.csv$",
  full.names = TRUE
)
benchmark_files <- benchmark_files[
  !grepl("benchmark_(combined|runtime_summary|label_warnings|old_new_comparison)", basename(benchmark_files))
]

if (length(benchmark_files) == 0) {
  stop("No raw benchmark_*.csv files found in ", output_dir, call. = FALSE)
}

read_benchmark <- function(path) {
  x <- utils::read.csv(path, stringsAsFactors = FALSE)
  x$benchmark_file <- basename(path)
  x
}

bench <- do.call(rbind, lapply(benchmark_files, read_benchmark))

infer_version <- function(git_commit, git_branch) {
  if (!is.na(git_commit) && git_commit == "86fcb29") {
    return("old_master")
  }
  if ((!is.na(git_commit) && git_commit == "fe959c1") ||
      (!is.na(git_branch) && git_branch == "local-validation-framework")) {
    return("new_branch")
  }
  "unknown"
}

bench$version <- mapply(infer_version, bench$git_commit, bench$git_branch)
bench$label_mentions_old <- grepl("old", bench$run_label, ignore.case = TRUE)
bench$label_mentions_new <- grepl("new", bench$run_label, ignore.case = TRUE)
bench$label_warning <- ""
bench$label_warning[bench$label_mentions_old & bench$version == "new_branch"] <-
  "run_label says old but git metadata says new_branch"
bench$label_warning[bench$label_mentions_new & bench$version == "old_master"] <-
  "run_label says new but git metadata says old_master"

parse_seed <- function(label) {
  hit <- regmatches(label, regexpr("(seed|s)[_-]?[0-9]+", label, ignore.case = TRUE))
  if (length(hit) == 0 || hit == "") {
    return(NA_integer_)
  }
  as.integer(gsub("[^0-9]", "", hit))
}

parse_subset <- function(label) {
  known <- c("fullyear_knownlast", "halfyear", "short", "full", "smoke")
  for (candidate in known) {
    if (grepl(candidate, label, ignore.case = TRUE)) {
      return(candidate)
    }
  }
  NA_character_
}

bench$seed <- vapply(bench$run_label, parse_seed, integer(1))
bench$seed_group <- ifelse(is.na(bench$seed), "NA", as.character(bench$seed))
bench$subset <- vapply(bench$run_label, parse_subset, character(1))
bench$n_twilights <- if ("n_twilights" %in% names(bench)) bench$n_twilights else NA_integer_
bench$nTwilights <- if ("nTwilights" %in% names(bench)) bench$nTwilights else bench$n_twilights

runtime_summary <- aggregate(
  cbind(
    get_tags_seconds,
    calibration_seconds,
    grid_seconds,
    prerun_seconds,
    particle_filter_seconds,
    total_seconds
  ) ~ version + nParticles + nTwilights + seed_group,
  data = bench,
  FUN = function(x) mean(x, na.rm = TRUE),
  na.action = NULL
)
names(runtime_summary)[names(runtime_summary) == "seed_group"] <- "seed"

label_warnings <- bench[bench$label_warning != "", c(
  "run_label", "version", "git_commit", "git_branch", "benchmark_file", "label_warning"
)]
if (nrow(label_warnings) == 0) {
  label_warnings <- data.frame(
    run_label = character(),
    version = character(),
    git_commit = character(),
    git_branch = character(),
    benchmark_file = character(),
    label_warning = character()
  )
}

combined_path <- file.path(output_dir, "benchmark_combined.csv")
runtime_path <- file.path(output_dir, "benchmark_runtime_summary.csv")
warnings_path <- file.path(output_dir, "benchmark_label_warnings.csv")

utils::write.csv(bench, combined_path, row.names = FALSE)
utils::write.csv(runtime_summary, runtime_path, row.names = FALSE)
utils::write.csv(label_warnings, warnings_path, row.names = FALSE)

message("Benchmark rows: ", nrow(bench))
message("Label warnings: ", nrow(label_warnings))
print(runtime_summary)
message("Wrote: ", combined_path)
message("Wrote: ", runtime_path)
message("Wrote: ", warnings_path)
