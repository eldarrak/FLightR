# Summarise local paired FLightR benchmark runs.
# Reads benchmark and optional GPS summary files from outputs only.

output_dir <- "D:/GitHub/FLightR/data-raw/local_validation/outputs"

if (!dir.exists(output_dir)) {
  stop("Output directory not found: ", output_dir, call. = FALSE)
}

benchmark_files <- list.files(output_dir, pattern = "^benchmark_.*\\.csv$", full.names = TRUE)
benchmark_files <- benchmark_files[!grepl("^benchmark_(combined|old_new_comparison|runtime_summary|label_warnings)", basename(benchmark_files))]

read_one <- function(path) {
  out <- tryCatch(utils::read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(out) || nrow(out) == 0 || !"run_label" %in% names(out)) return(NULL)
  out$benchmark_file <- path
  out
}

bench <- do.call(rbind, lapply(benchmark_files, read_one))
if (is.null(bench) || nrow(bench) == 0) {
  stop("No benchmark CSV rows with run_label found in: ", output_dir, call. = FALSE)
}

parse_label <- function(label) {
  m <- regexec("^(baseline|s2)_halfyear_([^_]+)_(.+)_seed([0-9]+)$", label)
  hit <- regmatches(label, m)[[1]]
  if (length(hit) == 0) {
    return(data.frame(
      paired_version = NA_character_,
      particle_label = NA_character_,
      pair_id = NA_character_,
      paired_seed = NA_integer_
    ))
  }
  data.frame(
    paired_version = hit[2],
    particle_label = hit[3],
    pair_id = hit[4],
    paired_seed = as.integer(hit[5]),
    stringsAsFactors = FALSE
  )
}

parsed <- do.call(rbind, lapply(bench$run_label, parse_label))
bench <- cbind(bench, parsed)
bench <- bench[bench$paired_version %in% c("baseline", "s2"), , drop = FALSE]

if (nrow(bench) == 0) {
  stop("No paired benchmark labels found. Expected labels like baseline_halfyear_1e5_pair1_seed123.",
       call. = FALSE)
}

num_col <- function(x, name) {
  if (name %in% names(x)) as.numeric(x[[name]]) else rep(NA_real_, nrow(x))
}

runtime_summary <- data.frame(
  run_label = bench$run_label,
  paired_version = bench$paired_version,
  pair_id = bench$pair_id,
  seed = bench$paired_seed,
  particle_label = bench$particle_label,
  nParticles = num_col(bench, "nParticles"),
  threads = num_col(bench, "threads"),
  n_twilights = if ("n_twilights" %in% names(bench)) num_col(bench, "n_twilights") else num_col(bench, "nTwilights"),
  grid_size = num_col(bench, "grid_size"),
  LL = num_col(bench, "LL"),
  particle_filter_seconds = num_col(bench, "particle_filter_seconds"),
  total_seconds = num_col(bench, "total_seconds"),
  git_branch = if ("git_branch" %in% names(bench)) bench$git_branch else NA_character_,
  git_commit = if ("git_commit" %in% names(bench)) bench$git_commit else NA_character_,
  benchmark_file = bench$benchmark_file,
  stringsAsFactors = FALSE
)

find_gps_summary <- function(run_label, nParticles) {
  pattern <- paste0("^gps_summary_", run_label, ".*Np", nParticles, ".*\\.csv$")
  files <- list.files(output_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) return(NULL)
  out <- tryCatch(utils::read.csv(files[1], stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(out) || nrow(out) == 0) return(NULL)
  out$gps_summary_file <- files[1]
  out
}

gps_rows <- lapply(seq_len(nrow(runtime_summary)), function(i) {
  gps <- find_gps_summary(runtime_summary$run_label[i], runtime_summary$nParticles[i])
  if (is.null(gps)) return(NULL)
  data.frame(
    run_label = runtime_summary$run_label[i],
    median_error_km = if ("median_error_km" %in% names(gps)) as.numeric(gps$median_error_km[1]) else NA_real_,
    mean_error_km = if ("mean_error_km" %in% names(gps)) as.numeric(gps$mean_error_km[1]) else NA_real_,
    q75_error_km = if ("q75_error_km" %in% names(gps)) as.numeric(gps$q75_error_km[1]) else NA_real_,
    q95_error_km = if ("q95_error_km" %in% names(gps)) as.numeric(gps$q95_error_km[1]) else NA_real_,
    gps_summary_file = gps$gps_summary_file[1],
    stringsAsFactors = FALSE
  )
})
gps_summary <- if (length(Filter(Negate(is.null), gps_rows)) > 0) {
  do.call(rbind, gps_rows)
} else {
  data.frame()
}

if (nrow(gps_summary) > 0) {
  runtime_summary <- merge(runtime_summary, gps_summary, by = "run_label", all.x = TRUE)
}

pair_key <- paste(runtime_summary$nParticles, runtime_summary$seed, runtime_summary$pair_id, sep = "|")
runtime_summary$pair_key <- pair_key

baseline <- runtime_summary[runtime_summary$paired_version == "baseline", , drop = FALSE]
s2 <- runtime_summary[runtime_summary$paired_version == "s2", , drop = FALSE]

pairs <- merge(baseline, s2, by = "pair_key", suffixes = c("_baseline", "_s2"))

comparison <- if (nrow(pairs) > 0) {
  data.frame(
    pair_id = pairs$pair_id_baseline,
    seed = pairs$seed_baseline,
    nParticles = pairs$nParticles_baseline,
    particle_label = pairs$particle_label_baseline,
    baseline_run_label = pairs$run_label_baseline,
    s2_run_label = pairs$run_label_s2,
    baseline_particle_filter_seconds = pairs$particle_filter_seconds_baseline,
    s2_particle_filter_seconds = pairs$particle_filter_seconds_s2,
    speedup_pf = pairs$particle_filter_seconds_baseline / pairs$particle_filter_seconds_s2,
    baseline_total_seconds = pairs$total_seconds_baseline,
    s2_total_seconds = pairs$total_seconds_s2,
    speedup_total = pairs$total_seconds_baseline / pairs$total_seconds_s2,
    baseline_LL = pairs$LL_baseline,
    s2_LL = pairs$LL_s2,
    delta_LL = pairs$LL_s2 - pairs$LL_baseline,
    delta_median_error_km = if ("median_error_km_s2" %in% names(pairs)) pairs$median_error_km_s2 - pairs$median_error_km_baseline else NA_real_,
    delta_mean_error_km = if ("mean_error_km_s2" %in% names(pairs)) pairs$mean_error_km_s2 - pairs$mean_error_km_baseline else NA_real_,
    delta_q75_error_km = if ("q75_error_km_s2" %in% names(pairs)) pairs$q75_error_km_s2 - pairs$q75_error_km_baseline else NA_real_,
    delta_q95_error_km = if ("q95_error_km_s2" %in% names(pairs)) pairs$q95_error_km_s2 - pairs$q95_error_km_baseline else NA_real_,
    stringsAsFactors = FALSE
  )
} else {
  data.frame()
}

summary_path <- file.path(output_dir, "paired_runtime_summary.csv")
comparison_path <- file.path(output_dir, "paired_runtime_comparison.csv")

utils::write.csv(runtime_summary, summary_path, row.names = FALSE)
utils::write.csv(comparison, comparison_path, row.names = FALSE)

message("Wrote: ", summary_path)
message("Wrote: ", comparison_path)

if (nrow(comparison) == 0) {
  message("No complete baseline/s2 pairs found.")
} else {
  print(comparison[, c(
    "pair_id",
    "seed",
    "nParticles",
    "speedup_pf",
    "speedup_total",
    "delta_LL",
    "delta_median_error_km",
    "delta_mean_error_km",
    "delta_q75_error_km",
    "delta_q95_error_km"
  ), drop = FALSE])
}
