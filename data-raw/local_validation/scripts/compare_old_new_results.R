output_dir <- "D:/GitHub/FLightR/data-raw/local_validation/outputs"

required <- c(
  "benchmark_combined.csv",
  "gps_error_by_result.csv",
  "gps_error_summary_by_result.csv"
)
missing <- required[!file.exists(file.path(output_dir, required))]
if (length(missing) > 0) {
  stop("Missing required summary files: ", paste(missing, collapse = ", "),
       call. = FALSE)
}

bench <- utils::read.csv(file.path(output_dir, "benchmark_combined.csv"), stringsAsFactors = FALSE)
gps_by <- utils::read.csv(file.path(output_dir, "gps_error_by_result.csv"), stringsAsFactors = FALSE)
gps_sum <- utils::read.csv(file.path(output_dir, "gps_error_summary_by_result.csv"), stringsAsFactors = FALSE)

if (!"subset" %in% names(bench)) bench$subset <- NA_character_
if (!"seed" %in% names(bench)) bench$seed <- NA_integer_
if (!"label_warning" %in% names(bench)) bench$label_warning <- ""
bench$label_warning[is.na(bench$label_warning)] <- ""

valid_bench <- bench[bench$label_warning == "" & bench$version %in% c("old_master", "new_branch"), ]
valid_gps <- gps_sum[gps_sum$version %in% c("old_master", "new_branch"), ]

pair_key <- function(x) paste(x$nParticles, x$seed, x$subset, x$nTwilights, sep = "|")
valid_bench$key <- pair_key(valid_bench)

old_bench <- valid_bench[valid_bench$version == "old_master", ]
new_bench <- valid_bench[valid_bench$version == "new_branch", ]
pairs <- merge(old_bench, new_bench, by = "key", suffixes = c("_old", "_new"))

if (nrow(pairs) == 0) {
  old_bench$obvious <- sub("^old_master_", "", old_bench$run_label)
  new_bench$obvious <- sub("^new_branch_", "", new_bench$run_label)
  pairs <- merge(old_bench, new_bench, by = "obvious", suffixes = c("_old", "_new"))
}

gps_key_cols <- c("run_label", "version", "git_commit", "git_branch", "nParticles",
                  "seed", "subset", "n_twilights")
summary_rows <- list()
twilight_rows <- list()

if (nrow(pairs) > 0) {
  for (i in seq_len(nrow(pairs))) {
    old_label <- pairs$run_label_old[i]
    new_label <- pairs$run_label_new[i]
    old_gps <- valid_gps[valid_gps$run_label == old_label, , drop = FALSE]
    new_gps <- valid_gps[valid_gps$run_label == new_label, , drop = FALSE]
    if (nrow(old_gps) == 0 || nrow(new_gps) == 0) next
    pair_id <- paste(old_label, "vs", new_label)
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      pair_id = pair_id,
      old_run_label = old_label,
      new_run_label = new_label,
      nParticles = pairs$nParticles_old[i],
      seed = pairs$seed_old[i],
      subset = pairs$subset_old[i],
      n_twilights = pairs$nTwilights_old[i],
      old_git_commit = pairs$git_commit_old[i],
      new_git_commit = pairs$git_commit_new[i],
      old_particle_filter_seconds = pairs$particle_filter_seconds_old[i],
      new_particle_filter_seconds = pairs$particle_filter_seconds_new[i],
      old_total_seconds = pairs$total_seconds_old[i],
      new_total_seconds = pairs$total_seconds_new[i],
      speedup_pf = pairs$particle_filter_seconds_old[i] / pairs$particle_filter_seconds_new[i],
      speedup_total = pairs$total_seconds_old[i] / pairs$total_seconds_new[i],
      delta_LL = pairs$LL_new[i] - pairs$LL_old[i],
      delta_median_error_km = new_gps$median_error_km[1] - old_gps$median_error_km[1],
      delta_mean_error_km = new_gps$mean_error_km[1] - old_gps$mean_error_km[1],
      delta_q75_error_km = new_gps$q75_error_km[1] - old_gps$q75_error_km[1],
      delta_q95_error_km = new_gps$q95_error_km[1] - old_gps$q95_error_km[1],
      delta_max_error_km = new_gps$max_error_km[1] - old_gps$max_error_km[1]
    )
    old_by <- gps_by[gps_by$run_label == old_label, ]
    new_by <- gps_by[gps_by$run_label == new_label, ]
    merged <- merge(
      old_by,
      new_by,
      by = "twilight_index",
      suffixes = c("_old", "_new")
    )
    if (nrow(merged) == 0) next
    merged$pair_id <- pair_id
    merged$twilight_time <- as.POSIXct(merged$twilight_time_old, tz = "GMT")
    merged$delta_error_km <- merged$error_km_new - merged$error_km_old
    merged$period <- ifelse(
      merged$twilight_time < as.POSIXct("2013-08-25", tz = "GMT"),
      "pre_movement",
      "movement_and_winter"
    )
    twilight_rows[[length(twilight_rows) + 1]] <- merged[, c(
      "pair_id", "twilight_index", "twilight_time", "period",
      "error_km_old", "error_km_new", "delta_error_km"
    )]
  }
}

comparison <- if (length(summary_rows)) do.call(rbind, summary_rows) else data.frame()
by_twilight <- if (length(twilight_rows)) do.call(rbind, twilight_rows) else data.frame()

by_period <- if (nrow(by_twilight) > 0) {
  do.call(rbind, lapply(split(by_twilight, list(by_twilight$pair_id, by_twilight$period), drop = TRUE), function(x) {
    data.frame(
      pair_id = x$pair_id[1],
      period = x$period[1],
      n = nrow(x),
      median_delta_error_km = stats::median(x$delta_error_km, na.rm = TRUE),
      mean_delta_error_km = mean(x$delta_error_km, na.rm = TRUE),
      q95_delta_error_km = unname(stats::quantile(x$delta_error_km, 0.95, na.rm = TRUE)),
      proportion_new_closer = mean(x$delta_error_km < 0, na.rm = TRUE)
    )
  }))
} else {
  data.frame()
}

comparison_path <- file.path(output_dir, "benchmark_old_new_comparison.csv")
twilight_path <- file.path(output_dir, "gps_old_new_by_twilight.csv")
period_path <- file.path(output_dir, "gps_old_new_by_period.csv")
notes_path <- file.path(output_dir, "diagnostic_notes.txt")

utils::write.csv(comparison, comparison_path, row.names = FALSE)
utils::write.csv(by_twilight, twilight_path, row.names = FALSE)
utils::write.csv(by_period, period_path, row.names = FALSE)

notes <- character()
add <- function(...) notes <<- c(notes, paste0(...))

add("FLightR local benchmark diagnostic notes")
add("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
add("")
if (nrow(comparison) == 0) {
  add("No valid old/new pairs were found. Pairing uses git_commit/git_branch as authoritative and excludes label/git inconsistencies.")
} else {
  add("Valid old/new pairs:")
  for (i in seq_len(nrow(comparison))) {
    add("- ", comparison$pair_id[i], " (nParticles=", comparison$nParticles[i],
        ", seed=", comparison$seed[i], ", subset=", comparison$subset[i], ")")
  }
  add("")
  add("Runtime:")
  for (i in seq_len(nrow(comparison))) {
    add("- ", comparison$pair_id[i], ": speedup_pf=", round(comparison$speedup_pf[i], 3),
        ", speedup_total=", round(comparison$speedup_total[i], 3))
  }
  add("")
  add("GPS accuracy deltas are new minus old; positive means old is closer.")
  for (i in seq_len(nrow(comparison))) {
    add("- ", comparison$pair_id[i], ": delta median=", round(comparison$delta_median_error_km[i], 3),
        " km, delta mean=", round(comparison$delta_mean_error_km[i], 3),
        " km, delta q75=", round(comparison$delta_q75_error_km[i], 3),
        " km, delta q95=", round(comparison$delta_q95_error_km[i], 3),
        " km, delta max=", round(comparison$delta_max_error_km[i], 3), " km")
  }
}

add("")
if (nrow(by_twilight) > 0) {
  add("Twilight-level delta distribution:")
  split_tw <- split(by_twilight, by_twilight$pair_id)
  for (nm in names(split_tw)) {
    x <- split_tw[[nm]]
    add("- ", nm, ": median delta=", round(stats::median(x$delta_error_km, na.rm = TRUE), 3),
        " km, q95 delta=", round(unname(stats::quantile(x$delta_error_km, 0.95, na.rm = TRUE)), 3),
        " km, proportion new closer=", round(mean(x$delta_error_km < 0, na.rm = TRUE), 3))
  }
  add("")
  add("Period-specific results:")
  if (is.data.frame(by_period) && nrow(by_period) > 0) {
    for (i in seq_len(nrow(by_period))) {
      add("- ", by_period$pair_id[i], " / ", by_period$period[i],
          ": n=", by_period$n[i],
          ", median delta=", round(by_period$median_delta_error_km[i], 3),
          " km, q95 delta=", round(by_period$q95_delta_error_km[i], 3),
          " km, proportion new closer=", round(by_period$proportion_new_closer[i], 3))
    }
  } else {
    add("- No period-specific paired rows were available.")
  }
} else {
  add("No twilight-level paired GPS comparisons were available.")
}

add("")
add("Particle-number/stochasticity diagnosis:")
if (nrow(valid_bench) == 0) {
  add("- No valid benchmark rows were available.")
} else if (any(valid_bench$nParticles == 10000, na.rm = TRUE)) {
  add("- 1e4 runs are available. Use them before drawing conclusions from 1e3 runs.")
} else if (any(valid_bench$nParticles == 1000, na.rm = TRUE)) {
  add("- 1e3 runs are available. If only one seed exists, evidence for accuracy differences is insufficient.")
} else {
  add("- No 1e3 or 1e4 benchmark rows were found.")
}

add("")
add("Likely cause of accuracy differences:")
add("- Dynamic transition encoding should not affect a 2591-cell grid unless decoding changed accidentally.")
add("- stop.point should not affect runs with known.last = FALSE.")
add("- threads=1 fix should not affect inference when threads=8 except through unchanged stochastic parallel behaviour.")
add("- The directional bearing fix is the main inference-changing candidate.")
add("- Low nParticles can produce particle-filter noise or particle impoverishment; one 1e3 run is not enough evidence.")

add("")
if (!any(valid_bench$nParticles == 10000, na.rm = TRUE)) {
  add("Recommended next experiment: run old/new halfyear with nParticles = 1e4, seed = 123.")
} else if (nrow(comparison) > 0 && any(abs(comparison$delta_median_error_km) > 10, na.rm = TRUE)) {
  add("Recommended next experiment: run old/new halfyear with nParticles = 1e3 across seeds 123, 456, 789 to assess stochasticity, then consider 1e5 if 1e4 still differs strongly.")
} else {
  add("Recommended next experiment: stay at 1e3 multi-seed if stochasticity remains unclear; do not jump directly to 1e6.")
}
add("Use 1e6 only after scripts are stable, 1e4/1e5 indicate the effect direction, expected runtime is known, and publication-quality validation is needed.")

writeLines(notes, notes_path)

print(comparison)
print(by_period)
message("Wrote: ", comparison_path)
message("Wrote: ", twilight_path)
message("Wrote: ", period_path)
message("Wrote: ", notes_path)
