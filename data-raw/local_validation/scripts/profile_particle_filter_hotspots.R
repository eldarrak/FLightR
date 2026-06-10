# Local-only particle-filter hotspot profiling for FLightR.
# Reads private validation inputs only through benchmark_double_tag_subset.R.
# Writes aggregate profiling outputs under data-raw/local_validation/outputs/.

repo_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
output_root <- "D:/GitHub/FLightR/data-raw/local_validation/outputs"
benchmark_script <- "D:/GitHub/FLightR/data-raw/local_validation/scripts/benchmark_double_tag_subset.R"

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_label <- paste0("hotspot_profile_halfyear_1e3_threads1_seed123_", timestamp)
rprof_path <- file.path(output_root, paste0("particle_filter_hotspot_profile_", timestamp, ".out"))
summary_path <- file.path(output_root, "particle_filter_hotspot_summary.csv")
notes_path <- file.path(output_root, "particle_filter_hotspot_notes.txt")

Sys.setenv(
  FLIGHTR_RUN_LABEL = run_label,
  FLIGHTR_NPARTICLES = "1000",
  FLIGHTR_PRERUN_THREADS = "4",
  FLIGHTR_PF_THREADS = "1",
  FLIGHTR_SEED = "123",
  FLIGHTR_PROC_START_DATE = "NA",
  FLIGHTR_PROC_END_DATE = "2014-01-15",
  FLIGHTR_KNOWN_LAST = "false"
)

if (!file.exists(benchmark_script)) {
  stop("Benchmark script not found: ", benchmark_script, call. = FALSE)
}

utils::Rprof(rprof_path, interval = 0.01, memory.profiling = FALSE)
profile_status <- "completed"
profile_error <- NULL
tryCatch(
  source(benchmark_script),
  error = function(e) {
    profile_status <<- "error"
    profile_error <<- conditionMessage(e)
  },
  finally = utils::Rprof(NULL)
)

profile <- utils::summaryRprof(rprof_path)
flat <- as.data.frame(profile$by.self)
flat$function_name <- rownames(flat)
rownames(flat) <- NULL

targets <- c(
  "pf.run.parallel.SO.resample", "generate.points.dirs",
  "sf::st_distance", "st_distance", "geosphere::bearing", "bearing",
  "get.transition.rle", "sort.int", "sort", "rle", "sample", "sample.int",
  "rowProds", "cbind", "inverse.rle", "get.coordinates.PF",
  "get_ZI_distances", "dist.fun", "estimate.movement.parameters"
)

flat$target_group <- NA_character_
for (target in targets) {
  hit <- if (target %in% c("sample", "sort", "rle", "bearing")) {
    flat$function_name %in% paste0("\"", target, "\"")
  } else {
    grepl(target, flat$function_name, fixed = TRUE)
  }
  flat$target_group[hit & is.na(flat$target_group)] <- target
}

summary <- flat[!is.na(flat$target_group), c("target_group", "function_name", "self.time", "self.pct", "total.time", "total.pct")]
summary <- summary[order(summary$self.time, decreasing = TRUE), ]

utils::write.csv(summary, summary_path, row.names = FALSE)

notes <- c(
  "Particle-filter hotspot profile",
  paste0("Date/time: ", Sys.time()),
  paste0("Repo: ", repo_root),
  paste0("Run label: ", run_label),
  "Setup: half-year subset, nParticles=1000, threads=1, seed=123, known.last=FALSE, check.outliers=FALSE via benchmark script.",
  paste0("Status: ", profile_status),
  if (!is.null(profile_error)) paste0("Error: ", profile_error) else NULL,
  paste0("Raw Rprof output: ", rprof_path),
  paste0("Aggregate summary CSV: ", summary_path),
  "Only aggregate profiling summaries are written; raw GPS/GLS rows are not printed."
)
writeLines(notes, notes_path)

if (profile_status == "error") {
  stop(profile_error, call. = FALSE)
}

print(utils::head(summary, 25))
message("Profile summary saved to: ", summary_path)
message("Profile notes saved to: ", notes_path)
