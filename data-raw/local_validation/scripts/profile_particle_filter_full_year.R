# Local-only profiled full-year benchmark runner.
# This reads private validation data only through benchmark_double_tag_subset.R and
# writes aggregate diagnostic outputs under data-raw/local_validation/outputs/.

repo_root <- "D:/GitHub/FLightR"
output_root <- "D:/GitHub/FLightR/data-raw/local_validation/outputs"
benchmark_script <- "D:/GitHub/FLightR/data-raw/local_validation/scripts/benchmark_double_tag_subset.R"

run_label <- Sys.getenv("FLIGHTR_RUN_LABEL", unset = "pf_profile_full_year_1e4_threads4_seed123")
nParticles <- Sys.getenv("FLIGHTR_NPARTICLES", unset = "10000")
prerun_threads <- Sys.getenv("FLIGHTR_PRERUN_THREADS", unset = "4")
pf_threads <- Sys.getenv("FLIGHTR_PF_THREADS", unset = "1")
seed <- Sys.getenv("FLIGHTR_SEED", unset = "123")
known_last <- Sys.getenv("FLIGHTR_KNOWN_LAST", unset = "true")
rerun_existing <- tolower(Sys.getenv("FLIGHTR_PROFILE_RERUN", unset = "false")) %in% c("1", "true", "yes")

canonical_label <- "full_year_1e4_threads4_seed123"
rprof_file <- file.path(output_root, paste0("pf_rprof_", canonical_label, ".out"))
rprof_summary_file <- file.path(output_root, paste0("pf_rprof_summary_", canonical_label, ".txt"))
notes_file <- file.path(output_root, paste0("pf_profile_diagnostic_notes_", canonical_label, ".txt"))

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

if (!rerun_existing && file.exists(rprof_file)) {
  stop("Profile output already exists. Set FLIGHTR_PROFILE_RERUN=true to rerun: ", rprof_file)
}

old_wd <- getwd()
old_env <- Sys.getenv(c(
  "FLIGHTR_RUN_LABEL",
  "FLIGHTR_NPARTICLES",
  "FLIGHTR_PRERUN_THREADS","FLIGHTR_PF_THREADS","FLIGHTR_THREADS",
  "FLIGHTR_SEED",
  "FLIGHTR_PROC_START_DATE",
  "FLIGHTR_PROC_END_DATE",
  "FLIGHTR_KNOWN_LAST",
  "FLIGHTR_PROFILE_PHASES"
), unset = NA)

restore_env <- function() {
  for (nm in names(old_env)) {
    if (is.na(old_env[[nm]])) {
      Sys.unsetenv(nm)
    } else {
      do.call(Sys.setenv, as.list(stats::setNames(old_env[[nm]], nm)))
    }
  }
  setwd(old_wd)
}
on.exit(restore_env(), add = TRUE)

setwd(repo_root)
Sys.setenv(
  FLIGHTR_RUN_LABEL = run_label,
  FLIGHTR_NPARTICLES = nParticles,
  FLIGHTR_PRERUN_THREADS = prerun_threads,
    FLIGHTR_PF_THREADS = pf_threads,
  FLIGHTR_SEED = seed,
  FLIGHTR_PROC_START_DATE = "NA",
  FLIGHTR_PROC_END_DATE = "NA",
  FLIGHTR_KNOWN_LAST = known_last,
  FLIGHTR_PROFILE_PHASES = "true"
)

message("Starting profiled benchmark: ", run_label)
message("Outputs: ", output_root)
message("Rprof: ", rprof_file)

Rprof(rprof_file, interval = 0.02, memory.profiling = TRUE)
profile_result <- tryCatch(
  source(benchmark_script),
  finally = Rprof(NULL)
)

summary_total <- summaryRprof(rprof_file, memory = "both")
summary_self <- summaryRprof(rprof_file, memory = "both")

capture.output({
  cat("Rprof summary for ", run_label, "\n", sep = "")
  cat("\nBy total time:\n")
  print(summary_total$by.total)
  cat("\nBy self time:\n")
  print(summary_self$by.self)
  cat("\nMemory by total, if available:\n")
  print(summary_total$by.total)
}, file = rprof_summary_file)

returned <- if (is.list(profile_result) && "value" %in% names(profile_result)) profile_result$value else profile_result

phase_src <- if (is.list(returned)) returned$phase_profile_path else NA_character_
top_src <- if (is.list(returned)) returned$top_level_profile_path else NA_character_
benchmark_path <- if (is.list(returned)) returned$benchmark_path else NA_character_
result_path <- if (is.list(returned)) returned$result_path else NA_character_

phase_dst <- file.path(output_root, paste0("pf_phase_profile_", canonical_label, ".csv"))
top_dst <- file.path(output_root, paste0("pf_top_level_timing_", canonical_label, ".csv"))
if (!is.na(phase_src) && file.exists(phase_src) && !identical(normalizePath(phase_src, winslash = "/", mustWork = FALSE), normalizePath(phase_dst, winslash = "/", mustWork = FALSE))) {
  file.copy(phase_src, phase_dst, overwrite = TRUE)
}
if (!is.na(top_src) && file.exists(top_src) && !identical(normalizePath(top_src, winslash = "/", mustWork = FALSE), normalizePath(top_dst, winslash = "/", mustWork = FALSE))) {
  file.copy(top_src, top_dst, overwrite = TRUE)
}

notes <- c(
  paste0("Profiled benchmark notes"),
  paste0("timestamp=", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  paste0("run_label=", run_label),
  paste0("repo_root=", repo_root),
  paste0("nParticles=", nParticles),
  paste0("threads=", threads),
  paste0("seed=", seed),
  paste0("known.last=", known_last),
  "proc_start_date=NA",
  "proc_end_date=NA",
  "profile.phases=true",
  paste0("benchmark_path=", benchmark_path),
  paste0("result_path=", result_path),
  paste0("phase_profile_path=", if (file.exists(phase_dst)) phase_dst else phase_src),
  paste0("top_level_timing_path=", if (file.exists(top_dst)) top_dst else top_src),
  paste0("rprof_file=", rprof_file),
  paste0("rprof_summary_file=", rprof_summary_file),
  "Private GPS/GLS rows were not printed by this script."
)
writeLines(notes, notes_file)

message("Profiled benchmark finished.")
message("Benchmark: ", benchmark_path)
message("Phase profile: ", if (file.exists(phase_dst)) phase_dst else phase_src)
message("Top-level timing: ", if (file.exists(top_dst)) top_dst else top_src)
message("Rprof summary: ", rprof_summary_file)

invisible(list(
  benchmark_path = benchmark_path,
  result_path = result_path,
  phase_profile_path = if (file.exists(phase_dst)) phase_dst else phase_src,
  top_level_profile_path = if (file.exists(top_dst)) top_dst else top_src,
  rprof_file = rprof_file,
  rprof_summary_file = rprof_summary_file,
  notes_file = notes_file
))
