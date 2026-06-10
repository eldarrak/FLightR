# Local-only paired benchmark runner for FLightR private validation data.
# Dry-run by default. Do not commit private outputs.

# ----------------------------
# Editable settings
# ----------------------------

nParticles <- as.integer(Sys.getenv("FLIGHTR_PAIR_NPARTICLES", "100000"))
threads_per_run <- as.integer(Sys.getenv("FLIGHTR_PAIR_THREADS_PER_RUN", "4"))
seed <- as.integer(Sys.getenv("FLIGHTR_PAIR_SEED", "123"))
pair_id <- Sys.getenv("FLIGHTR_PAIR_ID", "pair1")

run_baseline <- tolower(Sys.getenv("FLIGHTR_PAIR_RUN_BASELINE", "true")) %in% c("1", "true", "yes")
run_s2 <- tolower(Sys.getenv("FLIGHTR_PAIR_RUN_S2", "true")) %in% c("1", "true", "yes")
actually_run <- tolower(Sys.getenv("FLIGHTR_PAIR_ACTUALLY_RUN", "false")) %in% c("1", "true", "yes")
rerun_existing <- tolower(Sys.getenv("FLIGHTR_PAIR_RERUN_EXISTING", "false")) %in% c("1", "true", "yes")

repo_paths <- c(
  baseline = "D:/GitHub/FLightR_fixed_baseline",
  s2 = "D:/GitHub/FLightR"
)

output_dir <- "D:/GitHub/FLightR/data-raw/local_validation/outputs"
benchmark_script <- "D:/GitHub/FLightR/data-raw/local_validation/scripts/benchmark_double_tag_subset.R"

# ----------------------------
# Helpers
# ----------------------------

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

particle_label <- function(n) {
  if (n <= 0 || is.na(n)) stop("nParticles must be a positive integer.", call. = FALSE)
  pow <- log10(n)
  if (abs(pow - round(pow)) < .Machine$double.eps^0.5) {
    paste0("1e", round(pow))
  } else {
    as.character(n)
  }
}

label_for <- function(version) {
  paste(version, "halfyear", particle_label(nParticles), pair_id, paste0("seed", seed), sep = "_")
}

find_outputs <- function(run_label) {
  result_pattern <- paste0("^Result_", run_label, ".*Np", nParticles, ".*\\.rds$")
  benchmark_pattern <- paste0("^benchmark_", run_label, ".*Np", nParticles, ".*\\.csv$")
  list(
    result = list.files(output_dir, pattern = result_pattern, full.names = TRUE),
    benchmark = list.files(output_dir, pattern = benchmark_pattern, full.names = TRUE)
  )
}

run_exists <- function(run_label) {
  found <- find_outputs(run_label)
  length(found$result) > 0 && length(found$benchmark) > 0
}

run_one <- function(version, repo_path) {
  run_label <- label_for(version)
  exists <- run_exists(run_label)
  status <- if (exists && !rerun_existing) "skip_existing" else if (actually_run) "run" else "dry_run"

  message(version, ": ", status, " label=", run_label,
          " nParticles=", nParticles, " threads=", threads_per_run, " seed=", seed)

  if (!dir.exists(repo_path)) {
    if (!actually_run) {
      warning("Repository path not found for ", version, ": ", repo_path, call. = FALSE)
      return(data.frame(
        version = version,
        run_label = run_label,
        repo_path = repo_path,
        status = "dry_run_missing_repo",
        stringsAsFactors = FALSE
      ))
    }
    stop("Repository path not found for ", version, ": ", repo_path, call. = FALSE)
  }
  if (!file.exists(benchmark_script)) {
    stop("Benchmark script not found: ", benchmark_script, call. = FALSE)
  }

  if (!actually_run || (exists && !rerun_existing)) {
    return(data.frame(
      version = version,
      run_label = run_label,
      repo_path = repo_path,
      status = status,
      stringsAsFactors = FALSE
    ))
  }

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(repo_path)

  Sys.setenv(
    FLIGHTR_RUN_LABEL = run_label,
    FLIGHTR_NPARTICLES = as.character(nParticles),
    FLIGHTR_PRERUN_THREADS = as.character(threads_per_run),
    FLIGHTR_PF_THREADS = Sys.getenv("FLIGHTR_PAIR_PF_THREADS_PER_RUN", "1"),
    FLIGHTR_SEED = as.character(seed)
  )

  source(benchmark_script)

  data.frame(
    version = version,
    run_label = run_label,
    repo_path = repo_path,
    status = "completed",
    stringsAsFactors = FALSE
  )
}

message("FLightR paired benchmark runner")
message("actually_run=", actually_run, " rerun_existing=", rerun_existing)
message("outputs: ", output_dir)

plan <- list()
if (run_baseline) plan$baseline <- repo_paths["baseline"]
if (run_s2) plan$s2 <- repo_paths["s2"]

if (length(plan) == 0) {
  stop("Nothing selected: set run_baseline and/or run_s2 to TRUE.", call. = FALSE)
}

results <- do.call(rbind, Map(run_one, names(plan), unname(plan)))
print(results)

if (!actually_run) {
  message("Dry run only. Set actually_run=TRUE near the top or FLIGHTR_PAIR_ACTUALLY_RUN=true to run.")
}
