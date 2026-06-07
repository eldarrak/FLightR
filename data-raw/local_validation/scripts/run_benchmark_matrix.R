versions <- c("old_master", "new_branch")
repo_path <- c(
  old_master = "D:/GitHub/FLightR_old_master",
  new_branch = "D:/GitHub/FLightR"
)
nParticles <- as.integer(strsplit(Sys.getenv("FLIGHTR_MATRIX_NPARTICLES", "1000"), ",")[[1]])
seeds <- as.integer(strsplit(Sys.getenv("FLIGHTR_MATRIX_SEEDS", "123,456,789"), ",")[[1]])
total_cores <- as.integer(Sys.getenv("FLIGHTR_MATRIX_TOTAL_CORES", "8"))
parallel_jobs <- as.integer(Sys.getenv("FLIGHTR_MATRIX_PARALLEL_JOBS", "1"))
parallel_jobs <- max(1L, parallel_jobs)
threads <- max(1L, floor(total_cores / parallel_jobs))
rerun_existing <- FALSE
actually_run <- tolower(Sys.getenv("FLIGHTR_MATRIX_ACTUALLY_RUN", "false")) %in% c("1", "true", "yes")
subset_label <- Sys.getenv("FLIGHTR_MATRIX_SUBSET_LABEL", "halfyear")
known_last <- Sys.getenv("FLIGHTR_MATRIX_KNOWN_LAST", "false")
proc_start_date <- Sys.getenv("FLIGHTR_MATRIX_PROC_START_DATE", "NA")
proc_end_date <- Sys.getenv("FLIGHTR_MATRIX_PROC_END_DATE", "2014-01-15")

output_dir <- "D:/GitHub/FLightR/data-raw/local_validation/outputs"
benchmark_script <- "D:/GitHub/FLightR/data-raw/local_validation/scripts/benchmark_double_tag_subset.R"

existing_benchmarks <- list.files(output_dir, pattern = "^benchmark_.*\\.csv$", full.names = FALSE)

estimate_seconds <- function(version, np) {
  combined_path <- file.path(output_dir, "benchmark_combined.csv")
  if (!file.exists(combined_path)) return(NA_real_)
  b <- utils::read.csv(combined_path, stringsAsFactors = FALSE)
  b <- b[b$version == version & b$nParticles == np, ]
  if (nrow(b) == 0) return(NA_real_)
  mean(b$total_seconds, na.rm = TRUE)
}

plan <- expand.grid(
  version = versions,
  nParticles = nParticles,
  seed = seeds,
  stringsAsFactors = FALSE
)
plan$threads <- threads
plan$repo_path <- repo_path[plan$version]
plan$run_label <- paste(plan$version, subset_label, paste0("1e", nchar(plan$nParticles) - 1),
                        paste0("seed", plan$seed), sep = "_")
plan$estimated_total_seconds <- mapply(estimate_seconds, plan$version, plan$nParticles)

if (!actually_run) {
  message("Dry run only. Set FLIGHTR_MATRIX_ACTUALLY_RUN=true to launch benchmarks.")
  message("Parallel jobs: ", parallel_jobs, "; threads per run: ", threads,
          "; total core budget: ", total_cores)
  print(plan)
  quit(save = "no", status = 0)
}

if (!file.exists(benchmark_script)) {
  stop("Benchmark script not found: ", benchmark_script, call. = FALSE)
}

run_one_benchmark <- function(row) {
  result_pattern <- paste0("Result_", row$run_label, ".*Np", row$nParticles, ".*\\.rds$")
  benchmark_pattern <- paste0("benchmark_", row$run_label, ".*Np", row$nParticles, ".*\\.csv$")
  has_result <- length(list.files(output_dir, pattern = result_pattern)) > 0
  has_benchmark <- length(list.files(output_dir, pattern = benchmark_pattern)) > 0
  if (!rerun_existing && has_result && has_benchmark) {
    message("Skipping existing run: ", row$run_label)
    return(data.frame(run_label = row$run_label, status = "skipped"))
  }
  if (!dir.exists(row$repo_path)) {
    stop("Repository path not found: ", row$repo_path, call. = FALSE)
  }

  log_file <- file.path(output_dir, paste0("benchmark_matrix_", row$run_label, ".log"))
  err_file <- file.path(output_dir, paste0("benchmark_matrix_", row$run_label, ".err.log"))
  out_con <- file(log_file, open = "at")
  err_con <- file(err_file, open = "at")
  sink(out_con, type = "output")
  sink(err_con, type = "message")
  on.exit({
    sink(type = "message")
    sink(type = "output")
    close(err_con)
    close(out_con)
  }, add = TRUE)

  message("Running: ", row$run_label)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(row$repo_path)
  Sys.setenv(
    FLIGHTR_RUN_LABEL = row$run_label,
    FLIGHTR_NPARTICLES = as.character(row$nParticles),
    FLIGHTR_THREADS = as.character(row$threads),
    FLIGHTR_SEED = as.character(row$seed),
    FLIGHTR_KNOWN_LAST = known_last,
    FLIGHTR_PROC_START_DATE = proc_start_date,
    FLIGHTR_PROC_END_DATE = proc_end_date
  )
  source(benchmark_script)
  data.frame(run_label = row$run_label, status = "completed")
}

plan_rows <- split(plan, seq_len(nrow(plan)))
if (parallel_jobs <= 1L) {
  results <- do.call(rbind, lapply(plan_rows, run_one_benchmark))
} else {
  cl <- parallel::makeCluster(parallel_jobs)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterExport(
    cl,
    c(
      "benchmark_script",
      "output_dir",
      "rerun_existing",
      "run_one_benchmark",
      "known_last",
      "proc_start_date",
      "proc_end_date"
    ),
    envir = environment()
  )
  results <- do.call(rbind, parallel::parLapply(cl, plan_rows, run_one_benchmark))
}

print(results)
