# Local-only particle scaling benchmark for FLightR propagation backends.
# Reads private local validation data only through benchmark_double_tag_subset.R.
# Prints/saves aggregate timings only.

.libPaths(c('C:/Users/eldar/Documents/Codex/Rlibs-4.6', .libPaths()))
repo <- 'D:/GitHub/FLightR'
out_dir <- 'D:/GitHub/FLightR/data-raw/local_validation/outputs'
script <- 'D:/GitHub/FLightR/data-raw/local_validation/scripts/benchmark_double_tag_subset.R'
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(repo)

subset_start <- Sys.getenv('FLIGHTR_SCALING_START_DATE', unset = 'NA')
subset_end <- Sys.getenv('FLIGHTR_SCALING_END_DATE', unset = '2013-10-15')
seed <- Sys.getenv('FLIGHTR_SCALING_SEED', unset = '123')
legacy_1e6_pf_threshold_seconds <- as.numeric(Sys.getenv('FLIGHTR_LEGACY_1E6_PF_THRESHOLD_SECONDS', unset = '900'))

cases <- data.frame(
  case = c('cached_1e5', 'cached_1e6', 'legacy_1e5', 'legacy_1e6'),
  backend = c('cached', 'cached', 'legacy', 'legacy'),
  nParticles = c(100000L, 1000000L, 100000L, 1000000L),
  threads = 1L,
  stringsAsFactors = FALSE
)

run_one <- function(row) {
  run_label <- paste0('particle_scaling_short_movement_', row$backend, '_Np', row$nParticles, '_threads1_seed', seed)
  old <- Sys.getenv(c('FLIGHTR_RUN_LABEL','FLIGHTR_NPARTICLES','FLIGHTR_PRERUN_THREADS','FLIGHTR_PF_THREADS','FLIGHTR_THREADS','FLIGHTR_SEED','FLIGHTR_KNOWN_LAST','FLIGHTR_PROC_START_DATE','FLIGHTR_PROC_END_DATE','FLIGHTR_PROPAGATION_BACKEND','FLIGHTR_PROFILE_TOP_LEVEL','FLIGHTR_PROFILE_PHASES'))
  on.exit(do.call(Sys.setenv, as.list(old)), add = TRUE)
  Sys.setenv(
    FLIGHTR_RUN_LABEL = run_label,
    FLIGHTR_NPARTICLES = as.character(row$nParticles),
    FLIGHTR_PRERUN_THREADS = as.character(row$threads),
    FLIGHTR_PF_THREADS = "1",
    FLIGHTR_SEED = seed,
    FLIGHTR_KNOWN_LAST = 'true',
    FLIGHTR_PROC_START_DATE = subset_start,
    FLIGHTR_PROC_END_DATE = subset_end,
    FLIGHTR_PROPAGATION_BACKEND = row$backend,
    FLIGHTR_PROFILE_TOP_LEVEL = 'true',
    FLIGHTR_PROFILE_PHASES = 'false'
  )
  value <- source(script, local = new.env())$value
  b <- value$benchmark
  top <- value$Result$Results$top_level_profile
  pf_core <- top$elapsed_seconds[top$phase == 'pf.run.parallel.SO.resample'][1]
  run_pf <- top$elapsed_seconds[top$phase == 'total_run.particle.filter'][1]
  data.frame(
    case = row$case,
    backend = row$backend,
    nParticles = row$nParticles,
    threads = row$threads,
    subset_start = subset_start,
    subset_end = subset_end,
    includes_movement = TRUE,
    grid_size = b$grid_size,
    n_twilights = b$n_twilights,
    make_prerun_seconds = b$prerun_seconds,
    prerun_cache_enabled = if ("prerun_cache_enabled" %in% names(b)) b$prerun_cache_enabled else NA,
    prerun_cache_hit = if ("prerun_cache_hit" %in% names(b)) b$prerun_cache_hit else NA,
    prerun_cache_file = if ("prerun_cache_file" %in% names(b)) b$prerun_cache_file else NA_character_,
    prerun_cache_key = if ("prerun_cache_key" %in% names(b)) b$prerun_cache_key else NA_character_,
    prerun_build_seconds = if ("prerun_build_seconds" %in% names(b)) b$prerun_build_seconds else NA_real_,
    prerun_load_seconds = if ("prerun_load_seconds" %in% names(b)) b$prerun_load_seconds else NA_real_,
    prerun_save_seconds = if ("prerun_save_seconds" %in% names(b)) b$prerun_save_seconds else NA_real_,
    pf_core_seconds = pf_core,
    run_particle_filter_seconds = run_pf,
    total_benchmark_seconds = b$total_seconds,
    LL = b$LL,
    benchmark_path = value$benchmark_path,
    result_path = value$result_path,
    top_level_profile_path = value$top_level_profile_path,
    completed = TRUE,
    skipped = FALSE,
    skip_reason = NA_character_,
    stringsAsFactors = FALSE
  )
}

results <- list()
for (i in seq_len(nrow(cases))) {
  row <- cases[i, ]
  if (row$case == 'legacy_1e6') {
    legacy_1e5 <- do.call(rbind, results)[do.call(rbind, results)$case == 'legacy_1e5', , drop = FALSE]
    if (nrow(legacy_1e5) == 1) {
      projected <- legacy_1e5$pf_core_seconds * 10
      if (is.finite(projected) && projected > legacy_1e6_pf_threshold_seconds) {
        results[[length(results)+1]] <- data.frame(
          case = row$case, backend = row$backend, nParticles = row$nParticles, threads = row$threads,
          subset_start = subset_start, subset_end = subset_end, includes_movement = TRUE,
          grid_size = legacy_1e5$grid_size, n_twilights = legacy_1e5$n_twilights,
          make_prerun_seconds = NA_real_, prerun_cache_enabled = NA, prerun_cache_hit = NA, prerun_cache_file = NA_character_, prerun_cache_key = NA_character_, prerun_build_seconds = NA_real_, prerun_load_seconds = NA_real_, prerun_save_seconds = NA_real_, pf_core_seconds = NA_real_, run_particle_filter_seconds = NA_real_, total_benchmark_seconds = NA_real_, LL = NA_real_,
          benchmark_path = NA_character_, result_path = NA_character_, top_level_profile_path = NA_character_,
          completed = FALSE, skipped = TRUE,
          skip_reason = paste0('Projected PF seconds from legacy_1e5 exceeded threshold: ', round(projected, 2), ' > ', legacy_1e6_pf_threshold_seconds),
          stringsAsFactors = FALSE
        )
        next
      }
    }
  }
  results[[length(results)+1]] <- tryCatch(run_one(row), error = function(e) {
    data.frame(
      case = row$case, backend = row$backend, nParticles = row$nParticles, threads = row$threads,
      subset_start = subset_start, subset_end = subset_end, includes_movement = TRUE,
      grid_size = NA_integer_, n_twilights = NA_integer_, make_prerun_seconds = NA_real_, prerun_cache_enabled = NA, prerun_cache_hit = NA, prerun_cache_file = NA_character_, prerun_cache_key = NA_character_, prerun_build_seconds = NA_real_, prerun_load_seconds = NA_real_, prerun_save_seconds = NA_real_, pf_core_seconds = NA_real_, run_particle_filter_seconds = NA_real_, total_benchmark_seconds = NA_real_, LL = NA_real_,
      benchmark_path = NA_character_, result_path = NA_character_, top_level_profile_path = NA_character_, completed = FALSE, skipped = FALSE,
      skip_reason = conditionMessage(e), stringsAsFactors = FALSE
    )
  })
  utils::write.csv(do.call(rbind, results), file.path(out_dir, 'particle_scaling_short_movement_runtime.csv'), row.names = FALSE)
}
runtime <- do.call(rbind, results)
utils::write.csv(runtime, file.path(out_dir, 'particle_scaling_short_movement_runtime.csv'), row.names = FALSE)

correctness <- data.frame()
for (np in unique(runtime$nParticles[runtime$completed])) {
  pair <- runtime[runtime$nParticles == np & runtime$completed, ]
  if (all(c('cached','legacy') %in% pair$backend)) {
    cached <- pair[pair$backend == 'cached', ][1, ]
    legacy <- pair[pair$backend == 'legacy', ][1, ]
    correctness <- rbind(correctness, data.frame(
      nParticles = np,
      cached_LL = cached$LL,
      legacy_LL = legacy$LL,
      delta_LL = cached$LL - legacy$LL,
      same_threads = cached$threads == legacy$threads,
      note = 'Same seed/backend comparison; exact posterior comparison not loaded to avoid extra memory in scaling runner.',
      stringsAsFactors = FALSE
    ))
  }
}
utils::write.csv(correctness, file.path(out_dir, 'particle_scaling_short_movement_correctness.csv'), row.names = FALSE)

notes <- c(
  paste('Created:', format(Sys.time(), '%Y-%m-%d %H:%M:%S %Z')),
  paste('Subset start:', subset_start),
  paste('Subset end:', subset_end),
  'Subset was selected from an existing modelled full-year result as a movement-rich window; no raw GPS/GLS rows are printed.',
  '',
  'Runtime table:',
  paste(capture.output(print(runtime)), collapse = '\n'),
  '',
  'Correctness table:',
  paste(capture.output(print(correctness)), collapse = '\n'),
  '',
  'Raw private GPS/GLS rows were not printed. Nothing was staged or committed.'
)
writeLines(notes, file.path(out_dir, 'particle_scaling_short_movement_notes.txt'))
print(runtime)
print(correctness)
